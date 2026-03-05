import os
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import torch
from metvae.model import MetVAE
from metvae.sim import sim_data

# Function to run tasks in parallel and gather results
def _choose_rho(d: int) -> float:
    """Map feature dimension to SEC penalty."""
    if d == 50:
        return 0.9
    if d == 200:
        return 1.5
    if d == 500:
        return 2.2
    # Sensible default if other sizes sneak in
    return 1.5

def run_simulation(n: int, d: int, zero_prop: float, seed: int):
    
    try:
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
    except Exception:
        pass
    
    cor_pairs = int(0.2 * d)
    mu = list(range(10, 15))
    da_prop = 0.1
    rho = _choose_rho(d)

    # Simulate data
    np.random.seed(seed)
    
    sim = sim_data(
        n=n, d=d, cor_pairs=cor_pairs, mu=mu, da_prop=da_prop
    )
    y = sim['y']
    true_cor = sim['cor_matrix']
    
    # Apply log transformation and add biases
    log_y = np.log(y)
    log_sample_bias = np.log(np.random.uniform(1e-3, 1e-1, size=n))
    log_feature_bias = np.log(np.random.uniform(1e-1, 1, size=d))
    log_data = log_y + log_sample_bias[:, np.newaxis]  # Adding sample bias
    log_data = log_data + log_feature_bias.reshape(1, d)  # Adding feature bias
    data = np.exp(log_data)

    # Calculate thresholds and apply zeros
    thresholds = np.quantile(data, zero_prop, axis=0)
    data_miss = np.where(data<thresholds, 0, data)
    data_miss = pd.DataFrame(
        data_miss,
        index=y.index,
        columns=y.columns
    )

    # Run the MetVAE model
    max_epochs=1000
    learning_rate=1e-2

    try:
        model = MetVAE(
            data=data_miss,
            features_as_rows=False,
            meta=None,
            continuous_covariate_keys=None,
            categorical_covariate_keys=None,
            latent_dim=min(n, d),
            seed=0
        )
        
        model.train(
            batch_size=100,
            num_workers=0,
            max_epochs=max_epochs,
            learning_rate=learning_rate,
            log_every_n_steps=1
        )

        # SEC sparsification
        model.get_corr(num_sim=100, workers=1, seed=0)
        results_metvae = model.sparse_by_sec(rho=rho)
        est_cor = results_metvae['sparse_estimate'].values
    
        # Calculate summary statistics
        true_idx = true_cor[np.tril_indices_from(true_cor, k=-1)] != 0
        est_idx = est_cor[np.tril_indices_from(est_cor, k=-1)] != 0
        tpr = np.sum(est_idx & true_idx) / np.sum(true_idx)
        fpr = np.sum(est_idx & ~true_idx) / np.sum(~true_idx)
        fdr = np.sum(est_idx & ~true_idx) / np.sum(est_idx)
    except Exception as e:
        print(f"An error occurred: {e}")
        tpr, fpr, fdr = np.nan, np.nan, np.nan
    
    return d, zero_prop, seed, rho, tpr, fpr, fdr

def simulation_wrapper(params):
    """Wrapper function to unpack parameters and call run_simulation."""
    return run_simulation(*params)

def run_simulations_parallel(simparams: pd.DataFrame, max_workers: int) -> pd.DataFrame:
    # Convert DataFrame rows to list of tuples, each tuple representing parameters for one simulation
    params_list = [(int(n), int(d), float(zp), int(s)) 
                   for n, d, zp, s in simparams.to_numpy()]
    
    # ---- Prevent thread oversubscription in each process ----
    # Set in parent so children inherit
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    # Bound workers sensibly
    cpu = os.cpu_count() or 1
    W = min(max_workers, cpu, len(params_list))

    # Use a larger chunksize to reduce scheduling overhead
    chunksize = max(1, len(params_list) // (4 * W))
    with ProcessPoolExecutor(max_workers=W) as ex:
        results = list(ex.map(simulation_wrapper, params_list, chunksize=chunksize))

    cols = ['d', 'zero_prop', 'seed', 'rho', 'TPR', 'FPR', 'FDR']
    return pd.DataFrame(results, columns=cols)

def main():
    n = 100
    d_values = [50, 200, 500]
    zero_prop = np.arange(0.0, 0.4, 0.1)
    iter_num = 100
    seeds = np.arange(iter_num, dtype=int)
    max_workers = 14

    # Build full grid including d
    simparams_full = pd.DataFrame(
        [(n, d, zp, s) 
         for d in d_values
         for zp in zero_prop 
         for s in seeds],
        columns=["n", "d", "zero_prop", "seed"]
    )

    print(f"[INFO] total tasks: {len(simparams_full)}", flush=True)

    res = run_simulations_parallel(simparams_full, max_workers=max_workers)

    print(f"[INFO] ALL DONE ({len(res)} rows).", flush=True)

    res.to_csv("sim_metvae_unconfound.csv", index=False)

if __name__ == '__main__':
    main()



