from concurrent.futures import ProcessPoolExecutor
import random
import numpy as np
import pandas as pd
import torch
from metvae.model import MetVAE
from metvae.sim import sim_data

# Function to run tasks in parallel and gather results
def run_simulation(n, d, zero_prop, seed):
    cor_pairs = int(0.2 * d)
    mu = list(range(10, 15))
    da_prop = 0.1

    # Simulate data
    np.random.seed(seed)
    sim = sim_data(n=n, d=d, cor_pairs=cor_pairs, mu=mu, da_prop=da_prop)
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
    torch.manual_seed(123)
    np.random.seed(123)
    
    max_epochs=1000
    learning_rate=1e-2

    try:
        model = MetVAE(data=data_miss,
                       features_as_rows=False,
                       meta=None,
                       continuous_covariate_keys=None,
                       categorical_covariate_keys=None,
                       latent_dim=min(n, d))
        
        model.train(batch_size=100,
                    num_workers=0,
                    max_epochs=max_epochs,
                    learning_rate=learning_rate,
                    log_every_n_steps=1)
        
        model.get_corr(num_sim=1000)
        random.seed(123)
        results_metvae = model.sparse_by_thresholding(th_len=100, n_cv=5, soft=False, n_jobs=1)
        est_cor = results_metvae['sparse_estimate']
    
        # Calculate summary statistics
        true_idx = true_cor[np.tril_indices_from(true_cor, k=-1)] != 0
        est_idx = est_cor[np.tril_indices_from(est_cor, k=-1)] != 0
        tpr = np.sum(est_idx & true_idx) / np.sum(true_idx)
        fpr = np.sum(est_idx & ~true_idx) / np.sum(~true_idx)
        fdr = np.sum(est_idx & ~true_idx) / np.sum(est_idx)
    except Exception as e:
        print(f"An error occurred: {e}")
        tpr, fpr, fdr = np.nan, np.nan, np.nan
    
    return tpr, fpr, fdr

def simulation_wrapper(params):
    """Wrapper function to unpack parameters and call run_simulation."""
    return run_simulation(*params)

def run_simulations_parallel(simparams, max_workers):
    # Convert DataFrame rows to list of tuples, each tuple representing parameters for one simulation
    params_list = [(int(n), int(d), float(zp), int(s)) for n, d, zp, s in simparams.to_numpy()]
    
    # Use ProcessPoolExecutor to run simulations in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(simulation_wrapper, params_list))
    
    # Create a new DataFrame for results
    res_sim = pd.DataFrame(results, columns=['TPR', 'FPR', 'FDR'])
    
    return res_sim

def main():
    n = 100
    d = [50, 200, 500]
    zero_prop = np.arange(0, 0.4, 0.1)
    iter_num = 100
    seed = np.arange(iter_num)
    max_workers = 30

    simparams = pd.DataFrame([(n, di, zp, s) for di in d for zp in zero_prop for s in seed],
                             columns=["n", "d", "zero_prop", "seed"])
    res_sim = run_simulations_parallel(simparams=simparams, max_workers=max_workers)
    res_sim.to_csv('sim_metvae_unconfound.csv', index=False)

if __name__ == '__main__':
    main()



