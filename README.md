# MetVAE-Code-Archive
This is the repository archiving data and scripts for reproducing results presented in the paper "Variational Autoencoders for Metabolomics: Data Imputation, Deconfounding, and Correlation Discovery".

**For the corresponding Python package, refer to [MetVAE](https://github.com/FrederickHuangLin/MetVAE-PyPI) repository.**

## Notebooks

Run the notebooks from the `code/` directory; all paths are relative to it.

| Notebook | Contents |
| --- | --- |
| `code/01_quickstart.ipynb` | Minimal end-to-end example of the MetVAE workflow on simulated data. |
| `code/02_sim_study.ipynb` | Simulation studies and the figures comparing MetVAE with the benchmark methods. |
| `code/03_hcc.ipynb` | HCC worked example: preprocessing, the reference model, the sparse correlation network, and the exported GraphML. |
| `code/04_imputation_comparison.ipynb` | Comparison of the ways of handling zeros in a correlation analysis, under a latent-factor model, producing Figure S3. |
| `code/05_hcc_diagnostics.ipynb` | HCC diagnostics: reconstruction, sensitivity to the training settings, diet adjustment, and cutoff and subsampling stability, producing Figures S4 to S6. |

`code/simulation_code/` holds the standalone scripts (Python and R, with Slurm job files) that run the benchmark methods behind the simulation studies and write the `sim_*.csv` files read by `code/02_sim_study.ipynb`.

## Environment

`environment.yml` builds the conda environment, and `requirements-lock.txt` pins the exact package versions, including metvae 1.1.0.

## Checkpoint

`results/intermediate_results/hcc_checkpoint.pth` holds the trained weights of the reference HCC model, so `03_hcc.ipynb` and `05_hcc_diagnostics.ipynb` reproduce the reference network without retraining.
