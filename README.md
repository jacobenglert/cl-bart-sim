
# cl-bart-sim

<!-- badges: start -->
<!-- badges: end -->

This repository contains the R code for two simulation studies evaluating the performance of [CL-BART](https://github.com/jacobenglert/clbart). These simulations were originally ran on the [RSPH HPC cluster](https://scholarblogs.emory.edu/rsph-hpc/), so shell scripts are also included to execute the analysis on such a system that utilizes the SLURM job scheduler. Descriptions of the simulations are available at https://arxiv.org/abs/2311.12016.

In brief:

- R/`*-sim-example.R` provides an example script of what a single simulation run looks like.
- R/`set-params.R` sets the parameters for each simulation (stored in Params folder)
- R/`*-sim-production.R` contains the streamlined production version and is called by the shell scripts in the HPC folder. Simulation results are stored in dated subfolders within the Results folder.
- R/`*-sim-analysis.R` processes results to create interim datasets for tables and figures.

