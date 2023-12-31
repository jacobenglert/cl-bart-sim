#!/bin/bash

array_job0=$(sbatch --parsable HPC/friedman-sim.slurm 0)
array_job1=$(sbatch --parsable HPC/friedman-sim.slurm 1)

sbatch --depend=afterany:$array_job1:$array_job0 HPC/friedman-sim-combine.slurm
