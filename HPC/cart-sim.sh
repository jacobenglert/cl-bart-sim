#!/bin/bash

array_job0=$(sbatch --parsable HPC/cart-sim.slurm 0)
array_job1=$(sbatch --parsable HPC/cart-sim.slurm 1)
array_job2=$(sbatch --parsable HPC/cart-sim.slurm 2)
array_job3=$(sbatch --parsable HPC/cart-sim.slurm 3)
array_job4=$(sbatch --parsable HPC/cart-sim.slurm 4)
array_job5=$(sbatch --parsable HPC/cart-sim.slurm 5)

sbatch --depend=afterany:$array_job5:$array_job4:$array_job3:$array_job2:$array_job1:$array_job0 HPC/cart-combine.slurm
