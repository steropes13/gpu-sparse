#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=test
#SBATCH --output=test-%j.out
#SBATCH --error=test-%j.err

module load CUDA/11.8.0

make clean

make

rm res.txt 

## 1st parameter -> 0 = COO, 1 = CSR, it writes in a "res.txt" file in the current directory
## 1nd parameter -> sliceSize (optionnal)
./bin/spmvCusparse ../del-1/GPU_solution/mtx_matrix/cage15/cage15.mtx 2 1



