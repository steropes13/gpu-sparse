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

VAL=$(hostname) 
if [[ $VAL == "edu01" ]] ; then
		echo -e "\e[32m ============== NAME OF THE HOST (edu01/02) : $VAL ============= \e[0m"
	else  
		echo -e "\e[31m ============== NAME OF THE HOST (edu01/02) : $VAL ============= \e[0m"
fi

lscpu ## for each job we list the cpu used



 ./bin/spmv mtx_matrix/ASIC_680ks/ASIC_680ks.mtx 25 32  


