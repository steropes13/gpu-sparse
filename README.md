# Sparse-Matrix dense-Vector multiplication (SpMV) 
in this repository you will find files organized into diferent parts : 
  - """GPU_solution /""" : all the files necessary for the GPU computation and the final deliverable (time, matrix market reading, CPU comparison, cuSparese...)   
      - """src/""" : all the sources in """.c/.cu"""
          - """spmvCPU.c""" : CPU solution (_you can uncomment the printfs to see the results if necessary_) with time measurement in timeofday, you can e
          - """cuSparseComputation.cu""" : Cuda file that uses the cuSPARSE API provided by NVIDIA in order to compare time solutions between mine and them. It uses the cudaEvent and it provided also a comparison function (reltaive error) 
          - """mmio.c""" : some samples of matrix market reading provided by suiteSparse, modified for symmetric matrices.
          - """my_time_lib.c""" : functions that allows time reading using timeofday to compare with Cuda Event
      - """include/""" : headers files in """.h/cuh""" ncessary to the main file 
      - """mtx_matrix/""" : some samples of matrix in """.mtx""" format, due to the size of one of them (cage15,ASIC...) I do not provide each one in this folder
      - """bin/""" : the place we find the binary file called """spmv""", it has to be used with **<path to matrix in .mtx>,** **sliceSize** number and **seed** (_optional_) if you want to generate a random vector .
      - """obj/""" : all the files .o used for the compilation with **nvcc**.
      - """spmv.cu""" : the main file that contains the GPU solutions (COO,CSR/CSR WARP, Ellpack, Sliced-Ellpack) it gives you each kernel in res_<kernel_type>.txt files and also the cuSPARSE results to compare with. 
      - """sbatch.sh""" : necessary for the job loading for Slurm, you can rename it and it will give you some information such as the cpu used, make clean,  load cuda 11.8.0 module, node (edu01/02), GPU properties. Finally it will start the """./bin/spmv """ with a matrix. 
      - """makefile""" : used for the compilation via nvcc.
      - """report/""" : all the files used for the report writing in LaTeX (""".tex, .pdf""")
   
## Bonus
As we can't verify with cuSPARSRE the time execution and as it is not mandatory for this deliverable I immplemented 2 other SpMV methods : ELLpack (ELL) and Sliced-ELLpack (SELL). It is faster than COO,CSR (WARP) in certain cases... 
