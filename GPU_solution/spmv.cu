#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include "include/mmio.h"
#include "include/my_time_lib.h"
#include "include/spmvCPU.h"

//CUDA HEADERS (Time computation, cusparse comparison)
#include "include/cuSparseComputation.cuh"




#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float


#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUDA(func)                                                       \
{                                                                              \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})



// Print device properties
void printDevProp(cudaDeviceProp devProp)
{
    printf("Major revision number:         %d\n",  devProp.major);
    printf("Minor revision number:         %d\n",  devProp.minor);
    printf("Name:                          %s\n",  devProp.name);
    printf("  Memory Clock rate:           %.0f Mhz\n", devProp.memoryClockRate * 1e-3f);

    printf("  Memory Bus Width:            %d bit\n",devProp.memoryBusWidth);

    printf("  Peak Memory Bandwidth:       %7.3f GB/s\n",2.0*devProp.memoryClockRate*(devProp.memoryBusWidth/8)/1.0e6);

    printf("  Multiprocessors:             %3d\n",devProp.multiProcessorCount);
    printf("  Maximum number of threads per multiprocessor:  %d\n",devProp.maxThreadsPerMultiProcessor);
    printf("  Maximum number of threads per block:           %d\n",devProp.maxThreadsPerBlock);
    printf("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
           devProp.maxThreadsDim[0], devProp.maxThreadsDim[1],devProp.maxThreadsDim[2]);
    printf("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
           devProp.maxGridSize[0], devProp.maxGridSize[1],devProp.maxGridSize[2]);
    printf("  Total amount of shared memory per block:       %zu bytes\n", devProp.sharedMemPerBlock);
    return;

}

__global__ 
void computeSpmvCOOGPU(float * res, int * rows_array, int * cols_array, float * vals_array , float * vect, int nnz, int rows) {

	//printf("computation of the COO for the GPU ------------------------------------------- \n");
	// computing the global index of the thread 
	// threadIdx.x : position of the thrad in the block 
	// blockIdx.x : position of the block in the grid 
	// blockDim.x : niumber of threads per block

	int idx = threadIdx.x + blockIdx.x * blockDim.x;

	
	//loop "grid-stride" 
	//each thread takes several elements spaced by : 
	//blockDim.x * gridDim.x = total number of thread started 
	// it allows us covering the array enven if n is huge

	for (int number = idx; number < nnz; number += blockDim.x * gridDim.x) {
			// coalescent memory access : each neighbour thread access to neighours index
			// we use atomicAdd for better computationa because
			// we have some race condition (multiple threads can access to the same line
			// and multiple thread writie in the same cell. 
    		atomicAdd(&res[rows_array[number]], vals_array[number] * vect[cols_array[number]]);
}
    }   

__global__
void computeSpmvCSRGPU(float * res, int * rows_array, int * cols_array, float * vals_array , float * vect, int nnz, int rows, int * row_ptr) {
	/*
	problems without the shared memory + warp 
	if a line has a lot of nnz -> thread slow 
	other threads are waiting 
	divergence important 
	*/
		
	int idx = threadIdx.x + blockIdx.x * blockDim.x;

    //csr to matrix + result for ones
    for (int i = idx;i<rows;i+=blockDim.x*gridDim.x) {
        for (int j = row_ptr[i]; j < row_ptr[i+1];j++) {
                    //printf("element (%d %d), value : %.2f\n",i,cols_array[j],vals_array[j]); 
            res[i] += vals_array[j]*vect[cols_array[j]];
        }
    }

    }

__global__
void computeSpmvCSRWarpGPU(float * res, int * rows_array, int * cols_array, float * vals_array , float * vect, int nnz, int rows, int * row_ptr) {
	/*
	before, 1 thread -> the whole line 
	now, 32 threads -> the same line, perfect if the line is long 
	no need for the shared memory here, warp faster than shared memory 
	for recent GPU 
	here we do not re use the vect so it is useless. 

	*/


	// global thread id
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // warp id global
    int warp_id = idx / warpSize;

    // lane id inside warp
    int lane = threadIdx.x % warpSize;
	
	int row = 0; 
	float sum = 0.0;
	int j,offset;

    // one warp = one row
    if (warp_id < rows) {

        row = warp_id;

        sum = 0.0;

        // each lane processes part of row
        for (j = row_ptr[row] + lane;
             j < row_ptr[row + 1];
             j += warpSize)
        {
            sum += vals_array[j] * vect[cols_array[j]];
        }

        // warp reduction
        for (offset = warpSize / 2;
             offset > 0;
             offset /= 2)
        {
            sum += __shfl_down_sync(0xffffffff, sum, offset);
        }

        // lane 0 writes result
        if (lane == 0) {
            res[row] = sum;
        }
    }

}



__global__ 
void computeSpmvSellGPU(int sliceSize,int rows, int * slice_offsets,int * column_indices,float* values_array,float * ones, float * res_array) {
	
		
    int row  = threadIdx.x + blockIdx.x * blockDim.x;
    // grid-stride on the lines 
	int idx = 0, col = 0; 
	if (row >= rows) return; //in the case there is an overlapp 
	

    int slice  = row/sliceSize;
	int local_row = row % sliceSize; 
	int start = slice_offsets[slice]; 
	int end = slice_offsets[slice+1]; 
	int max_nnz = (end - start) / sliceSize; 
	float sum = 0.0;
	int k; 
	for (k = 0;k<max_nnz;k++) {
		    idx = start + k * sliceSize + local_row;

    		col = column_indices[idx];
    		if (col != -1)
        		sum += values_array[idx] * ones[col];
	}
	res_array[row] = sum; 
}









int main(int argc, char * argv[]) {
    printf("======================================= Device properties ========================================\n");
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("CUDA Device Query...\n");
    printf("There are %d CUDA devices.\n", devCount);

    cudaDeviceProp devProp;
    for (int i = 0; i < devCount; ++i)
    {
        // Get device properties
        printf("\nCUDA Device #%d\n", i);
        cudaGetDeviceProperties(&devProp, i);
        printDevProp(devProp);
    }
	//only in case the user wants to read a file .mtx
	int * row_ptr_mtx; 
	int * cols_array_mtx; 
	float * vals_array_mtx;
	char filename[256];
	int user_matrix = 0;
	MM_typecode matcode;
	int return_val = 0;
	

	unsigned int seed = 19;
	unsigned int sliceSize = 0;
	int rows = 0; 
	int cols = 0;
	int nnz = 0;
	float min = 0.0, max = 5.0;
	int occurRow = 0, actuRow = 0, numberOccur = 0;
	float * ones;
	char * matrix; 
	COOvalue * cooarray; 

	float * vals_array; 
	int * cols_array;
	int * rows_array;
    float timers[NITER];

	// GPU PART (arrays) 
	int * GPU_cols; 
	int * GPU_rows; 
	float * GPU_vals;
	float * GPU_COOres; 
	float * GPU_CSRres; 
	int * GPU_row_ptr; 
	float * GPU_vect; 
	
	int GPU_len = 0;
	int GPU_resLen = 0;

	int maxThreadsPerBlock = devProp.maxThreadsPerBlock;
	int threadsPerBlock = min(1024,maxThreadsPerBlock); 
	printf("threads per block : %d \n",threadsPerBlock);
	
	int numBlocks = 255; //needs to be changed (train with oth    er values)
	int sharedMemSize = 0;
	
	//TIME MEASUREMENT 
	TIMER_DEF;
	float GPU_COO_time,GPU_CSR_time,GPU_SELL_time, CPU_COO_time;	

	// for writing some results 
	FILE * file;
	file =	fopen("res.txt","w");


	int randomRow = 0, randomCol = 0;
	float randomVal = 0.0;
	// so it can be a filename (.mtx)
	if (argc == 3) {
		sliceSize = atoi(argv[2]);
		sprintf(filename, "%s",argv[1]);
		printf("filename : %s \n",filename);
		user_matrix = 1;
		return_val = mm_read_mtx_crd(filename,&rows,&cols,&nnz,&row_ptr_mtx,&cols_array_mtx,&vals_array_mtx,&matcode);
		printf("does the function work ? %d \n",return_val);
		if (return_val != 0) return 1;
		cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 
		vals_array = (float*) calloc(nnz,sizeof(float));
		cols_array =  (int *) calloc(nnz,sizeof(int));
		rows_array = (int *) calloc(nnz,sizeof(int));

		for (int i = 0;i<nnz;i++) {
			//printf("value %d (%d, %d) : %.2f\n",i,row_ptr_mtx[i],cols_array_mtx[i],vals_array_mtx[i]);
			cooarray[i].col = cols_array_mtx[i]-1; 
			cooarray[i].val = vals_array_mtx[i];
			cooarray[i].row = row_ptr_mtx[i]-1;

			cols_array[i] = cols_array_mtx[i]-1; 
			vals_array[i] = vals_array_mtx[i];
			rows_array[i] = row_ptr_mtx[i]-1;

			
		}
		if (!return_val) { 
			if (cols_array_mtx) free(cols_array_mtx);
			if (vals_array_mtx) free(vals_array_mtx);
			if(row_ptr_mtx) free(row_ptr_mtx);
		}
		

	}
	//in case the matrix has to be randomized with sizes written by user.
	else {
		if (argc < 5) {
			printf("Usage: %s <rows> <cols> <nnz> <sliceSize> [SEEED] or %s <path-of-file.mtx> <sliceSize>\n", argv[0],argv[0]);
			return 1;
		}
		rows = atoi(argv[1]);
		cols = atoi(argv[2]);
		nnz = atoi(argv[3]);
		sliceSize = atoi(argv[4]);
		//the case the user added a fifth parameter 
		// and it is a number (seed)
		if (argc == 6) {
			seed = atoi(argv[5]);
			srand(seed);
		}
		else srand(time(NULL));
//		printf("Value of nnz : %d \n",nnz);
		
		matrix = (char*) calloc(rows*cols,sizeof(char)); 
		cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 
		vals_array = (float*) calloc(nnz,sizeof(float));
		cols_array =  (int *) calloc(nnz,sizeof(int));
		rows_array = (int *) calloc(nnz,sizeof(int));


		for (int i = 0;i<nnz;i++) {
			randomVal = random_float(min,max); 
//			printf("Value  %d : %.2f\n", i,randomVal); 
			do {
				randomRow = rand() % (rows); 
				randomCol = rand() % (cols); 
			} while (matrix[randomRow*cols+randomCol]);
			matrix[randomRow*cols+randomCol] = 1; 
			//printf("Row position  : %d\n", randomRow); 
			//printf("Col position  : %d\n", randomCol); 
			cooarray[i].row = randomRow; 
			cooarray[i].col = randomCol;		
			cooarray[i].val = randomVal;

			cols_array[i] = randomCol; 
			vals_array[i] = randomVal;
			rows_array[i] = randomRow;


		}
		free(matrix);
	}

	

	printf("Silce size : %d\n",sliceSize);
	//vector we use for the multiplication
	ones = (float *) malloc(cols*sizeof(float));
	for (int i=0;i<cols;i++) {ones[i] = 1.0;}

	//c) SpMV COO 
	float * cooRes = (float*) calloc(rows,sizeof(float));
	TIMER_START;
	computeSpmvCOO(cooRes,rows_array,cols_array,vals_array,ones,nnz,rows);
	TIMER_STOP; 
	CPU_COO_time = TIMER_ELAPSED; 
	printf("Time of COO (CPU) : %3.5f ms \n",CPU_COO_time*1000);

	// GPU part ==========================
	
	GPU_len = (nnz); GPU_resLen = (rows);
	cudaMalloc(&GPU_rows,GPU_len*sizeof(int));
	cudaMalloc(&GPU_cols,GPU_len*sizeof(int));
	cudaMalloc(&GPU_vals,GPU_len*sizeof(float));
	cudaMalloc(&GPU_vect, (cols)*sizeof(float));
	cudaMalloc(&GPU_COOres,GPU_resLen*sizeof(float)); 

	cudaMemcpy(GPU_rows,rows_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice);
 
	cudaMemcpy(GPU_cols,cols_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice); 
	
	cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(float),cudaMemcpyHostToDevice); 

	cudaMemcpy(GPU_vect,ones,cols*sizeof(float),cudaMemcpyHostToDevice); 

	cudaMemset(GPU_COOres,0,GPU_resLen*sizeof(float));


	//time measurement with cuda event
	numBlocks = (nnz + threadsPerBlock - 1) / threadsPerBlock; //we are working on nnz so the grid size has to be adapted to this unit
	TIMER_START; 
	computeSpmvCOOGPU<<<numBlocks,threadsPerBlock>>>(GPU_COOres, GPU_rows, GPU_cols, GPU_vals , GPU_vect, GPU_len, GPU_resLen);
	cudaDeviceSynchronize();
	TIMER_STOP;
	GPU_COO_time = TIMER_ELAPSED; 
	printf("Time of GPU in COO (host) : %3.5f ms\n",GPU_COO_time*1000);
	
  
	cudaDeviceSynchronize(); //waits for the end of GPU

	float * cooResGPU_Copy = (float *) calloc(rows,sizeof(float));  
  cudaMemcpy(cooResGPU_Copy,GPU_COOres , GPU_resLen*sizeof(float),cudaMemcpyDeviceToHost);	

printf("COO res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
    //   fprintf(file,"%f\n", cooResGPU_Copy[i]);
        }   	
	

	




	//time measurement with cuda event
	cudaMemset(GPU_COOres,0,GPU_resLen*sizeof(float));

	cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 
	cudaEvent_t start_coo, stop_coo;
	cudaEventCreate(&start_coo);	
	cudaEventCreate(&stop_coo);
	cudaEventRecord(start_coo);
	computeSpmvCOOGPU<<<numBlocks,threadsPerBlock>>>(GPU_COOres, GPU_rows, GPU_cols, GPU_vals , GPU_vect, GPU_len, GPU_resLen);
	cudaEventRecord(stop_coo);
	cudaEventSynchronize(stop_coo);
	cudaEventElapsedTime(&GPU_COO_time,start_coo,stop_coo);
	cudaEventDestroy(start_coo);
	cudaEventDestroy(stop_coo);
	printf("Time of GPU in COO (cuda event) : %3.5f ms\n",GPU_COO_time);
	






	// spMV CSR  
	
	int * row_ptr_array = (int *) calloc(rows+1,sizeof(int));;
	qsort(cooarray,nnz,sizeof(COOvalue),compare);
	
	// re-ordering arrays after sorting
	for (int i = 0;i<nnz;i++) {
		vals_array[i] = cooarray[i].val;
		cols_array[i] = cooarray[i].col;
		rows_array[i] = cooarray[i].row;
	}

	free(cooarray);
	
	//copies the re-orgered CPU arrays in the GPU ones
	cudaMemcpy(GPU_rows,rows_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice);
 
	cudaMemcpy(GPU_cols,cols_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice); 
	
	cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(float),cudaMemcpyHostToDevice); 




	float * csrRes = (float*) calloc(rows,sizeof(float));
	
	cudaMalloc(&GPU_CSRres,rows*sizeof(float));

	cudaMemset(GPU_CSRres,0,rows*sizeof(float));
	
	cudaMalloc(&GPU_row_ptr,(rows+1)*sizeof(int));

	
	computeSpmvCSR(csrRes,rows_array,cols_array,vals_array,ones,nnz,rows,row_ptr_array);

	

	cudaMemcpy(GPU_row_ptr,row_ptr_array,(rows+1)*sizeof(int),cudaMemcpyHostToDevice);

	
	numBlocks = (rows + threadsPerBlock - 1) / threadsPerBlock; // we are working on lines, so the grid size has to be adapted
//	numBlocks = ((rows*32) + threadsPerBlock - 1)/threadsPerBlock; // here 1 wrap = 32 threads = 1 row
	
	TIMER_START; 	
	computeSpmvCSRGPU<<<numBlocks,threadsPerBlock>>>(GPU_CSRres,GPU_rows,GPU_cols,GPU_vals,GPU_vect,nnz,rows,GPU_row_ptr);
	cudaDeviceSynchronize();
	TIMER_STOP; 
	GPU_CSR_time = TIMER_ELAPSED;
	printf("TIME OF csr gpu (host) : %3.5f ms \n",GPU_CSR_time*1000);

	//time measurement in cuda event of CSR with cuda event for GPU 
	cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 
	cudaMemset(GPU_CSRres,0,rows*sizeof(float));
	cudaEvent_t start_csr, stop_csr;
	cudaEventCreate(&start_csr);	
	cudaEventCreate(&stop_csr);
	cudaEventRecord(start_csr);

	computeSpmvCSRGPU<<<numBlocks,threadsPerBlock>>>(GPU_CSRres,GPU_rows,GPU_cols,GPU_vals,GPU_vect,nnz,rows,GPU_row_ptr);
	cudaEventRecord(stop_csr);
	cudaEventSynchronize(stop_csr);
	cudaEventElapsedTime(&GPU_CSR_time,start_csr,stop_csr);
	cudaEventDestroy(start_csr);
	cudaEventDestroy(stop_csr);
	printf("Time of GPU in CSR (cuda event) : %3.5f ms\n",GPU_CSR_time);
	


	


	cudaDeviceSynchronize();

	float * csrResGPU_Copy = (float *) calloc(rows,sizeof(float));  
  	cudaMemcpy(csrResGPU_Copy,GPU_CSRres , rows*sizeof(float),cudaMemcpyDeviceToHost);	

	printf("CSR res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
       fprintf(file,"%f\n", csrResGPU_Copy[i]);
        }   	



	
	//spMV SELL 
	
	  float * sellRes = (float*) calloc(rows,sizeof(float));

    //computeSpmvSELL(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes);

    int sizeSellVect = 0;
    int sizeSliceOffset = 0;
    int * slice_offsetsSell;
    int * column_indicesSell;
    float * values_arraySell;
    computeSpmvSELLv2(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes,&column_indicesSell, &values_arraySell,&slice_offsetsSell,&sizeSellVect,&sizeSliceOffset);



    printf("sizeSellVect : %d \n",sizeSellVect);
    printf("sizeSliceOffset : %d \n",sizeSliceOffset);

	// sell GPU 
	int * GPU_sliceOffset; 
	cudaMalloc(&GPU_sliceOffset,sizeSliceOffset*sizeof(int)); 

	int * GPU_column_indicesSell; 
	cudaMalloc(&GPU_column_indicesSell,sizeSellVect*sizeof(int));

	float * GPU_values_arraySell; 
	cudaMalloc(&GPU_values_arraySell,sizeSellVect*sizeof(float));

	float * GPU_SELLres;
	cudaMalloc(&GPU_SELLres,rows*sizeof(float));

	cudaMemset(GPU_SELLres,0,rows*sizeof(float));

	cudaMemcpy(GPU_sliceOffset,slice_offsetsSell,sizeSliceOffset*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(GPU_column_indicesSell,column_indicesSell,sizeSellVect*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(GPU_values_arraySell,values_arraySell,sizeSellVect*sizeof(float),cudaMemcpyHostToDevice);
	
	numBlocks = (rows + threadsPerBlock - 1) / threadsPerBlock; //based on rows, so we change the numBlocks

	TIMER_START;

	cudaError_t err; 
	computeSpmvSellGPU<<<numBlocks,threadsPerBlock>>>(sliceSize,rows, GPU_sliceOffset,GPU_column_indicesSell,GPU_values_arraySell,GPU_vect, GPU_SELLres);
	err = cudaGetLestError(); 
	if (err != cudaSuccess) {
   	 	fprintf(stderr, "Kernel launch error: %s\n", cudaGetErrorString(err));
   		 exit(EXIT_FAILURE);
	}
	err = cudaDeviceSynchronize();
	if (err != cudaSuccess) {
		fprintf(stderr, "Kernel execution error: %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}


	cudaDeviceSynchronize();
	TIMER_STOP;
	GPU_SELL_time = TIMER_ELAPSED;
	printf("Time of GPU on SELL (Host) : %3.5f ms \n",GPU_SELL_time*1000);

	cudaDeviceSynchronize();

	float * sellResGPU_Copy = (float *) calloc(rows,sizeof(float));  
  	cudaMemcpy(sellResGPU_Copy, GPU_SELLres, rows*sizeof(float),cudaMemcpyDeviceToHost);	

	printf("SELL res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
     //   fprintf(file,"%f\n", sellResGPU_Copy[i]);
        }   	

	cudaDeviceSynchronize();
	cudaMemset(GPU_SELLres,0,rows*sizeof(float));

	cudaDeviceSynchronize();
	cudaEvent_t start_sell, stop_sell;
	cudaEventCreate(&start_sell);	
	cudaEventCreate(&stop_sell);
	cudaEventRecord(start_sell);
	computeSpmvSellGPU<<<numBlocks,threadsPerBlock>>>(sliceSize,rows, GPU_sliceOffset,GPU_column_indicesSell,GPU_values_arraySell,GPU_vect, GPU_SELLres);
	cudaEventRecord(stop_sell);
	cudaEventSynchronize(stop_sell);
	cudaEventElapsedTime(&GPU_SELL_time,start_sell,stop_sell);
	cudaEventDestroy(start_sell);
	cudaEventDestroy(stop_sell);
	printf("Time of GPU in SELL (cuda event) : %f ms\n",GPU_SELL_time);
	
	

	//CPU final comparisons =========================
	computeDiffArrays(cooRes , csrRes, rows, "COO CPU","CSR CPU") ;
	computeDiffArrays(cooRes , sellRes, rows, "COO CPU","SELL CPU") ;

	//GPU final comparisons =========================
	computeDiffArrays(cooRes , cooResGPU_Copy, rows, "COO CPU","COO GPU") ;
	computeDiffArrays(csrRes , csrResGPU_Copy, rows, "CSR CPU","CSR GPU") ;
	computeDiffArrays(sellRes , sellResGPU_Copy, rows, "SELL CPU","SELL GPU") ;

	//======================== GPU COO cusparse and GPU COO personnal result comparison	
	cuSparseCOOComparison(cooRes, GPU_COOres,GPU_rows,GPU_cols,GPU_vals, GPU_vect,nnz,rows,cols);

	cuSparseCSRComparison(csrRes,GPU_CSRres,GPU_row_ptr,GPU_cols,GPU_vals,GPU_vect,nnz,rows,cols);


    /* Note, free the memory */
//	free(row_ptr_array);
	free(rows_array);
	free(cols_array);
	free(vals_array);
	free(sellRes);
	free(cooRes); 
	free(csrRes);
	free(ones);
 	free(slice_offsetsSell);
    free(column_indicesSell);
    free(values_arraySell);

	// GPU FREE 
	

	cudaFree(GPU_rows); 
	cudaFree(GPU_vals);
	cudaFree(GPU_cols);
	cudaFree(GPU_COOres);
	cudaFree(GPU_row_ptr);
	cudaFree(GPU_CSRres); 
	cudaFree(GPU_SELLres);
	cudaFree(GPU_vect);

	fclose(file);




	return 0;

}

