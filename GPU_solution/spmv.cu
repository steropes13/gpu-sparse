#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>

#include "include/mmio.h"
#include "include/my_time_lib.h"
#include "include/spmvCPU.h"



#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float



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
void computeSpmvCOOGPU(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows) {

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
void computeSpmvCSRGPU(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows, int * row_ptr) {

		
	int idx = threadIdx.x + blockIdx.x * blockDim.x;

    //csr to matrix + result for ones
    for (int i = idx;i<rows;i+=blockDim.x*gridDim.x) {
        for (int j = row_ptr[i]; j < row_ptr[i+1];j++) {
                    //printf("element (%d %d), value : %.2f\n",i,cols_array[j],vals_array[j]); 
            res[i] += vals_array[j]*vect[cols_array[j]];
        }
    }

    }














int main(int argc, char * argv[]) {
    printf("======================================= Device properties ========================================\n");
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("CUDA Device Query...\n");
    printf("There are %d CUDA devices.\n", devCount);

    for (int i = 0; i < devCount; ++i)
    {
        // Get device properties
        printf("\nCUDA Device #%d\n", i);
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, i);
        printDevProp(devProp);
    }
	//only in case the user wants to read a file .mtx
	int * row_ptr_mtx; 
	int * cols_array_mtx; 
	double * vals_array_mtx;
	char filename[256];
	int user_matrix = 0;
	MM_typecode matcode;
	int return_val = 0;
	

	unsigned int seed = 19;
	unsigned int sliceSize = 0;
	int rows = 0; 
	int cols = 0;
	int nnz = 0;
	double min = 0.0, max = 5.0;
	int occurRow = 0, actuRow = 0, numberOccur = 0;
	double * ones;
	char * matrix; 
	COOvalue * cooarray; 

	double * vals_array; 
	int * cols_array;
	int * rows_array;
    double timers[NITER];

	// GPU PART (arrays) 
	int * GPU_cols; 
	int * GPU_rows; 
	double * GPU_vals;
	double * GPU_COOres; 
	double * GPU_CSRres; 
	int * GPU_row_ptr; 
	double * GPU_vect; 
	
	int GPU_len = 0;
	int GPU_resLen = 0;
	int block_size = 256; 
	int grid_size = 1; //needs to be changed (train with oth    er values)
	
	


	int randomRow = 0, randomCol = 0;
	double randomVal = 0.0;
	// so it can be a filename (.mtx)
	if (argc == 3) {
		sliceSize = atoi(argv[2]);
		sprintf(filename, "%s",argv[1]);
		printf("filename : %s \n",filename);
		user_matrix = 1;
		return_val = mm_read_mtx_crd(filename,&rows,&cols,&nnz,&row_ptr_mtx,&cols_array_mtx,&vals_array_mtx,&matcode);
		printf("does the function work ? %d \n",return_val);
		cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 
		vals_array = (double*) calloc(nnz,sizeof(double));
		cols_array =  (int *) calloc(nnz,sizeof(int));
		rows_array = (int *) calloc(nnz,sizeof(int));

		for (int i = 0;i<nnz;i++) {
			printf("value %d (%d, %d) : %.2f\n",i,row_ptr_mtx[i],cols_array_mtx[i],vals_array_mtx[i]);
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
		vals_array = (double*) calloc(nnz,sizeof(double));
		cols_array =  (int *) calloc(nnz,sizeof(int));
		rows_array = (int *) calloc(nnz,sizeof(int));


		for (int i = 0;i<nnz;i++) {
			randomVal = random_double(min,max); 
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
	ones = (double *) malloc(cols*sizeof(double));
	for (int i=0;i<cols;i++) {ones[i] = 1.0;}

	//c) SpMV COO 
	double * cooRes = (double*) calloc(rows,sizeof(double));
		computeSpmvCOO(cooRes,rows_array,cols_array,vals_array,ones,nnz,rows);


	// GPU part ==========================
	
	GPU_len = (nnz); 
	GPU_resLen = (rows);
	cudaMalloc(&GPU_rows,GPU_len*sizeof(int));
	cudaMalloc(&GPU_cols,GPU_len*sizeof(int));
	cudaMalloc(&GPU_vals,GPU_len*sizeof(double));
	cudaMalloc(&GPU_vect, (cols<<1)*sizeof(double));
	cudaMalloc(&GPU_COOres,GPU_resLen*sizeof(double)); 

	cudaMemcpy(GPU_rows,rows_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice);
 
	cudaMemcpy(GPU_cols,cols_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice); 
	
	cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(double),cudaMemcpyHostToDevice); 


	cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(double),cudaMemcpyHostToDevice); 


	cudaMemcpy(GPU_vect,ones,cols*sizeof(double),cudaMemcpyHostToDevice); 

	cudaMemset(GPU_COOres,0,GPU_resLen*sizeof(double));




	computeSpmvCOOGPU<<<grid_size,block_size>>>(GPU_COOres, GPU_rows, GPU_cols, GPU_vals , GPU_vect, GPU_len, GPU_resLen);

  cudaMemcpy(cooRes,GPU_COOres , GPU_resLen*sizeof(double),cudaMemcpyDeviceToHost);	

printf("COO res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
        printf("y[%d] = %f\n", i, cooRes[i]);
        }   	


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
	
	cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(double),cudaMemcpyHostToDevice); 




	double * csrRes = (double*) calloc(rows,sizeof(double));
	
	cudaMalloc(&GPU_CSRres,rows*sizeof(double));
	
	cudaMalloc(&GPU_row_ptr,(rows+1)*sizeof(int));
	computeSpmvCSR(csrRes,rows_array,cols_array,vals_array,ones,nnz,rows,row_ptr_array);

	

	cudaMemcpy(GPU_row_ptr,row_ptr_array,(rows+1)*sizeof(int),cudaMemcpyHostToDevice);


	computeSpmvCSRGPU<<<grid_size,block_size>>>(GPU_CSRres,GPU_rows,GPU_cols,GPU_vals,GPU_vect,nnz,rows,GPU_row_ptr);


  	cudaMemcpy(csrRes,GPU_CSRres , rows*sizeof(double),cudaMemcpyDeviceToHost);	

	printf("CSR res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
        printf("y[%d] = %f\n", i, csrRes[i]);
        }   	



	
	//spMV SELL 
	
	///void computeSpmvSELL(int nbSlices, COOvalue * cooArray, int rows, int cols) {
//	double * sellRes = (double*) calloc(rows,sizeof(double));

//	computeSpmvSELL(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes);

	





  

    /* Note, free the memory */
//	free(row_ptr_array);
	free(rows_array);
	free(cols_array);
	free(vals_array);
	//free(sellRes);
	free(cooRes); 
//	free(csrRes);
	free(ones);

	// GPU FREE 
	

	cudaFree(GPU_rows); 
	cudaFree(GPU_vals);
	cudaFree(GPU_cols);
	cudaFree(GPU_COOres);
	cudaFree(GPU_vect);




	return 0;

}

