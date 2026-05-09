#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include "include/mmio.h"
#include "include/my_time_lib.h"
#include "include/spmvCPU.h"


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

enum spmvType {
	COO,
	CSR,
	ELL,
	SELL
};



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

	//TIME MEASUREMENT 
	TIMER_DEF;
	float GPU_COO_time,GPU_CSR_time,GPU_SELL_time, CPU_COO_time;	

	
	int randomRow = 0, randomCol = 0;
	float randomVal = 0.0;
	// so it can be a filename (.mtx)

	//writing in a file 
	FILE * file; 
	file = fopen("res.txt", "w");
	enum spmvType computationType = COO;



	if (argc == 4) {
		sliceSize = atoi(argv[2]);
		computationType = (spmvType) atoi(argv[3]);
		printf("spmV type : %d\n",computationType);
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
		if (argc < 6) {
			printf("Usage: %s <rows> <cols> <nnz> <sliceSize> <spmvType> [SEEED] or %s <path-of-file.mtx> <sliceSize> <spmvType>\n", argv[0],argv[0]);
			return 1;
		}
		rows = atoi(argv[1]);
		cols = atoi(argv[2]);
		nnz = atoi(argv[3]);
		sliceSize = atoi(argv[4]);
		computationType = (spmvType) atoi(argv[5]);	
		printf("spmV type : %d\n",computationType);
		//the case the user added a fifth parameter 
		// and it is a number (seed)
		if (argc == 7) {
			seed = atoi(argv[6]);
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

	// SpMV COO 
	float * cooRes = (float*) calloc(rows,sizeof(float));



// CUSPARSE APIs ====================================================================
		cusparseHandle_t     handle = NULL;
		cusparseSpMatDescr_t matA;
		cusparseDnVecDescr_t vecX, vecY;
		void*                dBuffer    = nullptr;
		size_t               bufferSize = 0;
		float alpha = 1.0f;
		float beta = 0.0f;
		float GPU_CUSPARSE_time = 0;
	

		//arrays allocation 
		GPU_len = (nnz); GPU_resLen = (rows);
		cudaMalloc(&GPU_rows,GPU_len*sizeof(int));
		cudaMalloc(&GPU_cols,GPU_len*sizeof(int));
		cudaMalloc(&GPU_vals,GPU_len*sizeof(float));
		cudaMalloc(&GPU_vect, (cols)*sizeof(float));
		cudaMalloc(&GPU_COOres,GPU_resLen*sizeof(float)); 



	if (computationType == COO) {
		computeSpmvCOO(cooRes,rows_array,cols_array,vals_array,ones,nnz,rows);

		cudaMemcpy(GPU_rows,rows_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice);
	 
		cudaMemcpy(GPU_cols,cols_array,GPU_len*sizeof(int),cudaMemcpyHostToDevice); 
		
		cudaMemcpy(GPU_vals,vals_array,GPU_len*sizeof(float),cudaMemcpyHostToDevice); 

		cudaMemcpy(GPU_vect,ones,cols*sizeof(float),cudaMemcpyHostToDevice); 
	
		
		cudaMemset(GPU_COOres,0,GPU_resLen*sizeof(float));


		CHECK_CUSPARSE( cusparseCreate(&handle) )
				// Create sparse matrix A in COO format
		CHECK_CUSPARSE( cusparseCreateCoo(&matA, rows, cols, nnz,
										  GPU_rows, GPU_cols, GPU_vals,
										  CUSPARSE_INDEX_32I,
										  CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F) )      


		// Create dense vector X
		CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, cols, GPU_vect, CUDA_R_32F) )
		// Create dense vector y
		CHECK_CUSPARSE( cusparseCreateDnVec(&vecY,rows, GPU_COOres, CUDA_R_32F) )

		// allocate an external buffer if needed
		CHECK_CUSPARSE( cusparseSpMV_bufferSize(
									 handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
									 &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
									 CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )

		printf("size of buffer : %d \n",bufferSize);
		printf("nnz = %d, rows = %d, cols = %d, alpha = %f, beta = %f \n", nnz, rows, cols,alpha,beta);
		if (bufferSize >0) CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) ) 


		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 

		// Allocate the buffer and execute the spMV
		CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) );




		// cleaning 
		if (dBuffer != nullptr) CHECK_CUDA( cudaFree(dBuffer) );

			// destroy matrix/vector descriptors
		CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )
		CHECK_CUSPARSE( cusparseDestroy(handle) )

		cudaDeviceSynchronize(); //waits for the end of GPU
	  cudaMemcpy(cooRes,GPU_COOres , GPU_resLen*sizeof(float),cudaMemcpyDeviceToHost);	
		

 printf("COO res (GPU cusparse) =========== :\n");
    for (int i = rows-10; i < rows; i++) {
        fprintf(file,"y[%d] = %f\n", i, cooRes[i]);
        }   	
//====================================================================
	}
	// spMV CSR  
	
	int * row_ptr_array = (int *) calloc(rows+1,sizeof(int));
	
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

	
	cudaMemcpy(GPU_vect,ones,cols*sizeof(float),cudaMemcpyHostToDevice); 
	if (computationType == CSR) {
		CHECK_CUSPARSE( cusparseCreate(&handle) )
		computeSpmvCSR(csrRes,rows_array,cols_array,vals_array,ones,nnz,rows,row_ptr_array);
		
		for (int i = 0; i < rows+1; i++) {
    		printf("row_ptr[%d] = %d\n", i, row_ptr_array[i]);
}
		cudaMemcpy(GPU_row_ptr,row_ptr_array,(rows+1)*sizeof(int),cudaMemcpyHostToDevice);

printf("row_ptr last = %d (should be nnz=%d)\n", row_ptr_array[rows], nnz);

for (int i=0;i<5;i++) {
    printf("row_ptr[%d]=%d\n", i, row_ptr_array[i]);
    printf("col[%d]=%d val=%f\n", i, cols_array[i], vals_array[i]);
}

    CHECK_CUSPARSE( cusparseCreateCsr(&matA,rows, cols, nnz,
                                      GPU_row_ptr, GPU_cols, GPU_vals,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F) )

		// Create dense vector X
		CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, cols, GPU_vect, CUDA_R_32F) )
		// Create dense vector y
		CHECK_CUSPARSE( cusparseCreateDnVec(&vecY,rows, GPU_CSRres, CUDA_R_32F) )


	// allocate an external buffer if needed
		CHECK_CUSPARSE( cusparseSpMV_bufferSize(
									 handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
									 &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
									 CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )

		printf("size of buffer : %d \n",bufferSize);
		printf("nnz = %d, rows = %d, cols = %d, alpha = %f, beta = %f \n", nnz, rows, cols,alpha,beta);
		if (bufferSize >0) CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) ) 

		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 


		// Allocate the buffer and execute the spMV
		CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) );
		// cleaning 
		if (dBuffer != nullptr) CHECK_CUDA( cudaFree(dBuffer) );

			// destroy matrix/vector descriptors
		CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )
		CHECK_CUSPARSE( cusparseDestroy(handle) )

		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 

	  cudaMemcpy(csrRes,GPU_CSRres , GPU_resLen*sizeof(float),cudaMemcpyDeviceToHost);	
		

 printf("CSR res (GPU cusparse) =========== :\n");
    for (int i = 0; i < rows; i++) {
        fprintf(file,"y[%d] = %f\n", i, csrRes[i]);
        }   	


	}

/*	
	//numBlocks = (rows + blockSize - 1) / blockSize; // we are working on lines, so the grid size has to be adapted
	numBlocks = ((rows*32) + blockSize - 1)/blockSize; // here 1 wrap = 32 threads = 1 row
	
	TIMER_START; 	
	computeSpmvCSRWarpGPU<<<numBlocks,blockSize>>>(GPU_CSRres,GPU_rows,GPU_cols,GPU_vals,GPU_vect,nnz,rows,GPU_row_ptr);
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

	computeSpmvCSRWarpGPU<<<numBlocks,blockSize>>>(GPU_CSRres,GPU_rows,GPU_cols,GPU_vals,GPU_vect,nnz,rows,GPU_row_ptr);
	cudaEventRecord(stop_csr);
	cudaEventSynchronize(stop_csr);
	cudaEventElapsedTime(&GPU_CSR_time,start_csr,stop_csr);
	cudaEventDestroy(start_csr);
	cudaEventDestroy(stop_csr);
	printf("Time of GPU in CSR (cuda event) : %3.5f ms\n",GPU_CSR_time);
	


	


	cudaDeviceSynchronize();
  	cudaMemcpy(csrRes,GPU_CSRres , rows*sizeof(float),cudaMemcpyDeviceToHost);	

	printf("CSR res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
      // printf("y[%d] = %f\n", i, csrRes[i]);
        }   	


*/
	
	//spMV SELL 
	
	///void computeSpmvSELL(int nbSlices, COOvalue * cooArray, int rows, int cols) {
//	float * sellRes = (float*) calloc(rows,sizeof(float));

//	computeSpmvSELL(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes);

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
	
	/*	
	numBlocks = (rows + blockSize - 1) / blockSize; //based on rows, so we change the numBlocks

	TIMER_START;
	computeSpmvSellGPU<<<numBlocks,blockSize>>>(sliceSize,rows, GPU_sliceOffset,GPU_column_indicesSell,GPU_values_arraySell,GPU_vect, GPU_SELLres);
	cudaDeviceSynchronize();
	TIMER_STOP;
	GPU_SELL_time = TIMER_ELAPSED;
	printf("Time of GPU on SELL (Host) : %3.5f ms \n",GPU_SELL_time*1000);

	cudaDeviceSynchronize();
  	cudaMemcpy(sellRes, GPU_SELLres, rows*sizeof(float),cudaMemcpyDeviceToHost);	

	printf("SELL res (GPU) =========== :\n");
    for (int i = 0; i < rows; i++) {
        //printf("y[%d] = %f\n", i, sellRes[i]);
        }   	

	cudaDeviceSynchronize();
	cudaMemset(GPU_SELLres,0,rows*sizeof(float));

	cudaDeviceSynchronize();
	cudaEvent_t start_sell, stop_sell;
	cudaEventCreate(&start_sell);	
	cudaEventCreate(&stop_sell);
	cudaEventRecord(start_sell);
	computeSpmvSellGPU<<<numBlocks,blockSize>>>(sliceSize,rows, GPU_sliceOffset,GPU_column_indicesSell,GPU_values_arraySell,GPU_vect, GPU_SELLres);
	cudaEventRecord(stop_sell);
	cudaEventSynchronize(stop_sell);
	cudaEventElapsedTime(&GPU_SELL_time,start_sell,stop_sell);
	cudaEventDestroy(start_sell);
	cudaEventDestroy(stop_sell);
	printf("Time of GPU in SELL (cuda event) : %f ms\n",GPU_SELL_time);
	
	



*/

  

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

