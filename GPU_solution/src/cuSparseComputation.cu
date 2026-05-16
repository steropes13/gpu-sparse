#include <stdio.h> 
#include <stdlib.h>
#include <cusparse.h>
#include "../include/cuSparseComputation.cuh"


void computeDiffArrays(float * cpuVect , float * gpuVect, int n, char * s1, char * s2) {
    float relativeError;
    int isDiff = 0;
    const float epsilon = 1e-4f;
	int counter = 0;

    for (int i = 0; i < n; i++) {
        float diff = fabs(gpuVect[i] - cpuVect[i]);
        float cpuValueAbs = fabs(cpuVect[i]);

        if (cpuValueAbs < 1e-12) {
            relativeError = diff;  // in case the cpuValue is really close to 0  
								   // we compute the absolute error
        } else {
            relativeError = diff / cpuValueAbs; // if the cpuValue is not close to 0 we can compute 
										  // the relative error 
        }

        if (relativeError > epsilon) {
            fprintf(stderr,
                "ERROR for %s VS %s  at i=%d : CPU=%f GPU=%f relErr=%e\n",
                s1,s2,i, cpuVect[i], gpuVect[i], relativeError);
            isDiff = 1;
            break; // stops the loop after the first error encountered
        }
    }

    if (!isDiff)
        fprintf(stderr, "OK: %s vs %s within tolerance %.1e\n", s1, s2, epsilon);
}

int cuSparseCOOComparison(float * cooRes, float * GPU_COOres, int * GPU_rows,int * GPU_cols,float * GPU_vals,float * GPU_vect, int nnz,int rows,int cols) {
		FILE * file; 
		file = fopen("res_COO_cusparse.txt","w");
        cusparseHandle_t     handle = NULL;
        cusparseSpMatDescr_t matA;
        cusparseDnVecDescr_t vecX, vecY;
        void*                dBuffer    = nullptr;
        size_t               bufferSize = 0;
        float alpha = 1.0;
        float beta = 0.0;
        float GPU_CUSPARSE_time = 0;

		cudaEvent_t start, stop;
    	cudaEventCreate(&start);
    	cudaEventCreate(&stop);


		cudaMemset(GPU_COOres,0,rows*sizeof(float));
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

		//printf("size of buffer : %d \n",bufferSize);
		//printf("nnz = %d, rows = %d, cols = %d, alpha = %f, beta = %f \n", nnz, rows, cols,alpha,beta);
		if (bufferSize >0) CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) ) 

		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 

    	// starting of the time measurement 
    	cudaEventRecord(start);


		// Allocate the buffer and execute the spMV
		CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) );

		// the end of time measuring
		cudaEventRecord(stop);
		cudaEventSynchronize(stop); // waits the end of time register

		// time computing in miliseconds
		cudaEventElapsedTime(&GPU_CUSPARSE_time, start, stop);
		printf("Time execution of cusparseSpMV COO : %f ms\n", GPU_CUSPARSE_time);

		// cleaning the cuda events
		cudaEventDestroy(start);
		cudaEventDestroy(stop);



		// cleaning 
		if (dBuffer != nullptr) CHECK_CUDA( cudaFree(dBuffer) );

			// destroy matrix/vector descriptors
		CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )
		CHECK_CUSPARSE( cusparseDestroy(handle) )

		cudaDeviceSynchronize(); //waits for the end of GPU
		float * cusparseCopy = (float *)  malloc(rows*sizeof(float));
	  cudaMemcpy(cusparseCopy,GPU_COOres , rows*sizeof(float),cudaMemcpyDeviceToHost);	
		computeDiffArrays(cusparseCopy, cooRes, rows, "GPU cuspare COO", "CPU COO"); 
		

 printf("COO res (GPU cusparse) =========== :\n");
    for (int i = 0; i < rows; i++) {
        fprintf(file,"%f\n", cusparseCopy[i]);
        }   	
	
//====================================================================
	
		fclose(file);

		free(cusparseCopy);

		return 0; //by default we return a value because each macro of cusparse return a value 
		// (CHECK_USPARSE and CHECK_CUDA) 

}


int cuSparseCSRComparison(float * csrRes, float * GPU_CSRres, int * GPU_row_ptr,int * GPU_cols,float * GPU_vals,float * GPU_vect, int nnz,int rows,int cols) {
		FILE * file; 
		file = fopen("res_CSR_cusparse.txt","w");

        cusparseHandle_t     handle = NULL;
        cusparseSpMatDescr_t matA;
        cusparseDnVecDescr_t vecX, vecY;
        void*                dBuffer    = nullptr;
        size_t               bufferSize = 0;
        float alpha = 1.0;
        float beta = 0.0;
        float GPU_CUSPARSE_time = 0;

    	// Initialisation des événements CUDA pour mesurer le temps
    	cudaEvent_t start, stop;
    	cudaEventCreate(&start);
    	cudaEventCreate(&stop);

		cudaMemset(GPU_CSRres,0,rows*sizeof(float));
		CHECK_CUSPARSE( cusparseCreate(&handle) )

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

		//printf("nnz = %d, rows = %d, cols = %d, alpha = %f, beta = %f \n", nnz, rows, cols,alpha,beta);
		if (bufferSize >0) CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) ) 

		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 

    	// starting of the time measurement 
    	cudaEventRecord(start);


		// Allocate the buffer and execute the spMV
		CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) );


		    // end of time measurement 
		cudaEventRecord(stop);
		cudaEventSynchronize(stop); // waits for the end of time register

		// time computation in miliseconds
		cudaEventElapsedTime(&GPU_CUSPARSE_time, start, stop);
		printf("execution time of cusparseSpMV CSR : %f ms\n", GPU_CUSPARSE_time);

		// cleaning cuda events
		cudaEventDestroy(start);
		cudaEventDestroy(stop);


		// cleaning 
		if (dBuffer != nullptr) CHECK_CUDA( cudaFree(dBuffer) );

			// destroy matrix/vector descriptors
		CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
		CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )
		CHECK_CUSPARSE( cusparseDestroy(handle) )

		cudaDeviceSynchronize(); //waits for the end of GPU and avoid last kernel "pollution" 


		float * cusparseCopy = (float *) malloc(rows*sizeof(float));
	  cudaMemcpy(cusparseCopy,GPU_CSRres , rows*sizeof(float),cudaMemcpyDeviceToHost);	
		computeDiffArrays(cusparseCopy, csrRes, rows, "GPU cuspare CSR", "CPU CSR"); 
		computeDiffArrays(cusparseCopy, csrRes, rows, "GPU cuspare CSR", "CPU CSR"); 
	    for (int i = 0; i < rows; i++) {
        	fprintf(file,"%f\n", cusparseCopy[i]);
        }   	

		
		fclose(file);	

		free(cusparseCopy);
		return 0; //by default we return a value because each macro of cusparse return a value 
		// (CHECK_USPARSE and CHECK_CUDA) 

}
