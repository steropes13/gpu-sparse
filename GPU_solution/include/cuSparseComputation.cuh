#ifndef __CUSPARSE__CUH
#define __CUSPARSE__CUH
#include <cusparse.h> 
#include <stdio.h>
#include <stdlib.h> 

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


void computeDiffArrays(float * a , float * b, int n,  char * s1, char * s2);


int cuSparseCOOComparison(float *cooRes, float *GPU_COOres, int *GPU_rows, int *GPU_cols, float *GPU_vals, float *GPU_vect, int nnz, int rows, int cols);
int cuSparseCSRComparison(float *csrRes, float *GPU_CSRres, int *GPU_row_ptr, int *GPU_cols, float *GPU_vals, float *GPU_vect, int nnz, int rows, int cols);


#endif


