#ifndef SPMV_H
#define SPMV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WARMUP 2
#define NITER 10

double random_double(double min, double max);
int random_int(int min, int max);
int compare(const void * a, const void * b);

typedef struct {
    int col;
    int row;
    double val;
} COOvalue;

void computeSpmvCOO(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows);

void computeSpmvCSR(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows, int * row_ptr);

void computeSpmvSELL(int sliceSize,int nnz, int * rows_array,int * cols_array, double * vals_array, int rows, int cols, int * row_ptr,double * ones,double * res_array);

void computeSpmvSELLv2(int sliceSize,int nnz, int * rows_array,int * cols_array, double * vals_array, int rows, int cols, int * row_ptr,double * ones,double * res_array,int ** column_indices, double ** values_array, int ** slice_offsets, int * sizeVect, int * sizeOffset);


#ifdef __cplusplus
}
#endif

#endif
