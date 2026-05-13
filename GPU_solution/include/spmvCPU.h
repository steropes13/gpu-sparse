#ifndef SPMV_H
#define SPMV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


float random_float(float min, float max);
int random_int(int min, int max);
int compare(const void * a, const void * b);

typedef struct {
    int col;
    int row;
    float val;
} COOvalue;

void computeSpmvCOO(float * res, int * rows_array, int * cols_array, float * vals_array , float * vect, int nnz, int rows);

void computeSpmvCSR(float * res, int * rows_array, int * cols_array, float * vals_array , float * vect, int nnz, int rows, int * row_ptr);

void computeSpmvSELL(int sliceSize,int nnz, int * rows_array,int * cols_array, float * vals_array, int rows, int cols, int * row_ptr,float * ones,float * res_array);

void computeSpmvSELLv2(int sliceSize,int nnz, int * rows_array,int * cols_array, float * vals_array, int rows, int cols, int * row_ptr,float * ones,float * res_array,int ** column_indices, float ** values_array, int ** slice_offsets, int * sizeVect, int * sizeOffset);


#ifdef __cplusplus
}
#endif

#endif
