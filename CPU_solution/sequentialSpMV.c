#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/mmio.h" 

#define WARMUP 2
#define NITER 10

#include "include/my_time_lib.h"

double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

int random_int(int min, int max) {
   	return min + (max - min) * rand();
}

typedef struct {
	int col;
	int row; 
	double val;	
} COOvalue;

	

int compare(const void * a, const void * b) {
	COOvalue * ae = (COOvalue *) a; 
	COOvalue * be = (COOvalue *) b; 
	if (ae->row < be->row) return -1; 
	if (ae->row > be->row) return 1;

	if (ae->col < be->col) return -1;
	if (ae->col > be->col) return 1; 
	
	return 0; 
	

}

void computeSpmvCOO(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows) {
	for (int number = 0 ; number < nnz ; number++) {
		res[rows_array[number]] += vals_array[number]*vect[cols_array[number]];
	}

	printf("COO res :\n");
	for (int i = 0; i < rows; i++) {
    	printf("y[%d] = %f\n", i, res[i]);
		}
	}


void computeSpmvCSR(double * res, int * rows_array, int * cols_array, double * vals_array , double * vect, int nnz, int rows, int * row_ptr) {
	int * occurenceArray = calloc(rows, sizeof(int));
//	printf("sorting the arrys : \n");
	for (int i = 0; i < nnz;i++) {
	//	printf("value  : %.2f\n",vals_array[i]); 
	//	printf("Row position  : %d\n", rows_array[i]); 
	//	printf("Col position  : %d\n", cols_array[i]); 
	}
	
	//number of elements for each row 
	for (int j = 0; j <nnz;j++) {
		occurenceArray[rows_array[j]] += 1;
	}

	//printing the occurence array
	for (int j = 0; j <rows;j++) {
	//	printf("row %d : %d \n",j, occurenceArray[j]);
	}

	
	//prefix sum 
	row_ptr[0] = 0; 
	for (int j = 1; j <rows+1;j++) {
		row_ptr[j] = (row_ptr[j-1]) + occurenceArray[j-1];
	//	printf("prefix sum  %d : %d \n",j, (row_ptr[j]));
	}	
	free(occurenceArray);	
	//printf("final value of prefix-sum : %d \n",(row_ptr[rows]));

	//csr to matrix + result for ones
	for (int i = 0;i<rows;i++) {
		for (int j = row_ptr[i]; j < row_ptr[i+1];j++) {
					//printf("element (%d %d), value : %.2f\n",i,cols_array[j],vals_array[j]); 
			res[i] += vals_array[j]*vect[cols_array[j]];
		}
	}

	printf("result csr : \n"); 
	for (int i = 0;i<rows;i++) {
    	printf("y[%d] = %f\n", i, res[i]);
	}

}

void computeSpmvSELLv2(int sliceSize,int nnz, int * rows_array,int * cols_array, double * vals_array, int rows, int cols, int * row_ptr,double * ones,double * res_array,int ** column_indices, double ** values_array, int ** slice_offsets, int * sizeVect, int * sizeOffset ) {
	int nbSlices = 0; 
	if (sliceSize > rows || sliceSize <= 0) sliceSize = rows; //ELLPACK by default
	nbSlices = (rows + sliceSize - 1) / sliceSize;
	printf("nb slices : %d\n",nbSlices);	
    // number of rows per slice
    int rowsPerSlice = sliceSize;

    //Using  CSR for SELL (sorted) 
	// easier to access to the NNZs per row

    int *col_ind = (int *) calloc(nnz, sizeof(int));
    double *val   = (double *) calloc(nnz, sizeof(double));
    //int *offset = (int *) calloc(rows, sizeof(int));
   // for (int i = 0; i < rows; i++) offset[i] = row_ptr[i];

    for (int i = 0; i < nnz; i++) {
        int r = rows_array[i];
    //    int pos = offset[r]++;
        col_ind[i] = cols_array[i];
        val[i]     = vals_array[i];
    }

   // free(offset);

    //compute slice_offsets 

	*slice_offsets = (int *) calloc(nbSlices+1, sizeof(int));
	(*slice_offsets)[0] = 0;
	*sizeOffset = nbSlices+1;

	for (int s = 0; s < nbSlices; s++) {

		int start = s * rowsPerSlice;
		int end   = start + rowsPerSlice;
		if (end > rows) end = rows;

		int max_nnz = 0;

		// compute max nnz in slice (lignes réelles seulement)
		for (int r = start; r < end; r++) {
			int nnz_row = row_ptr[r+1] - row_ptr[r];
			if (nnz_row > max_nnz)
				max_nnz = nnz_row;
		}

		(*slice_offsets)[s+1] = (*slice_offsets)[s] + max_nnz * rowsPerSlice;

		printf("slice %d : max_nnz=%d offset=%d\n",
			   s, max_nnz, (*slice_offsets)[s]);
	}
    //allocate SELL 
	
	int vectorSize = (*slice_offsets)[nbSlices]; 
printf("vector size : %d \n",(*slice_offsets)[nbSlices]);
    *column_indices = (int *) calloc(vectorSize, sizeof(int));
    *values_array = (double *) calloc(vectorSize, sizeof(double));
	*sizeVect = vectorSize;
   

    // fill SELL 

    int pos = 0;

  for (int s = 0; s < nbSlices; s++) {

    int start = s * rowsPerSlice;
    int slice_start = (*slice_offsets)[s];
    int slice_end   = (*slice_offsets)[s+1];

    int pos = slice_start;

    int max_nnz_slice = (slice_end - slice_start) / rowsPerSlice;

    for (int k = 0; k < max_nnz_slice; k++) {

        for (int i = 0; i < rowsPerSlice; i++) {

            if (pos >= slice_end) {
                printf("ERROR overflow slice %d\n", s);
                exit(1);
            }

            int r = start + i;

            if (r < rows) {
                int row_nnz = row_ptr[r+1] - row_ptr[r];

                if (k < row_nnz) {
                    int idx = row_ptr[r] + k;
                    (*values_array)[pos] = val[idx];
                    (*column_indices)[pos] = col_ind[idx];
                } else {
                    (*values_array)[pos] = 0.0;
                    (*column_indices)[pos] = -1;
                }
            } else {
                (*values_array)[pos] = 0.0;
                (*column_indices)[pos] = -1;
            }

            pos++;
        }
    }
}    //DEBUG PRINT

    printf("SELL result:\n");
    for (int i = 0; i < vectorSize; i++) {
        printf("[%d] col=%d val=%f\n", i, (*column_indices)[i], (*values_array)[i]);
    }

	printf("computation of the res : \n");
	int rowActu = 0;
	int columnActu = 0;
	int valuesIndex = 0;
	int start_line = 0, end_line = 0;
	int start = 0;
	int nbElmBlock = 0;
	for (int indexOffset = 1;indexOffset<nbSlices+1;indexOffset++) {
				nbElmBlock = (*slice_offsets)[indexOffset] - (*slice_offsets)[indexOffset-1];

				rowActu = start_line;
				end_line = (start_line + sliceSize - 1);
				
				//in the case the end_line is outside of the matrix due to the slice size
				//if (end_line >= rows) end_line = rows-1;
				//if (rowActu >= rows) rowActu = end_line-1;
				printf("start line --------- : %d \n",start_line);
				for (int nnzBlock = 0;nnzBlock < nbElmBlock; nnzBlock++) {
					printf("valuesIndex : %d \n",valuesIndex);
					if (rowActu > end_line) rowActu = start_line; 
					 if (rowActu <rows && (*column_indices)[valuesIndex] != -1)	{

					printf("line_actu : %d, end_line : %d \n",rowActu,end_line);
						res_array[rowActu] += ones[(*column_indices)[valuesIndex]]*(*values_array)[valuesIndex]; 
						rowActu++;
					}
					else rowActu++;
					
					valuesIndex++;
						
				}
				start_line += sliceSize;
	}

	int i =0; 
	printf("SELL SpmV Res \n");
	for (;i<rows;i++) {

		printf("y[%d] = %f\n",i,res_array[i]);
	}



    //FREE

    
    free(col_ind);
    free(val);
    //free(slice_offsets);
   // free(column_indices);
    //free(values_array);
}



void computeSpmvSELL(int sliceSize,int nnz, int * rows_array,int * cols_array, double * vals_array, int rows, int cols, int * row_ptr,double * ones,double * res_array) {
	int nbSlices = 0; 
	if (sliceSize > rows || sliceSize <= 0) sliceSize = rows; //ELLPACK by default
	nbSlices = (rows + sliceSize - 1) / sliceSize;
	printf("nb slices : %d\n",nbSlices);	
    // number of rows per slice
    int rowsPerSlice = sliceSize;

    //Using  CSR for SELL (sorted) 
	// easier to access to the NNZs per row

    int *col_ind = (int *) calloc(nnz, sizeof(int));
    double *val   = (double *) calloc(nnz, sizeof(double));
    //int *offset = (int *) calloc(rows, sizeof(int));
   // for (int i = 0; i < rows; i++) offset[i] = row_ptr[i];

    for (int i = 0; i < nnz; i++) {
        int r = rows_array[i];
    //    int pos = offset[r]++;
        col_ind[i] = cols_array[i];
        val[i]     = vals_array[i];
    }

   // free(offset);

    //compute slice_offsets 

	int *slice_offsets = (int *) calloc(nbSlices+1, sizeof(int));
	slice_offsets[0] = 0;

	for (int s = 0; s < nbSlices; s++) {

		int start = s * rowsPerSlice;
		int end   = start + rowsPerSlice;
		if (end > rows) end = rows;

		int max_nnz = 0;

		// compute max nnz in slice (lignes réelles seulement)
		for (int r = start; r < end; r++) {
			int nnz_row = row_ptr[r+1] - row_ptr[r];
			if (nnz_row > max_nnz)
				max_nnz = nnz_row;
		}

		slice_offsets[s+1] = slice_offsets[s] + max_nnz * rowsPerSlice;

		printf("slice %d : max_nnz=%d offset=%d\n",
			   s, max_nnz, slice_offsets[s]);
	}
    //allocate SELL 
	
	int vectorSize = slice_offsets[nbSlices]; 
	printf("vector size : %d \n",slice_offsets[nbSlices]);
    int *column_indices = (int *) calloc(vectorSize, sizeof(int));
    double *values_array = (double *) calloc(vectorSize, sizeof(double));
   

    // fill SELL 

    int pos = 0;

  for (int s = 0; s < nbSlices; s++) {

    int start = s * rowsPerSlice;
    int slice_start = slice_offsets[s];
    int slice_end   = slice_offsets[s+1];

    int pos = slice_start;

    int max_nnz_slice = (slice_end - slice_start) / rowsPerSlice;

    for (int k = 0; k < max_nnz_slice; k++) {

        for (int i = 0; i < rowsPerSlice; i++) {

            if (pos >= slice_end) {
                printf("ERROR overflow slice %d\n", s);
                exit(1);
            }

            int r = start + i;

            if (r < rows) {
                int row_nnz = row_ptr[r+1] - row_ptr[r];

                if (k < row_nnz) {
                    int idx = row_ptr[r] + k;
                    values_array[pos] = val[idx];
                    column_indices[pos] = col_ind[idx];
                } else {
                    values_array[pos] = 0.0;
                    column_indices[pos] = -1;
                }
            } else {
                values_array[pos] = 0.0;
                column_indices[pos] = -1;
            }

            pos++;
        }
    }
}    //DEBUG PRINT

    printf("SELL result:\n");
    for (int i = 0; i < vectorSize; i++) {
        printf("[%d] col=%d val=%f\n", i, column_indices[i], values_array[i]);
    }

	printf("computation of the res : \n");
	int rowActu = 0;
	int columnActu = 0;
	int valuesIndex = 0;
	int start_line = 0, end_line = 0;
	int start = 0;
	int nbElmBlock = 0;
	for (int indexOffset = 1;indexOffset<nbSlices+1;indexOffset++) {
				nbElmBlock = slice_offsets[indexOffset] - slice_offsets[indexOffset-1];

				rowActu = start_line;
				end_line = (start_line + sliceSize - 1);
				
				//in the case the end_line is outside of the matrix due to the slice size
				//if (end_line >= rows) end_line = rows-1;
				//if (rowActu >= rows) rowActu = end_line-1;
				printf("start line --------- : %d \n",start_line);
				for (int nnzBlock = 0;nnzBlock < nbElmBlock; nnzBlock++) {
					printf("valuesIndex : %d \n",valuesIndex);
					if (rowActu > end_line) rowActu = start_line; 
					 if (rowActu <rows && column_indices[valuesIndex] != -1)	{

					printf("line_actu : %d, end_line : %d \n",rowActu,end_line);
						res_array[rowActu] += ones[column_indices[valuesIndex]]*values_array[valuesIndex]; 
						rowActu++;
					}
					else rowActu++;
					
					valuesIndex++;
						
				}
				start_line += sliceSize;
	}

	int i =0; 
	printf("SELL SpmV Res \n");
	for (;i<rows;i++) {

		printf("y[%d] = %f\n",i,res_array[i]);
	}



    //FREE

    
    free(col_ind);
    free(val);
    free(slice_offsets);
    free(column_indices);
    free(values_array);
}


int main(int argc, char *argv[]) {
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

	double * csrRes = (double*) calloc(rows,sizeof(double));
	computeSpmvCSR(csrRes,rows_array,cols_array,vals_array,ones,nnz,rows,row_ptr_array);

	//spMV SELL 
	
	///void computeSpmvSELL(int nbSlices, COOvalue * cooArray, int rows, int cols) {
	double * sellRes = (double*) calloc(rows,sizeof(double));

	//computeSpmvSELL(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes);
	
	int sizeSellVect = 0; 
	int sizeSliceOffset = 0;	
	int * slice_offsetsSell;	
	int * column_indicesSell; 
	double * values_arraySell;
	computeSpmvSELLv2(sliceSize,nnz,rows_array,cols_array,vals_array,rows,cols,row_ptr_array,ones,sellRes,&column_indicesSell, &values_arraySell,&slice_offsetsSell,&sizeSellVect,&sizeSliceOffset); 


	printf("sizeSellVect : %d \n",sizeSellVect);
	printf("sizeSliceOffset : %d \n",sizeSliceOffset);

  
    /* EX1: Banchmark the runtime of spmv_COO and spmv_CSR separately. With 2 warmups and 10 iterations, calculate their arithmetic mean
    */

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    /* EX1: Banchmark the effective bandwidths of spmv_COO and spmv_CSR separately.
    */

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */



    /* Note, free the memory */
	free(row_ptr_array);
	free(rows_array);
	free(cols_array);
	free(vals_array);
	free(sellRes);
	//free(cooarray);	
 free(slice_offsetsSell);
    free(column_indicesSell);
    free(values_arraySell);

	free(cooRes); 
	free(csrRes);
	free(ones);
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    return 0;
}

