#include "../include/spmvCPU.h"

double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

int random_int(int min, int max) {
   	return min + (max - min) * rand();
}


	

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
void computeSpmvSELLv2(int sliceSize,int nnz, int * rows_array,int * cols_array, double * vals_array, int rows, int cols, int * row_ptr,double * ones,double * res_array,int ** column_indices, double ** values_array, int ** slice_offsets, int * sizeVect, int * sizeOffset) {
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

    pos = slice_start;

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
	}

    //DEBUG PRINT

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
                     if (rowActu <rows && (*column_indices)[valuesIndex] != -1) {

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


