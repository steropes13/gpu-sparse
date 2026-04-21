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

void computeSpmvCOO(double * res, COOvalue *  cooArray, double * vect, int nnz, int rows) {
	for (int number = 0 ; number < nnz ; number++) {
		res[cooArray[number].row] += cooArray[number].val*vect[cooArray[number].col];
	}

	printf("COO res :\n");
	for (int i = 0; i < rows; i++) {
    	printf("y[%d] = %.2f\n", i, res[i]);
		}
	}


void computeSpmvCSR(double * res, COOvalue *  cooArraySorted, double * vect, int nnz, int rows) {
	int * occurenceArray = calloc(rows, sizeof(int));
	printf("sorting the arrys : \n");
	for (int i = 0; i < nnz;i++) {
		printf("value  : %.2f\n",cooArraySorted[i].val); 
		printf("Row position  : %d\n", cooArraySorted[i].row); 
		printf("Col position  : %d\n", cooArraySorted[i].col); 
	}
	
	//number of elements for each row 
	for (int j = 0; j <nnz;j++) {
		occurenceArray[cooArraySorted[j].row] += 1;
	}

	//printing the occurence array
	for (int j = 0; j <rows;j++) {
		printf("row %d : %d \n",j, occurenceArray[j]);
	}

	
	//prefix sum 
	int * row_ptr = calloc(rows+1,sizeof(int)); 
	row_ptr[0] = 0; 
	for (int j = 1; j <rows+1;j++) {
		row_ptr[j] = row_ptr[j-1] + occurenceArray[j-1];
		printf("prefix sum  %d : %d \n",j, row_ptr[j]);
	}	
	free(occurenceArray);	
	printf("final value of prefix-sum : %d \n",row_ptr[rows]);

	//csr to matrix + result for ones
	for (int i = 0;i<rows;i++) {
		for (int j = row_ptr[i]; j < row_ptr[i+1];j++) {
					printf("element (%d %d), value : %.2f\n",i,cooArraySorted[j].col,cooArraySorted[j].val); 
			res[i] += cooArraySorted[j].val*vect[cooArraySorted[j].col];
		}
	}

	printf("result csr : \n"); 
	for (int i = 0;i<rows;i++) {
    	printf("y[%d] = %.2f\n", i, res[i]);
	}

	free(row_ptr);
}

void computeSpmvSELL(int nbSlices,int nnz, COOvalue * cooArray, int rows, int cols) {

	//number of rows per slice computed 
	int rowsPerSlice = !(rows%2)?(rows/nbSlices):((rows+1)/nbSlices); 
	

	//construct the matrix by the cooArray 
	double * matrix;
	//in the case we have an even number of rows
	if(rows%2 == 1) rows++;	
	matrix = (double*) calloc((rows)*cols,sizeof(double));
	
	int i = 0, j = 0;
	int rowCurrent = 0, colCurrent = 0;
	int nnzLineMaxi = 0;

	//printf("nb rows %d \n",rows);

	int * slice_offsets = (int *) calloc(nbSlices+1,sizeof(int));
	int * column_indices;
	double * values_array;
	for (;i<nnz;i++) {
		rowCurrent = cooArray[i].row;
		colCurrent = cooArray[i].col; 
		matrix[rowCurrent*cols+colCurrent] = cooArray[i].val;
		printf("element (%d, %d) : %.2f\n",rowCurrent,colCurrent,matrix[rowCurrent*cols+colCurrent]); 
	}
	


		
	//compute the number of block 

	printf("rowsPerSlice : %d \n", rowsPerSlice);



	//create the different slices 
	//int nnzLineMaxi = 0;
	int newBlock = 0;
	float curElm = 0;
	int counterNnz = 0;
	int currentBlock = 1;
	slice_offsets[0] = 0; 
	
	//compute the slice offsets array	
	for (i = 0;i<rows;i++) {
			//printf("row actu : %d \n",i);
			//printf("nnzLinemaxi %d \n",nnzLineMaxi);
			

			counterNnz = 0;
			for (j = 0;j<cols;j++) {
				curElm = matrix[i*cols+j];
				if (curElm != 0.0) {
						counterNnz++; 
						if (counterNnz > nnzLineMaxi){
						//	printf("counternnz sup : %d (%d %d) \n",counterNnz,i,j);
							nnzLineMaxi = counterNnz;	
						}
							
					}
			}

			if (i != 0 && (i)%rowsPerSlice == 1) {
				slice_offsets[currentBlock] = nnzLineMaxi*rowsPerSlice + slice_offsets[currentBlock-1];	
			printf("bloc %d : %d (line %d), max -> %d \n",currentBlock,slice_offsets[currentBlock],i,nnzLineMaxi);
				newBlock = 1;
				nnzLineMaxi = 0;
				currentBlock++;
				
			}
			
	}

	int vectorSize = slice_offsets[nbSlices];

	column_indices = (int*) calloc(vectorSize,sizeof(int));
	values_array = (double*) calloc(vectorSize,sizeof(double));
	printf("vectorSize of column and values arrays : %d \n", vectorSize);
	int currSellIndx = 0;	
	currentBlock = 1;
	int blockSize = 0;
	nnzLineMaxi = 0;
	int nbNnz = 0;
	int firstNnzLine = 0;
	int startLine = 0,endLine = 0;
	int nbelm = 0;
	int nnzFound = 0;
	//todo : remplacer par une boucle while
	while (currSellIndx<vectorSize) {
		blockSize = slice_offsets[currentBlock] - slice_offsets[currentBlock-1];
		nnzLineMaxi = blockSize/rowsPerSlice;
		nbelm = 0;
		endLine = startLine+rowsPerSlice;
		nnzFound = 0;
		nbNnz = 0;
		
		firstNnzLine = endLine;
		printf("end line : %d : \n",endLine);
		for (int col = 0;col<cols && nbelm < blockSize;col++) {
	
			for (startLine=endLine-rowsPerSlice;startLine<endLine && nbelm < blockSize;startLine++) {

				curElm = matrix[startLine*cols+col];
				if (curElm != 0.0) {
					if (!nnzFound) {
						nnzFound = 1; 
						firstNnzLine = startLine;
					}	
					column_indices[currSellIndx] = col;
					nbNnz++;

					values_array[currSellIndx++] = curElm;
					nbelm++;
				}
				else {
			
					//todo ? rajouter nbNnz >= nnzLineMaxi	
					if (nnzFound && firstNnzLine != startLine ){
						 nbelm++;
						 column_indices[currSellIndx] = -1.0;
						values_array[currSellIndx++] = curElm;
						}

				}
				//todo ? rajouter une autre condition (voir le if juste au dessus) 
			}
		}
		startLine = endLine;
		

		currentBlock++;
		//currSellIndx++;

	}

	printf("tableau final : \n"); 
	for (int i = 0;i<vectorSize;i++) {
		printf("column_indices[%d] = %d \n",i,column_indices[i]);
		
		printf("values_array[%d] = %.2f \n",i,values_array[i]);
		puts("");
	}


	


	
	



	free(slice_offsets);	

	free(column_indices);
	free(values_array);
	
	free(matrix);	
} 



int main(int argc, char *argv[]) {
	//only in case the user wants to read a file .mtx
	int * row_ptr; 
	int * cols_array; 
	double * vals_array;
	char filename[256];
	int user_matrix = 0;
	MM_typecode matcode;
	int return_val = 0;
	

	unsigned int seed = 19;
	int rows = 0; 
	int cols = 0;
	int nnz = 0;
	double min = 0.0, max = 5.0;
	int occurRow = 0, actuRow = 0, numberOccur = 0;
	double * ones;
	char * matrix; 
	COOvalue * cooarray; 
    double timers[NITER];


	int randomRow = 0, randomCol = 0;
	double randomVal = 0.0;




	// so it can be a filename (.mtx)
	if (argc == 2) {
		sprintf(filename, "%s",argv[1]);
		printf("filename : %s \n",filename);
		user_matrix = 1;
		return_val = mm_read_mtx_crd(filename,&rows,&cols,&nnz,&row_ptr,&cols_array,&vals_array,&matcode);
		printf("does the function work ? %d \n",return_val);
		cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 
		for (int i = 0;i<nnz;i++) {
			printf("value %d (%d, %d) : %.2f\n",i,row_ptr[i],cols_array[i],vals_array[i]);
			cooarray[i].col = cols_array[i]-1; 
			cooarray[i].val = vals_array[i];
			cooarray[i].row = row_ptr[i]-1;
			
		}
		if (!return_val) { 
			if (cols_array) free(cols_array);
			if (vals_array) free(vals_array);
			if(row_ptr) free(row_ptr);
		}
		

	}
	//in case the matrix has to be randomized with sizes written by user.
	else {
		if (argc < 4) {
			printf("Usage: %s <rows> <cols> <nnz> or %s <path-of-file.mtx>\n", argv[0],argv[0]);
			return 1;
		}
		rows = atoi(argv[1]);
		cols = atoi(argv[2]);
		nnz = atoi(argv[3]);
		//the case the user added a fifth parameter 
		// and it is a number (seed)
		if (argc == 5) {
			seed = atoi(argv[4]);
			srand(seed);
		}
		else srand(time(NULL));
		printf("Value of nnz : %d \n",nnz);
		
		matrix = (char*) calloc(rows*cols,sizeof(char)); 
		cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 



		for (int i = 0;i<nnz;i++) {
			randomVal = random_double(min,max); 
			printf("Value  %d : %.2f\n", i,randomVal); 
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
		}
		free(matrix);
	}
	//vector we use for the multiplication
	ones = (double *) malloc(cols*sizeof(double));
	for (int i=0;i<cols;i++) {ones[i] = 1.0;}

	//c) SpMV COO 
	double * cooRes = (double*) calloc(rows,sizeof(double));
	computeSpmvCOO(cooRes,cooarray,ones,nnz,rows);

	// spMV CSR  
	qsort(cooarray,nnz,sizeof(COOvalue),compare);

	double * csrRes = (double*) calloc(rows,sizeof(double));
	computeSpmvCSR(csrRes,cooarray,ones,nnz,rows);

	//spMV SELL 
	
	//void computeSpmvSELL(int nbSlices, COOvalue * cooArray, int rows, int cols) {
	computeSpmvSELL(3,nnz,cooarray,rows,cols);




  
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
	free(cooarray);	
	free(cooRes); 
	free(csrRes);
	free(ones);
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    return 0;
}

