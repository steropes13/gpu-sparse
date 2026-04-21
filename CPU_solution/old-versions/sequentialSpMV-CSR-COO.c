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

