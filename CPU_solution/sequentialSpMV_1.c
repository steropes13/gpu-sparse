#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

	
// todo: rewrite the functions with cooarrow (coovalues) 

int compare(const void * a, const void * b) {
	COOvalue * ae = (COOvalue *) a; 
	COOvalue * be = (COOvalue *) b; 
	if (ae->row < be->row) return -1; 
	if (ae->row > be->row) return 1;

	if (ae->col < be->col) return -1;
	if (ae->col > be->col) return 1; 
	
	return 0; 
	

}

int main(int argc, char *argv[]) {
	unsigned int seed = 19;
	double min = 0.0, max = 5.0;
    if (argc < 4) {
        printf("Usage: %s <rows> <cols> <nnz>\n", argv[0]);
        return 1;
    }

    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    int nnz = atoi(argv[3]);
	if (argc == 5) {
		seed = atoi(argv[4]);
		srand(seed);
	}
	else srand(time(NULL));
	printf("Value of nnz : %d \n",nnz);
	int occurRow = 0, actuRow = 0, numberOccur = 0;

	double * ones = (double *) malloc(cols*sizeof(double));
	char * matrix = (char*) calloc(rows*cols,sizeof(char)); 
	COOvalue * cooarray = (COOvalue *) calloc(nnz,sizeof(COOvalue)); 
    double timers[NITER];


		int randomRow = 0, randomCol = 0;
		double randomVal = 0.0;


		for (int i=0;i<cols;i++) {ones[i] = 1.0;}


		for (int i = 0;i<nnz;i++) {
				randomVal = random_double(min,max); 
				printf("Value  %d : %.2f\n", i,randomVal); 
				randomRow = rand() % (rows); 
				do {
					randomCol = rand() % (cols); 
					printf("col %d found for row %d \n",randomCol,randomRow);
				} while (matrix[randomRow*cols+randomCol]);
				matrix[randomRow*cols+randomCol] = 1; 
				printf("Row position  : %d\n", randomRow); 
				printf("Col position  : %d\n", randomCol); 
				cooarray[i].row = randomRow; 
				cooarray[i].col = randomCol;		
				cooarray[i].val = randomVal;
}
		free(matrix);
		//c) SpMV COO 
		double * cooRes = (double*) calloc(rows,sizeof(double));
		for (int number = 0 ; number < nnz ; number++) {
				cooRes[cooarray[number].row] += cooarray[number].val*ones[cooarray[number].col];
		}
		printf("COO res :\n");
		for (int i = 0; i < rows; i++) {
    printf("y[%d] = %.2f\n", i, cooRes[i]);
}


	 
	qsort(cooarray,nnz,sizeof(COOvalue),compare);
	int * occurenceArray = calloc(rows, sizeof(int));
	printf("sorting the arrys : \n");
	for (int i = 0; i < nnz;i++) {
		printf("value  : %.2f\n",cooarray[i].val); 
		printf("Row position  : %d\n", cooarray[i].row); 
		printf("Col position  : %d\n", cooarray[i].col); 
}
	
	//number of elements for each row 
	for (int j = 0; j <nnz;j++) {
		occurenceArray[cooarray[j].row] += 1;
}
	
	//prefix sum 
	int * row_ptr = calloc(rows+1,sizeof(int)); 
	row_ptr[0] = 0; 
	for (int j = 1; j <rows+1;j++) {
		row_ptr[j] = row_ptr[j-1] + occurenceArray[j-1];
		printf("prefix sum  %d : %d \n",j, row_ptr[j]);
	}	
	
	printf("final value of prefix-sum : %d \n",row_ptr[rows]);


	for (int j = 0; j <rows;j++) {
		printf("row %d : %d \n",j, occurenceArray[j]);
}
	//csr to matrix + result for ones
	double * csrRes = (double*) calloc(rows,sizeof(double));
	for (int i = 0;i<rows;i++) {
			for (int j = row_ptr[i]; j < row_ptr[i+1];j++) {
					printf("element (%d %d), value : %.2f",i,cooarray[j].col,cooarray[j].val); 
					csrRes[i] += cooarray[j].val*ones[cooarray[j].col];
}
	}

	printf("result csr : \n"); 
	for (int i = 0;i<rows;i++) {
    		printf("y[%d] = %.2f\n", i, csrRes[i]);
}


	



    /* EX1: Check if spmv_COO and spmv_CSR have the same output
    */

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */


  
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
	free(occurenceArray);
	free(row_ptr);
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    return 0;
}

