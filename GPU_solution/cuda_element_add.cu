#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda_runtime.h>

#define TIMER_DEF     struct timeval temp_1, temp_2

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)+(temp_2.tv_usec-temp_1.tv_usec)/1000000.0)


#define STR(s) #s
#define XSTR(s) STR(s)
#define dtype float

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

__global__ 
void kernel_element_add(int n, dtype *a, dtype* b, dtype* c)
{
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

	//each thread has to cover a part ot the array 
	//one block, multiple thread and covers all n
	for (int i = threadIdx.x;i<n;i+=blockDim.x) {
		c[i] = a[i] + b[i]; 
	}

}


int main(int argc, char *argv[])
{
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



    printf("====================================== Problem computations ======================================\n");
    // =========================================== Set-up the problem ============================================

    if (argc < 2) {
        printf("Usage: lab6_ex3 n\n");
        return(1);
    }
    printf("argv[1] = %s\n", argv[1]);

    // ---------------- set-up the problem size -------------------

    int n = atoi(argv[1]), len = (1<<n), i;

    printf("n = %d --> len = 2^(n) = %d\n", n, len);
    printf("dtype = %s\n", XSTR(dtype));

    // ------------------ set-up the timers ---------------------

    TIMER_DEF;
    float error, cputime, gputime_host, gputime_event;

    // ------------------- set-up the problem -------------------

    dtype *a, *b, *CPU_c, *GPU_c;
    a = (dtype*)malloc(sizeof(dtype)*len);
    b = (dtype*)malloc(sizeof(dtype)*len);
    CPU_c = (dtype*)malloc(sizeof(dtype)*len);
    GPU_c = (dtype*)malloc(sizeof(dtype)*len);
    time_t t;
    srand((unsigned) time(&t));
    int blk_size = 256;
    int grd_size = 1; //needs to be changed (train with other values)


    // ------------------------gererate random float -----------------------
    	 
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

    for (int truc = 0; truc < len;truc++) {
	a[truc] = (float)(rand()) / RAND_MAX; 
	b[truc] = (float)(rand()) / RAND_MAX; 
    }
    
    // ======================================== Running the computations on CPU =========================================
    
    TIMER_START;
    for (i=0; i<len; i++)
        CPU_c[i] = a[i] + b[i];
    TIMER_STOP;
    //   errors[0] = 0.0;
    cputime = TIMER_ELAPSED;

   

    // ---------------- allocing GPU vectors -------------------
    dtype *dev_a, *dev_b, *dev_c;

    cudaMalloc(&dev_a, len*sizeof(dtype));
    cudaMalloc(&dev_b, len*sizeof(dtype));
    cudaMalloc(&dev_c, len*sizeof(dtype));


    // ------------ copy date from host to device --------------

    cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, len*sizeof(dtype), cudaMemcpyHostToDevice);
    cudaMemset(dev_c, 0, len*sizeof(dtype));



    // ------------ Running the computation on GPU, measure time on host side --------------
    
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */
    TIMER_START;
    kernel_element_add<<<grd_size,blk_size>>>(len,dev_a,dev_b,dev_c);
    cudaDeviceSynchronize();
    TIMER_STOP; 
    gputime_host = TIMER_ELAPSED; 

    // ----------- copy results from device to host ------------

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */


    cudaMemcpy(GPU_c, dev_c, len*sizeof(dtype),cudaMemcpyDeviceToHost);

    // ------------- Compare GPU and CPU solution --------------

    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */
    // float error, cputime, gputime_host, gputime_event;
	
   // error = (float) fabs(cputime - gputime_host); 
    //compare the real error 
    error = 0.0;
    for (int i = 0; i < len; i++) {
    		error += fabs(CPU_c[i] - GPU_c[i]);
	}



    // ------------------- reset gpu buffers ---------------------
    cudaMemset(dev_a, 0, len*sizeof(dtype));
    cudaMemset(dev_b, 0, len*sizeof(dtype));
    cudaMemcpy(dev_a, a, len*sizeof(dtype), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, len*sizeof(dtype), cudaMemcpyHostToDevice);
    cudaMemset(dev_c, 0, len*sizeof(dtype));
    
    // ------------------- measure time with cuda event -------------------
   
    /* |========================================| */
    /* |           Put here your code           | */
    /* |========================================| */

cudaEvent_t start, stop; 
cudaEventCreate(&start); 
cudaEventCreate(&stop); 
cudaEventRecord(start); 
kernel_element_add<<<grd_size,blk_size>>>(len,dev_a,dev_b,dev_c);
cudaEventRecord(stop);
cudaEventSynchronize(stop); 
cudaEventElapsedTime(&gputime_event,start,stop);

cudaEventDestroy(start);
cudaEventDestroy(stop);


    // ============================================ Print the results ============================================
    printf("================================== Times and results of my code ==================================\n");
    printf("Error between CPU and GPU is %5.3f\n", error);
    printf("\nVector len = %d, CPU time = %5.3f\n", len, cputime);
    printf("\nblk_size = %d, grd_size = %d, GPU time (gettimeofday): %5.3f sec\n", blk_size, grd_size, gputime_host);
    printf("\nblk_size = %d, grd_size = %d, GPU time (CUDA Event): %5.3f ms\n", blk_size, grd_size,gputime_event);
 
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }


    // ----------------- free GPU variable ---------------------

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
    free(a);
    free(b);
    free(CPU_c);
    free(GPU_c);

    // ---------------------------------------------------------
    return 0;
}

