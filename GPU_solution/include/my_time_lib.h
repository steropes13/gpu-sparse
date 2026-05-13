#ifndef LAB1_EX2_LIB
#define LAB1_EX2_LIB

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/time.h>
#include <math.h>

#define STR(s) #s
#define XSTR(s) STR(s)

#define TIMER_DEF     struct timeval temp_1, temp_2

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)+(temp_2.tv_usec-temp_1.tv_usec)/1000000.0)





float geometric_mean(float *v, int len);
float arithmetic_mean(float *v, int len);
float sigma_fn_sol(float *v, float mu, int len);

#ifdef __cplusplus
}
#endif

#endif
