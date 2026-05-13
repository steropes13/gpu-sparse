#include "../include/my_time_lib.h"


float arithmetic_mean(float *v, int len) {

    float mu = 0.0;
    for (int i=0; i<len; i++)
        mu += (float)v[i];
    mu /= (float)len;

    return(mu);
}

float geometric_mean(float *v, int len) {
    
    float mu = 1.0;
    for (int i=0; i<len; i++) {
        mu *= (v[i] > 0) ? ((float)v[i]) : 1;
    }
    mu = pow(mu, 1.0 / len);
    
    return(mu);
}

float sigma_fn_sol(float *v, float mu, int len) {

    float sigma = 0.0;
    for (int i=0; i<len; i++) {
        sigma += ((float)v[i] - mu)*((float)v[i] - mu);
    }
    sigma /= (float)len;

    return(sigma);
}

// -------------------------------------------------
