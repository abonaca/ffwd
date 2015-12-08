#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// functions
double compare(double **x_obs, double **x_mod, double **err_obs, double ***sigma_inv, double *logdet, int N, int k, int K);
double median(int n, double x[]);

// void sort(unsigned long n, double arr[]);
// 
// // Numerical Recipes definitions
// #define NRANSI
// #include "nrutil.h"
// #define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
// #define M 7
// #define NSTACK 50