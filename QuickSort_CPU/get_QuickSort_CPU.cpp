#include <memory.h>
#include "mex.h"
	
void quickSort(double *arr, int left, int right) {
    int i = left, j = right;
    double tmp;
    double pivot = arr[(left + right)>>1];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
              i++;
        while (arr[j] > pivot)
              j--;
        if (i <= j) {
              tmp = arr[i];
              arr[i] = arr[j];
              arr[j] = tmp;
              i++;
              j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

void mexFunction(int nlhs,  mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *arr_i, *arr_o;
    int m, n;


    arr_i = mxGetPr(prhs[0]);
    m = (int)mxGetM(prhs[0]); 
    n = (int)mxGetN(prhs[0]); 

    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    arr_o = mxGetPr(plhs[0]);
    memcpy(arr_o, arr_i, n*m*sizeof(double));

    quickSort(arr_o, 0, n*m-1);
}