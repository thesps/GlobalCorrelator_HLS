#include "func.h"


ap_int<30> myfunc(const ap_int<16> a[NDATA], const ap_int<16> b[NDATA], ap_int<24> c[NDATA]) {
    #pragma HLS array_partition variable=a complete
    #pragma HLS array_partition variable=b complete
    #pragma HLS array_partition variable=c complete
    #pragma HLS pipeline II=1

    int sum = 0;
    for (unsigned int i = 0; i < NDATA; ++i) {
        int prod = (a[i] * b[i]) >> 8;
        c[i] = prod;
        sum += prod;
    }
    return sum;
}
