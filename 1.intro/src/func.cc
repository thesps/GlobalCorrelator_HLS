#include "func.h"


ap_int<30> myfunc(const ap_int<16> a[NDATA], const ap_int<16> b[NDATA], ap_int<24> c[NDATA]) {
    // Uncomment the following to tell Vivado to pipeline the function
    //#pragma HLS pipeline II=1

    // Uncomment the following to tell Vivado to fully split a,b,c into registers instead of using BRAMs
    //#pragma HLS array_partition variable=a complete
    //#pragma HLS array_partition variable=b complete
    //#pragma HLS array_partition variable=c complete

    int sum = 0;
    for (unsigned int i = 0; i < NDATA; ++i) {
        int prod = (a[i] * b[i]) >> 8;
        c[i] = prod;
        sum += prod;
    }
    return sum;
}
