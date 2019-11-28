#include "src/func.h"
#include <cstdio>
#include <cstdlib>

int main() {
    srand(125);
    ap_int<16> a[NDATA], b[NDATA]; ap_int<24> c[NDATA]; ap_int<30> fret;
    for (unsigned int itest = 0, ntest = 20; itest <= ntest; ++itest) {
        // create some input data
        int check_sum = 0;
        for (unsigned int i = 0; i < NDATA; ++i) {
            a[i] = ap_int<16>(rand() & 0xFFFF);
            b[i] = ap_int<16>(rand() & 0xFFFF);
            check_sum += (a[i] * b[i]) >> 8;
        }
        // call the function
        fret = myfunc(a,b,c);
        // check the results (here we check only the total, for lazyness)
        if (fret != check_sum) {
            printf("ERROR: test %d: sum = %d , while %d was expected\n", itest, int(fret), check_sum); 
            return 1;
        } else {
            printf("INFO: test %d: sum = %d , expected %d ---> ok \n", itest, int(fret), check_sum); 
        }
    }
}
