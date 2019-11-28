#ifndef my_func_h
#define my_func_h

#include <ap_int.h>

#define NDATA 12

ap_int<30> myfunc(const ap_int<16> a[NDATA], const ap_int<16> b[NDATA], ap_int<24> c[NDATA]) ;

#endif
