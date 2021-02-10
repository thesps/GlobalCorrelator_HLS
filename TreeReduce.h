#ifndef TREEREDUCE_H_
#define TREEREDUCE_H_

/* ---
 * Tree Reduce implementation borrowed from hls4ml
 * --- */

#include "ap_fixed.h"

constexpr int ceillog2(int x){
    return (x <= 2) ? 1 : 1 + ceillog2((x+1) / 2);
}

constexpr int floorlog2(int x){
    return (x < 2) ? 0 : 1 + floorlog2(x / 2);
}

template<int B>
constexpr int pow(int x){
    return x == 0 ? 1 : B * pow<B>(x - 1);
}

constexpr int pow2(int x){
    return pow<2>(x);
}

/*constexpr int pow2(int x){
    return x == 0 ? 1 : 2 * pow2(x - 1);
}*/

/* ---
 * Balanced tree reduce implementation.
 * For use in scenarios where Vivado cannot expression balance
 * Reduces an array of inputs to a single value using the template binary operator 'Op',
 * for example summing all elements with Op_add, or finding the maximum with Op_max
 * Use only when the input array is fully unrolled. Or, slice out a fully unrolled section
 * before applying and accumulate the result over the rolled dimension.
 * --- */
template<class T, int N, class Op>
T reduce(const T* x, Op op){
    #pragma HLS pipeline II=1
    static constexpr int leftN = pow2(floorlog2(N - 1)) > 0 ? pow2(floorlog2(N - 1)) : 0;
    static constexpr int rightN = N - leftN > 0 ? N - leftN : 0;
    if(N == 1){
        return x[0];
    }else if(N == 2){
        return op(x[0],x[1]);
    }else{
        T left[leftN];
        T right[rightN];
        #pragma HLS array_partition variable=left complete
        #pragma HLS array_partition variable=right complete
        for(int i = 0; i < leftN; i++){
            #pragma HLS unroll
            left[i] = x[i];
        }
        for(int i = 0; i < rightN; i++){
            #pragma HLS unroll
            right[i] = x[i+leftN];
        }
        return op(reduce<T,leftN,Op>(left, op), reduce<T,rightN,Op>(right, op));
    }
}

template<class T>
class Op_add{
    public:
    T operator()(T a, T b){
        return a + b;
    }
};

template<class T>
class Op_and{
    public:
    T operator()(T a, T b){
        return a && b;
    }
};

template<class T>
class Op_or{
    public:
    T operator()(T a, T b){
        return a || b;
    }
};

template<class T>
class Op_max{
    public:
    T operator()(T a, T b){
        return a >= b ? a : b;
    }
};

template<class T>
class Op_min{
    public:
    T operator()(T a, T b){
        return a <= b ? a : b;
    }
};

#endif
