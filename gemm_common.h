#ifndef _GEMM_COMMON_H
#define _GEMM_COMMON_H

#define FLOAT float
//typedef float FLOAT;
typedef long BLASLONG;
typedef unsigned int BLASUINT;

/* or gcc build in function: int  __builtin_popcount(unsigned int) */
inline int popcount( BLASUINT i )
{
    /*
     * http://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer?page=1&tab=votes#tab-top
     * C or C++: use uint32_t
     */
    i	= i - ( (i >> 1) & 0x55555555);
    i	= (i & 0x33333333) + ( (i >> 2) & 0x33333333);
    return( ( ( (i + (i >> 4) ) & 0x0F0F0F0F) * 0x01010101) >> 24);
}

#endif