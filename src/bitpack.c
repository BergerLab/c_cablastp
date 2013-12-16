#include <stdlib.h>
#include "bitpack.h"

/*Takes in as input a number of bits and returns a 64-bit mask with that
  number of bits on the right of the mask set to 1 and all bits to the left
  set to 0.*/
uint64_t make_mask(int bits){
    uint64_t mask = (uint64_t)1;

    mask <<= bits;
    mask = mask & (mask-1);
}

/*Left shift for 64-bit integers*/
uint64_t shift_left(uint64_t x, int bits){
    return x<<bits;
}

/*Right shift for 64-bit integers*/
uint64_t shift_right(uint64_t x, int bits){
    int mask = make_mask(64-bits);
    x >>= bits;
    return x & mask;
}
