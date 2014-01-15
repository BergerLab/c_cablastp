#include <stdint.h>

#ifndef __CABLASTP_BITPACK_H__
#define __CABLASTP_BITPACK_H__

uint64_t make_mask(int bytes);
uint64_t shift_left(uint64_t x, int bits);
uint64_t shift_right(uint64_t x, int bits);
char *read_int_to_bytes(uint64_t number, int bytes);
void output_int_to_file(uint64_t number, int length, FILE *f);

#endif
