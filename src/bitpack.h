#include <stdint.h>

#ifndef __CABLASTP_EDITSCRIPTS_H__
#define __CABLASTP_EDITSCRIPTS_H__

uint64_t make_mask(int bytes);
uint64_t shift_left(uint64_t x, int bits);
uint64_t shift_right(uint64_t x, int bits);
char *read_int_to_bytes(uint64_t number, int bytes);

#endif
