#ifndef __CABLASTP_READ_COARSE_H__
#define __CABLASTP_READ_COARSE_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "coarse.h"
#include "stdbool.h"

char *get_coarse_header(FILE *f);
struct cbp_link_to_compressed *read_coarse_link(FILE *f);
struct cbp_coarse_seq **read_coarse(FILE *f);

#endif
