#ifndef __CABLASTP_READ_COMPRESSED_H__
#define __CABLASTP_READ_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "compressed.h"
#include "stdbool.h"

char *get_compressed_header(FILE *f);
struct cbp_link_to_coarse *read_link(FILE *f);
struct cbp_compressed_seq *get_compressed_sequence(FILE *f);
struct cbp_compressed_seq **read_compressed(FILE *f);

#endif
