#ifndef __CABLASTP_READ_COMPRESSED_H__
#define __CABLASTP_READ_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "coarse.h"
#include "compressed.h"
#include "stdbool.h"

char *get_compressed_header(FILE *f);
struct cbp_link_to_coarse *read_link(FILE *f);
struct cbp_compressed_seq *get_compressed_seq(FILE *f, int id);
struct cbp_compressed_seq **read_compressed(FILE *f);
int64_t cbp_compressed_link_offset(struct cbp_compressed *comdb, int id);
struct cbp_seq *cbp_compressed_read_seq(struct cbp_compressed *com_db,
                                        struct cbp_coarse *coarse_db, int id);
#endif
