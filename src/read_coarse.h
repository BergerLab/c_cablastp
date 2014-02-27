#ifndef __CABLASTP_READ_COARSE_H__
#define __CABLASTP_READ_COARSE_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "coarse.h"
#include "compressed.h"
#include "fasta.h"

#include "stdbool.h"

char *get_coarse_header(FILE *f);
struct cbp_link_to_compressed *read_coarse_link(FILE *f);
struct DSVector *get_coarse_sequence_links(FILE *f);
struct DSVector *get_coarse_sequence_links_at(FILE *links, FILE *index,
                                                           int32_t id);


struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length);

int64_t cbp_coarse_find_offset(FILE *index_file, int id);
struct fasta_seq *cbp_coarse_read_fasta_seq(struct cbp_coarse *coarsedb,
                                            int id);
#endif
