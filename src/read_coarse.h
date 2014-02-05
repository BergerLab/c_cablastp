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
/*void read_coarse(struct cbp_coarse *coarse_db, FILE *links_file,
                 struct fasta_seq_gen *fsg);*/

struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end);

int64_t cbp_coarse_link_offset(struct cbp_coarse *coarsedb, int id);

#endif
