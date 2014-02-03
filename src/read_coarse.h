#ifndef __CABLASTP_READ_COARSE_H__
#define __CABLASTP_READ_COARSE_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "coarse.h"
#include "fasta.h"

#include "stdbool.h"

char *get_coarse_header(FILE *f);
struct cbp_link_to_compressed *read_coarse_link(FILE *f);
struct DSVector *get_sequence_links(FILE *f);
/*void read_coarse(struct cbp_coarse *coarse_db, FILE *links_file,
                 struct fasta_seq_gen *fsg);*/

#endif
