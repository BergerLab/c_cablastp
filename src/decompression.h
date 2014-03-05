#ifndef __CABLASTP_DECOMPRESSION_H__
#define __CABLASTP_DECOMPRESSION_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "stdbool.h"

#include "coarse.h"
#include "compressed.h"
#include "seq.h"

struct cbp_seq *cbp_decompress_seq(struct cbp_compressed_seq *cseq,
                                   struct cbp_coarse *coarsedb);
struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length);

#endif
