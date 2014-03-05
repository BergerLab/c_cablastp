#ifndef __CABLASTP_READ_COARSE_H__
#define __CABLASTP_READ_COARSE_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "coarse.h"
#include "compressed.h"
#include "fasta.h"

#include "stdbool.h"

struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length);
#endif
