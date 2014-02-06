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

#endif