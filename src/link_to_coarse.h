#ifndef __CABLASTP_LINK_TO_COARSE_H__
#define __CABLASTP_LINK_TO_COARSE_H__

#include <stdint.h>

struct cbp_link_to_coarse {
    char *diff;
    uint16_t coarse_seq_id;
    uint64_t original_start;
    uint64_t original_end;
    uint16_t coarse_start;
    uint16_t coarse_end;
    struct cbp_link_to_coarse *next;
};
#endif
