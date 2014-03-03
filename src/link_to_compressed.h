#ifndef __CABLASTP_LINK_TO_COMPRESSED_H__
#define __CABLASTP_LINK_TO_COMPRESSED_H__

#include <stdbool.h>
#include <stdint.h>

struct cbp_link_to_compressed {
    bool dir;
    int32_t org_seq_id;
    int16_t coarse_start;
    int16_t coarse_end;
    uint64_t original_start;
    uint64_t original_end;
    struct cbp_link_to_compressed *next;
};
#endif
