#ifndef __CABLASTP_COMPRESSED_H__
#define __CABLASTP_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "align.h"
#include "bitpack.h"
#include "edit_scripts.h"

#include "stdbool.h"

struct cbp_link_to_coarse {
    char *diff;
    uint16_t coarse_seq_id;
    uint64_t original_start;
    uint64_t original_end;
    uint16_t coarse_start;
    uint16_t coarse_end;
    struct cbp_link_to_coarse *next;
};

struct cbp_link_to_coarse *
cbp_link_to_coarse_init(int32_t coarse_seq_id,
                        uint64_t original_start, uint64_t original_end,
                        uint16_t coarse_start, uint16_t coarse_end,
                        struct cbp_alignment alignment, bool dir);

struct cbp_link_to_coarse *
cbp_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               uint64_t original_start, uint64_t original_end,
                               uint16_t coarse_start, uint16_t coarse_end,
                               bool dir);

void
cbp_link_to_coarse_free(struct cbp_link_to_coarse *link);

struct cbp_compressed_seq {
    uint64_t id;
    char *name;
    struct cbp_link_to_coarse *links;
};

struct cbp_compressed_seq *
cbp_compressed_seq_init(int32_t id, char *name);

void
cbp_compressed_seq_free(struct cbp_compressed_seq *seq);

void
cbp_compressed_seq_addlink(struct cbp_compressed_seq *seq,
                           struct cbp_link_to_coarse *link);

struct cbp_compressed {
    struct DSVector *seqs;
    FILE *file_compressed;
    FILE *file_index;
};

struct cbp_compressed *
cbp_compressed_init(FILE *file_compressed, FILE *file_index);

void
cbp_compressed_free(struct cbp_compressed *com_db);

int32_t
cbp_compressed_size(struct cbp_compressed *com_db);

void
cbp_compressed_add(struct cbp_compressed *com_db,
                   struct cbp_compressed_seq *seq);

void
cbp_compressed_save_binary(struct cbp_compressed *com_db);

void
cbp_compressed_save_plain(struct cbp_compressed *com_db);

void
cbp_compressed_write(struct cbp_compressed *com_db,
                     struct cbp_compressed_seq *seq);

void
cbp_compressed_write_binary(struct cbp_compressed *com_db,
                            struct cbp_compressed_seq *seq);

static struct cbp_compressed_seq *
cbp_compressed_seq_at(struct cbp_compressed *com_db, int32_t i);
#endif
