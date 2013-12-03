#ifndef __CABLASTP_ALIGN_H__
#define __CABLASTP_ALIGN_H__

/*
This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the biogo
package: code.google.com/p/biogo.

It's mostly copied from its original form, but it is optimized specifically
for cablastp to limit allocations and to absolve the need for the biogo/seq.Seq
type.
*/

#include <stdint.h>

#define CABLASTP_ALIGN_SEQ_SIZE 10000

int32_t
cbp_align_ungapped(char *rseq, int32_t rstart, int32_t rend, int32_t dir1, int32_t i1,
                   char *oseq, int32_t ostart, int32_t oend, int32_t dir2, int32_t i2,
                   bool *matches, bool *matches_past_clump, int matches_index);

int32_t
cbp_align_identity(char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend);

struct cbp_align_nw_memory {
    int32_t *table;
    int32_t *zeroes;
    char *ref;
    char *org;
};

typedef struct cbp_nw_tables{
    int **dp_score;
    int **dp_from;
};

struct cbp_nw_tables
make_nw_tables(char *rseq, int dp_len1, int i1, int dir1,
               char *oseq, int dp_len2, int i2, int dir2);


struct cbp_align_nw_memory *
cbp_align_nw_memory_init();

void
cbp_align_nw_memory_free(struct cbp_align_nw_memory *mem);

struct cbp_alignment {
    char *ref;
    char *org;
    int32_t length;
};

struct cbp_alignment
cbp_align_nw(struct cbp_align_nw_memory *mem,
             char *rseq, int dp_len1, int i1, int dir1,
             char *oseq, int dp_len2 ,int i2, int dir2);

int32_t
cbp_align_length_nogaps(char *residues);

int32_t
attempt_ext(int32_t i1, const int32_t dir1, const char *s1, int32_t len1,
            int32_t start1, int32_t i2, const int32_t dir2, const char *s2,
            int32_t len2, int32_t start2);

int check_and_update(bool *matches, int *matches_index, int *num_matches, bool *temp, int temp_index);

int max_dp_len(int i, int dir, int len);
#endif
