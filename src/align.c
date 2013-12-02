#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DNAutils.h"
#include "align.h"
#include "blosum62.h"
#include "flags.h"

int32_t
cbp_align_ungapped(char *rseq, int32_t rstart, int32_t rend, int32_t dir1, int32_t i1,
                   char *oseq, int32_t ostart, int32_t oend, int32_t dir2, int32_t i2,
                   bool *matches, bool *matches_past_clump, int matches_index)
{
    int32_t length, scanned, successive;
    int32_t rlen, olen;
    int32_t id;
    int32_t i, rs, os;
    int32_t consec_matches;
    int32_t matches_since_last_consec;
    int consec_match_clump_size;
    int32_t dir_prod;
    int matches_count;
    int temp_index;

    length = 0;
    scanned = 0;
    consec_match_clump_size = compress_flags.consec_match_clump_size;
    successive = consec_match_clump_size;
    rlen = rend - rstart;
    olen = oend - ostart;
    dir_prod = dir1 * dir2;
    temp_index = 0;
    matches_count = 0;
    matches_since_last_consec = 0;
    for(i = matches_index - 100; i < matches_count; i++)
        if(matches[i])
            matches_count++;

    while(i1 >= rstart && i1 < rend && i2 >= ostart && i2 < oend){
        int cur_ismatch;
        i1 += dir1;
        i2 += dir2;
        scanned++;
        cur_ismatch = bases_match(rseq[i1], oseq[i2], dir_prod);
        if(cur_ismatch == 1){
            matches_past_clump[temp_index] = true;
            temp_index++;
            successive++;
            if(successive >= consec_match_clump_size){
                int update = check_and_update(matches, &matches_index, 
                                              &matches_count,
                                              matches_past_clump, temp_index);
                    length += update;
                    if(update != temp_index)
                        return -1 * length;
                matches_since_last_consec = 0;
            }
            else
                matches_since_last_consec++;
        }
        else {
            successive = 0;
            if(scanned - length >= compress_flags.btwn_match_min_dist_check)
                if((double)matches_since_last_consec <
                   (scanned-length)*compress_flags.btwn_match_ident_thresh)
                    return length;
        }
    }
    return scanned;
}

int32_t
cbp_align_identity(char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend)
{
    int32_t rlen, olen, same, i;

    rlen = rend - rstart;
    olen = oend - ostart;
    assert(rlen == olen);

    if (rlen == 0 && olen == 0)
        return 0;

    same = 0;
    for (i = 0; i < rlen; i++)
        if (rseq[rstart + i] == oseq[ostart + i])
            same++;

    return (same * 100) / rlen;
}

struct cbp_align_nw_memory *
cbp_align_nw_memory_init()
{
    struct cbp_align_nw_memory *mem;
    int seq_size = CABLASTP_ALIGN_SEQ_SIZE;

    mem = malloc(sizeof(*mem));
    assert(mem);

    mem->table = malloc(seq_size * seq_size * sizeof(*mem->table));
    assert(mem->table);

    mem->zeroes = malloc(seq_size * seq_size * sizeof(*mem->zeroes));
    assert(mem->zeroes);

    mem->ref = malloc(seq_size * sizeof(*mem->ref));
    assert(mem->ref);

    mem->org = malloc(seq_size * sizeof(*mem->org));
    assert(mem->org);

    return mem;
}

void
cbp_align_nw_memory_free(struct cbp_align_nw_memory *mem)
{
    free(mem->zeroes);
    free(mem->table);
    free(mem->ref);
    free(mem->org);
    free(mem);
}

struct cbp_alignment
cbp_align_nw(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend)
{
    struct cbp_alignment alignment;
    int32_t *table;
    int32_t gapi;
    int32_t i, i2, i3, j, r, c, off, tablen;
    int32_t sdiag, sup, sleft, rval, oval;
    int32_t constraint;
    bool constrained;
    bool allocated = false;
    char tmp;

    gapi = BLOSUM62_SIZE - 1;
    r = rend - rstart + 1; /* include gap penalty */
    c = oend - ostart + 1; /* include gap penalty */
    tablen = r * c;
    off = 0;

    constrained = true;
    constraint = r / 4;
    if (r <= 11 || c <= 11)
        constrained = false;

    if (tablen > CABLASTP_ALIGN_SEQ_SIZE * CABLASTP_ALIGN_SEQ_SIZE) {
        table = malloc(tablen * sizeof(*table));
        assert(table);
        allocated = true;
    } else {
        table = mem->table;
        for (i = 0; i < tablen; i++)
            table[i] = 0;

        /* for (i = 0; i < r; i++) */
            /* for (j = 0; j < c; j++) { */
                /* We need to be slightly more relaxed on our contraint
                 * here as opposed to below. Namely, we need to make sure to
                 * zero out the boundary cells, which could still be used. */
                /* if (constrained && */
                    /* ((i-j-1) > constraint || (j-i-1) > constraint)) */
                    /* continue; */
                /* table[i * c + j] = 0; */
            /* } */
    }

    for (i = 1; i < r; i++) {
        i2 = (i - 1) * c;
        i3 = i * c;
        for (j = 1; j < c; j++) {
            if (constrained && ((i-j) > constraint || (j-i) > constraint))
                continue;

            rval = BLOSUM62_RESIDUE_TO_INDEX[rseq[rstart + i - 1] - 'A'];
            oval = BLOSUM62_RESIDUE_TO_INDEX[oseq[ostart + j - 1] - 'A'];

            off = i2 + (j - 1);
            sdiag = table[off] + BLOSUM62_MATRIX[rval][oval];
            sup = table[off+1] + BLOSUM62_MATRIX[rval][gapi];
            sleft = table[off+c] + BLOSUM62_MATRIX[gapi][oval];
            
            if (sdiag > sup && sdiag > sleft)
                table[i3+j] = sdiag;
            else if (sup > sleft)
                table[i3+j] = sup;
            else
                table[i3+j] = sleft;
        }
    }

    alignment.ref = mem->ref;
    alignment.org = mem->org;
    alignment.length = 0;

    i = r - 1;
    j = c - 1;
    while (i > 0 && j > 0) {
        rval = BLOSUM62_RESIDUE_TO_INDEX[rseq[rstart + i - 1] - 'A'];
        oval = BLOSUM62_RESIDUE_TO_INDEX[oseq[ostart + j - 1] - 'A'];

        sdiag = table[(i-1) * c + (j-1)] + BLOSUM62_MATRIX[rval][oval];
        sup = table[(i-1) * c + j] + BLOSUM62_MATRIX[gapi][oval];
        sleft = table[i * c + (j-1)] + BLOSUM62_MATRIX[rval][gapi];

        if (sdiag > sup && sdiag > sleft) {
            i--;
            j--;
            alignment.ref[alignment.length] = rseq[rstart + i];
            alignment.org[alignment.length] = oseq[ostart + j];
            alignment.length++;
        } else if (sup > sleft) {
            i--;
            alignment.ref[alignment.length] = rseq[rstart + i];
            alignment.org[alignment.length] = '-';
            alignment.length++;
        } else {
            j--;
            alignment.ref[alignment.length] = '-';
            alignment.org[alignment.length] = oseq[ostart + j];
            alignment.length++;
        }
    }
    for ( ; i > 0; i--) {
        alignment.ref[alignment.length] = rseq[rstart + i - 1];
        alignment.org[alignment.length] = '-';
        alignment.length++;
    }
    for ( ; j > 0; j--) {
        alignment.ref[alignment.length] = '-';
        alignment.org[alignment.length] = oseq[ostart + j - 1];
        alignment.length++;
    }

    for (i = 0, j = alignment.length - 1; i < j; i++, j--) {
        tmp = alignment.ref[i];
        alignment.ref[i] = alignment.ref[j];
        alignment.ref[j] = tmp;

        tmp = alignment.org[i];
        alignment.org[i] = alignment.org[j];
        alignment.org[j] = tmp;
    }

    if (allocated)
        free(table);

    alignment.ref[alignment.length] = '\0';
    alignment.org[alignment.length] = '\0';

    return alignment;
}

/*Returns the number of non-gap characters in a string*/
int32_t
cbp_align_length_nogaps(char *residues)
{
    int i, len, rlen;

    len = 0;
    rlen = strlen(residues);
    for (i = 0; i < rlen; i++)
        if (residues[i] != '-')
            len++;

    return len;
}

/* returns extension distance (but does not move pointers)*/
int32_t
attempt_ext(int32_t i1, const int32_t dir1, const char *s1, int32_t len1,
            int32_t start1, int32_t i2, const int32_t dir2, const char *s2,
            int32_t len2, int32_t start2)
{
    const int32_t dir_prod = dir1*dir2;
    int32_t progress = 0;
    int32_t consec_mismatch = 0;
    i1 += dir1;
    i2 += dir2;
    /*Replace this 3 with the flag for max_consec_mismatch*/
    while (consec_mismatch < 3 &&
        i1 >= start1 && i1 < start1+len1 && i2 >= start2 && i2 < start2+len2) {
        if (!bases_match(s1[i1], s2[i2], dir_prod))
            consec_mismatch++;
        else
            consec_mismatch = 0;
        i1 += dir1; i2 += dir2;
        progress++;
    }
    return progress;
}

/*Takes in as input an array of bools representing data on whether or not previous
 *pairs of bases matched, an index into that array, the current number of matches,
 *an array of bools to add to the array of matches, and the number of bools to add.
 *check_and_update adds bools from the temp array to the matches array until it
 *either has added all of the matches or a bad window was found.  If a bad window
 *was found, check_and_update returns false.  Otherwise, check_and_update returns true.
 */
int check_and_update(bool *matches, int *matches_index, int *num_matches, bool *temp, int temp_index){
    int i;
    for(i = 0; i < temp_index; i++){
        int hundred_bases_ago = *matches_index - 100;
        (*matches_index)++;
        matches[(*matches_index)] = temp[i];
        
        if(matches[hundred_bases_ago])
            num_matches--;
        if(*num_matches < 85)
            return i;
    }
    return temp_index;
}
