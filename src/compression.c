#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"
 
#include "compression.h"
#include "flags.h"
#include "DNAutils.h"
#include "edit_scripts.h"

struct worker_args {
    struct cbp_database *db;
    struct DSQueue *jobs;
};

struct extend_match {
    int32_t rlen;
    int32_t olen;
};

struct extend_match_with_res {
    char *rseq;
    char *oseq;
    int32_t rlen;
    int32_t olen;
};

struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind,
             int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
             int32_t current, int32_t dir2);

struct extend_match_with_res
extend_match_with_res(struct cbp_align_nw_memory *mem,
                char *rseq, int32_t rstart, int32_t rend, int32_t resind,
                int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
                int32_t current, int32_t dir2);

static int32_t
add_without_match(struct cbp_coarse *coarse_db,
                  struct cbp_seq *org_seq, int32_t ostart, int32_t oend);

static void *
cbp_compress_worker(void *data);

static int32_t
min(int32_t a, int32_t b);

struct cbp_compress_workers *
cbp_compress_start_workers(struct cbp_database *db, int32_t num_workers)
{
    struct DSQueue *jobs;
    struct cbp_compress_workers *workers;
    struct worker_args *wargs;
    int32_t i, errno;

    jobs = ds_queue_create(20);

    wargs = malloc(sizeof(*wargs));
    assert(wargs);
    wargs->db = db;
    wargs->jobs = jobs;

    workers = malloc(sizeof(*workers));
    assert(workers);
    workers->threads = malloc(num_workers * sizeof(*workers->threads));
    assert(workers->threads);
    workers->num_workers = num_workers;
    workers->jobs = jobs;
    workers->args = (void*) wargs;

    for (i = 0; i < num_workers; i++) {
        errno = pthread_create(&workers->threads[i], NULL,
            cbp_compress_worker, (void*) wargs);
        if (errno != 0) {
            fprintf(stderr,
                "cbp_compress_start_workers: Could not start thread. Errno: %d",
                errno);
            exit(1);
        }
    }

    return workers;
}

void
cbp_compress_join_workers(struct cbp_compress_workers *workers)
{
    int i, errno;

    ds_queue_close(workers->jobs);

    for (i = 0; i < workers->num_workers; i++) {
        errno = pthread_join(workers->threads[i], NULL);
        if (errno != 0) {
            fprintf(stderr,
                "cbp_compress_join_workers: Could not join thread. Errno: %d",
                errno);
            exit(1);
        }
    }
}

void
cbp_compress_free_workers(struct cbp_compress_workers *workers)
{
    ds_queue_free(workers->jobs);
    free(workers->args);
    free(workers->threads);
    free(workers);
}

void
cbp_compress_send_job(struct cbp_compress_workers *workers,
                      struct cbp_seq *org_seq)
{
    ds_queue_put(workers->jobs, (void*) org_seq);
}

static void *
cbp_compress_worker(void *data)
{
    struct worker_args *args;
    struct cbp_align_nw_memory *mem;
    struct cbp_seq *s;
    struct cbp_compressed_seq *cseq;

    args = (struct worker_args *) data;
    mem = cbp_align_nw_memory_init();
    while (NULL != (s = (struct cbp_seq *) ds_queue_get(args->jobs))) {
        cseq = cbp_compress(args->db->coarse_db, s, mem);
        cbp_compressed_write_binary(args->db->com_db, cseq);

        args->db->coarse_db->dbsize += s->length;

        cbp_seq_free(s);
        cbp_compressed_seq_free(cseq);
    }

    cbp_align_nw_memory_free(mem);

    return NULL;
}


struct cbp_compressed_seq *
cbp_compress(struct cbp_coarse *coarse_db, struct cbp_seq *org_seq,
             struct cbp_align_nw_memory *mem)
{
    fprintf(stderr, "Starting compression      %d\n", org_seq->id);
    struct extend_match mlens_fwd, mlens_rev;
    struct cbp_coarse_seq *coarse_seq;
    struct cbp_compressed_seq *cseq;
    struct cbp_seed_loc *seeds, *seeds_r, *seedLoc;
    struct cbp_alignment alignment;
    char *kmer, *revcomp;
    int32_t seed_size, ext_seed, resind, mext, new_coarse_seq_id, min_progress;
    int32_t last_match, current;
    int32_t i;
    int chunks;

    int32_t max_chunk_size, max_section_size;
    int32_t overlap;

    int32_t start_of_section;
    int32_t end_of_chunk;
    int32_t end_of_section;

    bool changed, found_match;
    bool *matches, *matches_temp;

    cseq = cbp_compressed_seq_init(org_seq->id, org_seq->name);
    seed_size = coarse_db->seeds->seed_size;
    mext = compress_flags.match_extend;
    ext_seed = compress_flags.ext_seed_size;
    min_progress = compress_flags.min_progress;

    max_chunk_size = compress_flags.max_chunk_size;
    max_section_size = 2 * max_chunk_size;
    overlap = compress_flags.overlap;

    last_match = 0;
    current = 0;
    start_of_section = 0;
    end_of_chunk = start_of_section + max_chunk_size;
    end_of_section = start_of_section + max_section_size;
    chunks = 0;

    /*Initialize the matches and matches_temp arrays*/
    matches = malloc(max_section_size*sizeof(*matches));
    matches_temp = malloc(max_section_size*sizeof(*matches_temp));
    for (i = 0; i < max_section_size; i++) {
        matches[i] = true;
        matches_temp[i] = true;
    }

    for (current = 0; current <= org_seq->length-seed_size - ext_seed;
                                                             current++) {
        found_match = false;

        /*If we are at the beginning of the first chunk of the first sequence,
         *add the first chunk without a match and skip ahead to the start of
         *the second chunk.
         */
        if (current == 0 && coarse_db->seqs->size == 0) {
            new_coarse_seq_id = add_without_match(coarse_db, org_seq, 0,
                                                  max_chunk_size);
            cbp_compressed_seq_addlink(cseq, cbp_link_to_coarse_init_nodiff(
                                                 new_coarse_seq_id, 0,
                                                 end_of_chunk - 1, 0,
                                                 end_of_chunk - 1, true));

            if (end_of_chunk < org_seq->length - seed_size - ext_seed) {
                start_of_section += max_chunk_size - overlap;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq->length - ext_seed);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq->length - ext_seed);
                current = start_of_section-1;
            }
            chunks++;
            continue;
        }

        /*Get the k-mer and allocate a copy of its reverse complement*/
        kmer = org_seq->residues + current;
        revcomp = kmer_revcomp(kmer);

/*char *base = kmer;
for(;base < kmer + 10; base++)
    printf("%c", *base);
printf("\n");*/

        /*The locations of all seeds in the database that start with the
          current k-mer.*/
        seeds = cbp_seeds_lookup(coarse_db->seeds, kmer);

        /*The locations of all seeds in the database that start with the
          current k-mer's reverse complement.*/
        seeds_r = cbp_seeds_lookup(coarse_db->seeds, revcomp);

        free(revcomp);

        for (seedLoc = seeds; seedLoc != NULL; seedLoc = seedLoc->next) {
            int coarse_align_len, original_align_len,
                start_coarse_align, end_coarse_align,
                start_original_align, end_original_align;
            char *cor_match, *org_match, *new_cor, *new_org;
            int i1, i2;

            if (found_match)
                break;

            resind = seedLoc->residue_index;
            coarse_seq = cbp_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            if (resind + seed_size + ext_seed > coarse_seq->seq->length)
                continue;

            if ((attempt_ext(current, -1, org_seq->residues,
                            end_of_section - start_of_section,
                            start_of_section, resind, -1,
                            coarse_seq->seq->residues,
                            coarse_seq->seq->length, 0) +
                  attempt_ext(current+seed_size-1, 1, org_seq->residues,
                              end_of_section - start_of_section,
                              start_of_section, resind+seed_size-1, 1,
                              coarse_seq->seq->residues,
                              coarse_seq->seq->length, 0)) >
                  compress_flags.attempt_ext_len) {
                mlens_rev = extend_match(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq->residues, start_of_section,
                                         end_of_section, current, -1);
                            extend_match_with_res(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq->residues, start_of_section,
                                         end_of_section, current, -1);


                mlens_fwd = extend_match(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind + seed_size - 1, 1,
                                         org_seq->residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);

                            extend_match_with_res(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind + seed_size - 1, 1,
                                         org_seq->residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);


/*printf("->  extend_match: %d %d %d\n", mlens_rev.olen, seed_size, mlens_fwd.olen);*/

                /*If the match was too short, try the next seed*/                
                if (mlens_rev.olen + seed_size + mlens_fwd.olen <
                      compress_flags.min_match_len)
                    continue;

                found_match = true;
                changed = false;
                coarse_align_len = mlens_rev.rlen + seed_size + mlens_fwd.rlen;
                original_align_len=mlens_rev.olen + seed_size + mlens_fwd.olen;

                start_coarse_align = resind - mlens_rev.rlen;
                end_coarse_align = resind + seed_size + mlens_fwd.rlen;
                start_original_align = current - mlens_rev.olen;
                end_original_align = current + seed_size + mlens_fwd.olen;

                cor_match = malloc(coarse_align_len * sizeof(*cor_match));
                org_match = malloc(original_align_len * sizeof(*org_match));

                /*Copy the matching parts of the coarse and original sequences*/
                for (i1 = start_coarse_align; i1 < end_coarse_align; i1++)
                    cor_match[i1-start_coarse_align] =
                      coarse_seq->seq->residues[i1];
                for (i2 = start_original_align; i2 < end_original_align; i2++)
                    org_match[i2-start_original_align] = org_seq->residues[i2];

                /*Get an alignment of the matching sequences*/
                alignment = cbp_align_nw(mem,
                                         cor_match, coarse_align_len, 0, 1,
                                         org_match, original_align_len, 0, 1,
                                         matches, NULL);
                free(cor_match);
                free(org_match);

                /*If the length of either the coarse or original match
                  was changed from backtracking in the alignment, update
                  the lengths.*/
                new_cor = no_dashes(alignment.ref);
                new_org = no_dashes(alignment.org);
                for (i1 = 0; new_cor[i1] != '\0'; i1++);
                for (i2 = 0; new_org[i2] != '\0'; i2++);
                if (i1 < coarse_align_len)
                    mlens_fwd.rlen -= (coarse_align_len - i1);
                if (i2 < original_align_len)
                    mlens_fwd.olen -= (original_align_len - i2);
                free(new_cor);
                free(new_org);
                
                /*Make a new chunk for the parts of the chunk before the
                  match.*/
                if (current - mlens_rev.olen - start_of_section > 0) {
                    new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                    start_of_section,
                                                    current - mlens_rev.olen +
                                                      compress_flags.overlap);
                    cbp_compressed_seq_addlink(cseq,
                        cbp_link_to_coarse_init_nodiff(
                            new_coarse_seq_id,
                            start_of_section,
                            current - mlens_rev.olen
                                    + compress_flags.overlap - 1,
                            0, current - mlens_rev.olen
                                       + compress_flags.overlap
                                       - start_of_section - 1,
                            true));
                    chunks++;
                }

                /*Add a link to the coarse sequence in the compressed
                  sequence.*/
                cbp_compressed_seq_addlink(cseq,
                    cbp_link_to_coarse_init(coarse_seq->id,
                                            current - mlens_rev.olen,
                                            current+seed_size+mlens_fwd.olen-1,
                                            resind - mlens_rev.rlen,
                                            resind+seed_size+mlens_fwd.rlen-1,
                                            alignment, true));

                /*Add a link to the compressed sequence in the coarse
                  sequence.*/
                cbp_coarse_seq_addlink(coarse_seq,
                                       cbp_link_to_compressed_init(
                                       org_seq->id,
                                       resind - mlens_rev.rlen,
                                       resind + seed_size + mlens_fwd.rlen-1,
                                       current - mlens_rev.olen,
                                       current + seed_size + mlens_fwd.olen-1,
                                       true));
                /*Update the current position in the sequence*/
                if (current + mlens_fwd.olen <
                      org_seq->length - seed_size - ext_seed - 1)
                    start_of_section = current + mlens_fwd.olen -
                                       compress_flags.overlap + seed_size;
                else
                    start_of_section = current + mlens_fwd.olen + seed_size;

                current = start_of_section-1;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq->length-ext_seed);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq->length-ext_seed);

                chunks++;
            }
        }
        for (seedLoc = seeds_r; seedLoc != NULL; seedLoc = seedLoc->next) {
            int coarse_align_len, original_align_len,
                start_coarse_align, end_coarse_align,
                start_original_align, end_original_align;
            char *cor_match, *org_match, *new_cor, *new_org;
            int i1, i2;

            /*If we found a match in the seed locations for the k-mer, then
             *there is no need to check the locations for the reverse
             *complement.
             */
            if (found_match)
                break;

            resind = seedLoc->residue_index;
            coarse_seq = cbp_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            if (resind + seed_size + ext_seed > coarse_seq->seq->length)
                continue;

            if ((attempt_ext(current, -1, org_seq->residues,
                            end_of_section - start_of_section,
                            start_of_section, resind + seed_size - 1, 1,
                            coarse_seq->seq->residues,
                            coarse_seq->seq->length, 0) +
                  attempt_ext(current+seed_size-1, 1, org_seq->residues, 
                              end_of_section - start_of_section,
                              start_of_section, resind, -1,
                              coarse_seq->seq->residues,
                              coarse_seq->seq->length, 0)) >
                  compress_flags.attempt_ext_len) {
                mlens_rev = extend_match(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq->residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);

                            extend_match_with_res(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq->residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);


                mlens_fwd = extend_match(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind+seed_size-1, 1,
                                         org_seq->residues, start_of_section,
                                         end_of_section, current, -1);

                            extend_match_with_res(mem, coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind+seed_size-1, 1,
                                         org_seq->residues, start_of_section,
                                         end_of_section, current, -1);

/*printf("<-  extend_match: %d %d %d\n", mlens_fwd.olen, seed_size, mlens_rev.olen);*/

                /*If the match was too short, try the next seed*/                
                if (mlens_rev.olen+seed_size+mlens_fwd.olen <
                      compress_flags.min_match_len)
                    continue;

                found_match = true;
                changed = false;

                coarse_align_len = mlens_rev.rlen + seed_size + mlens_fwd.rlen;
                original_align_len=mlens_rev.olen + seed_size + mlens_fwd.olen;

                start_coarse_align = resind - mlens_rev.rlen;
                end_coarse_align = resind + seed_size + mlens_fwd.rlen;
                start_original_align = current - mlens_fwd.olen;
                end_original_align = current + seed_size + mlens_rev.olen;
                
                cor_match = malloc(coarse_align_len * sizeof(*cor_match));
                org_match = malloc(original_align_len * sizeof(*org_match));

                /*Copy the matching parts of the coarse and original sequences*/
                for (i1 = start_coarse_align; i1 < end_coarse_align; i1++)
                    cor_match[i1-start_coarse_align] =
                      coarse_seq->seq->residues[i1]; 
                for (i2 = start_original_align; i2 < end_original_align; i2++)
                    org_match[i2-start_original_align] = org_seq->residues[i2];

                /*Get an alignment of the matching sequences*/
                alignment = cbp_align_nw(mem,
                                         cor_match, coarse_align_len, 0, 1,
                                         org_match, original_align_len,
                                         original_align_len-1, -1,
                                         matches, NULL);
                free(cor_match);
                free(org_match);

                /*If the length of either the coarse or original match
                  was changed from backtracking in the alignment, update
                  the lengths.*/
                new_cor = no_dashes(alignment.ref);
                new_org = no_dashes(alignment.org);
                for (i1 = 0; new_cor[i1] != '\0'; i1++);
                for (i2 = 0; new_org[i2] != '\0'; i2++);
                if (i1 < coarse_align_len)
                    mlens_fwd.rlen -= (coarse_align_len - i1);
                if (i2 < original_align_len)
                    mlens_fwd.olen -= (original_align_len - i2);
                free(new_cor);
                free(new_org);

                /*Make a new chunk for the parts of the chunk before the
                  match.*/
                if (current - mlens_fwd.olen - start_of_section > 0) {
                    new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                    start_of_section,
                                                    current - mlens_fwd.olen +
                                                      compress_flags.overlap);
                    cbp_compressed_seq_addlink(cseq,
                        cbp_link_to_coarse_init_nodiff(
                            new_coarse_seq_id,
                            start_of_section,
                            current - mlens_fwd.olen
                                    + compress_flags.overlap - 1,
                            0, current - mlens_fwd.olen
                                       + compress_flags.overlap
                                       - start_of_section - 1,
                        true));
                    chunks++;
                }

                /*Add a link to the coarse sequence in the compressed
                  sequence.*/
                cbp_compressed_seq_addlink(cseq,
                    cbp_link_to_coarse_init(coarse_seq->id,
                                            current - mlens_fwd.olen,
                                            current + seed_size +
                                                      mlens_rev.olen - 1,
                                            resind - mlens_rev.rlen,
                                            resind + seed_size +
                                                     mlens_fwd.rlen - 1,
                                            alignment, false));

                /*Add a link to the compressed sequence in the coarse
                  sequence.*/
                cbp_coarse_seq_addlink(coarse_seq,
                                       cbp_link_to_compressed_init(
                                       org_seq->id,
                                       resind - mlens_rev.rlen,
                                       resind + seed_size + mlens_fwd.rlen-1,
                                       current - mlens_fwd.olen,
                                       current + seed_size + mlens_rev.olen - 1,
                                       false));

                /*Update the current position in the sequence*/
                if (current + mlens_rev.olen <
                      org_seq->length - seed_size - ext_seed - 1)
                    start_of_section = current + mlens_rev.olen -
                                       compress_flags.overlap + seed_size;
                else
                    start_of_section = current + mlens_rev.olen + seed_size;
                current = start_of_section - 1;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq->length-ext_seed);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq->length-ext_seed);

                chunks++;
            }
        }
        cbp_seed_loc_free(seeds);
        cbp_seed_loc_free(seeds_r);

        /*If we have traversed an entire chunk of bases without finding a match,
         *then add the whole chunk as a sequence in the database and update
         *start_of_section, end_of_chunk, and end_of_section
         */
        if (current >= end_of_chunk - seed_size && !found_match) {
            new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                  start_of_section,
                                                  end_of_chunk);

            cbp_compressed_seq_addlink(cseq, cbp_link_to_coarse_init_nodiff(
                                                 new_coarse_seq_id,
                                                 start_of_section,
                                                 end_of_chunk - 1, 0,
                                                 end_of_chunk -
                                                   start_of_section - 1,
                                                 true));

            if (end_of_chunk < org_seq->length - seed_size - ext_seed - 1) {
                start_of_section = end_of_chunk - overlap;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq->length - ext_seed);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq->length - ext_seed);
                current = start_of_section - 1;
            }
            chunks++;
        }
if(chunks >= 4)break;
    }
    fprintf(stderr, "Compress finished       %d\n", org_seq->id);
    free(matches);
    free(matches_temp);
    return cseq;
}

struct extend_match_with_res
extend_match_with_res(struct cbp_align_nw_memory *mem,
                char *rseq, int32_t rstart, int32_t rend, int32_t resind,
                int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
                int32_t current, int32_t dir2)
{
fprintf(stderr, "%d %d-%d %c, %d %d-%d %c\n", resind, rstart, rend, (dir1>0?'+':'-'), current, ostart, oend, (dir2>0?'+':'-'));
    struct cbp_alignment alignment;
    struct extend_match_with_res mseqs;
    int32_t rlen, olen;
    struct ungapped_alignment ungapped;
    int32_t m;
    bool *matches;
    bool *matches_past_clump;
    int matches_count;
    int matches_index;
    int max_section_size;
    int i, j;
    int rseq_len = 0, oseq_len = 0;
    bool found_bad_window;
    struct DSVector *rseq_segments = ds_vector_create();
    struct DSVector *oseq_segments = ds_vector_create();

    max_section_size = 2 * compress_flags.max_chunk_size;

    /*Initialize the matches and matches_past_clump arrays.*/
    matches = malloc(2*compress_flags.max_chunk_size*sizeof(*matches));
    matches_past_clump = malloc(2 * compress_flags.max_chunk_size
                                  * sizeof(*matches_past_clump));
    matches_index = compress_flags.gapped_window_size;
    for (i = 0; i < max_section_size; i++) {
        matches[i] = true;
        matches_past_clump[i] = true;
    }

    resind += dir1;
    current += dir2;

    rlen = rend - rstart;
    olen = oend - ostart;

    mseqs.rlen = 0;
    mseqs.olen = 0;
    mseqs.rseq = "";
    mseqs.oseq = "";
fprintf(stderr, "!!!");
    while (true) {
        int dp_len1, dp_len2, i, r_align_len, o_align_len;
        char *r_segment, *o_segment;
        if (mseqs.rlen == rlen || mseqs.olen == olen)
            break;

        /*Get the maximum length for ungapped alignment and extend the match
          by that distance.*/
        ungapped = cbp_align_ungapped(rseq, rstart, rend, dir1, resind,
                               oseq, ostart, oend, dir2, current,
                               matches, matches_past_clump, &matches_index);

        m = ungapped.length;
        found_bad_window = ungapped.found_bad_window;

        mseqs.rlen += m;
        mseqs.olen += m;

        r_segment = malloc((m+1)*sizeof(*r_segment));
        o_segment = malloc((m+1)*sizeof(*o_segment));
        for (i = 0; i < m; i++) {
            r_segment[i] = dir1 ? rseq[resind+i] :
                                  base_complement(rseq[resind-i]);
            o_segment[i] = dir2 ? oseq[current+i] :
                                  base_complement(oseq[current-i]);
        }
        r_segment[m] = '\0';
        o_segment[m] = '\0';

        ds_vector_append(rseq_segments, (void *)r_segment);
        ds_vector_append(oseq_segments, (void *)o_segment);

        resind += m * dir1;
        current += m * dir2;

        /*End the extension if we found a bad window in ungapped alignment.*/
        if (found_bad_window)
            break;

        /*Carry out Needleman-Wunsch alignment and end the extension if we
          found a bad window or couldn't find a 4-mer match in the alignment.*/
        dp_len1 = max_dp_len(resind-rstart, dir1, rend-rstart);
        dp_len2 = max_dp_len(current-ostart, dir2, oend-ostart);

        alignment = cbp_align_nw(mem, rseq, dp_len1, resind, dir1,
                                      oseq, dp_len2, current, dir2,
                                 matches, &matches_index);

        if (alignment.length == -1)
            break;

        matches_count = 0;

        /*Check for a bad window manually and end the extension if a bad
          window is found.*/
        for (i = matches_index - 100; i < matches_index; i++)
            if (matches[i])
                matches_count++;
        if (matches_count < compress_flags.window_ident_thresh)
            break;

        ds_vector_append(rseq_segments, (void *)alignment.ref);
        ds_vector_append(oseq_segments, (void *)alignment.org);

        r_align_len = cbp_align_length_nogaps(alignment.ref);
        o_align_len = cbp_align_length_nogaps(alignment.org);

        /*Update the lengths of the alignments and the indices of the
          sequences.*/
        mseqs.rlen += r_align_len;
        mseqs.olen += o_align_len;
        resind += r_align_len * dir1;
        current += o_align_len * dir2;
    }

    mseqs.rseq = malloc((mseqs.rlen+1)*sizeof(*(mseqs.rseq)));
    mseqs.oseq = malloc((mseqs.olen+1)*sizeof(*(mseqs.oseq)));
    for (i = 0; i < rseq_segments->size; i++) {
        char *r_segment = (char *)ds_vector_get(rseq_segments, i);
        char *o_segment = (char *)ds_vector_get(oseq_segments, i);
        for (j = 0; r_segment[j] != '\0'; j++)
            mseqs.rseq[rseq_len++] = r_segment[j];
        for (j = 0; o_segment[j] != '\0'; j++)
            mseqs.oseq[oseq_len++] = o_segment[j];
    }
    mseqs.rseq[mseqs.rlen] = '\0';
    mseqs.oseq[mseqs.olen] = '\0';

    free(matches);
    free(matches_past_clump);
    return mseqs;
}


struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind,
             int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
             int32_t current, int32_t dir2)
{
    struct cbp_alignment alignment;
    struct extend_match mlens;
    int32_t rlen, olen;
    struct ungapped_alignment ungapped;
    int32_t m;
    bool *matches;
    bool *matches_past_clump;
    int matches_count;
    int matches_index;
    int max_section_size;
    int i;
    bool found_bad_window;

    max_section_size = 2 * compress_flags.max_chunk_size;

    /*Initialize the matches and matches_past_clump arrays.*/
    matches = malloc(2*compress_flags.max_chunk_size*sizeof(*matches));
    matches_past_clump = malloc(2 * compress_flags.max_chunk_size
                                  * sizeof(*matches_past_clump));
    matches_index = compress_flags.gapped_window_size;
    for (i = 0; i < max_section_size; i++) {
        matches[i] = true;
        matches_past_clump[i] = true;
    }

    resind += dir1;
    current += dir2;

    rlen = rend - rstart;
    olen = oend - ostart;

    mlens.rlen = 0;
    mlens.olen = 0;
    while (true) {
        int dp_len1, dp_len2, i, r_align_len, o_align_len;
        if (mlens.rlen == rlen || mlens.olen == olen)
            break;

        /*Get the maximum length for ungapped alignment and extend the match
          by that distance.*/
        ungapped = cbp_align_ungapped(rseq, rstart, rend, dir1, resind,
                               oseq, ostart, oend, dir2, current,
                               matches, matches_past_clump, &matches_index);
/*printf("align_ungapped finished: %d\n", ungapped.length);*/
        m = ungapped.length;
        found_bad_window = ungapped.found_bad_window;
        mlens.rlen += m;
        mlens.olen += m;
        resind += m * dir1;
        current += m * dir2;

        /*End the extension if we found a bad window in ungapped alignment.*/
        if (found_bad_window)
            break;

        /*Carry out Needleman-Wunsch alignment and end the extension if we
          found a bad window or couldn't find a 4-mer match in the alignment.*/
        dp_len1 = max_dp_len(resind-rstart, dir1, rend-rstart);
        dp_len2 = max_dp_len(current-ostart, dir2, oend-ostart);
/*printf("dp_len1: %d, dp_len2: %d\n", dp_len2, dp_len1);*/

        alignment = cbp_align_nw(mem, rseq, dp_len1, resind, dir1,
                                      oseq, dp_len2, current, dir2,
                                 matches, &matches_index);

        if (alignment.length == -1)
            break;

        matches_count = 0;

        /*Check for a bad window manually and end the extension if a bad
          window is found.*/
        for (i = matches_index - 100; i < matches_index; i++)
            if (matches[i])
                matches_count++;
        if (matches_count < compress_flags.window_ident_thresh)
            break;

        r_align_len = cbp_align_length_nogaps(alignment.ref);
        o_align_len = cbp_align_length_nogaps(alignment.org);
/*printf("r_align_len: %d, o_align_len: %d\n", r_align_len, o_align_len);*/
        /*Update the lengths of the alignments and the indices of the
          sequences.*/
        mlens.rlen += r_align_len;
        mlens.olen += o_align_len;
        resind += r_align_len * dir1;
        current += o_align_len * dir2;
        free(alignment.org);
        free(alignment.ref);
    }
    free(matches);
    free(matches_past_clump);
    return mlens;
}

/*Creates a new coarse sequence for a section of the original DNA sequence that
 *does not have a match, such as if the maximum chunk length was traversed in
 *the sequence without finding a match or if there is any DNA in the sequence
 *before the latest match that isn't part of the match.
 *
 *In addition to creating a new coarse sequence, add_without_match also adda a
 *a link to the sequence being compressed to the new coarse sequence.
 */
static int32_t
add_without_match(struct cbp_coarse *coarse_db,
                  struct cbp_seq *org_seq, int32_t ostart, int32_t oend)
{
    struct cbp_link_to_compressed *link = NULL;
    struct cbp_coarse_seq *coarse_seq =
        cbp_coarse_add(coarse_db, org_seq->residues, ostart, oend);
    cbp_coarse_seq_addlink(
        coarse_seq,
        cbp_link_to_compressed_init(org_seq->id, 0, oend - ostart - 1,
                                    ostart, oend - 1, true));
    link = coarse_seq->links;
    return coarse_seq->id;
}

static int32_t
min(int32_t a, int32_t b)
{
    if (a < b)
        return a;
    return b;
}
