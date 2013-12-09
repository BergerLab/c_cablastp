#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compression.h"
#include "flags.h"
#include "DNAutils.h"

struct worker_args {
    struct cbp_database *db;
    struct DSQueue *jobs;
};

struct extend_match {
    int32_t rlen;
    int32_t olen;
};

struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind,  int32_t dir1,
             char *oseq, int32_t ostart, int32_t oend, int32_t current, int32_t dir2);

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

struct cbp_compressed_seq *
cbp_compress(struct cbp_coarse *coarse_db, struct cbp_seq *org_seq,
             struct cbp_align_nw_memory *mem)
{
    struct extend_match mlens_fwd, mlens_rev;
    struct cbp_coarse_seq *coarse_seq;
    struct cbp_compressed_seq *cseq;
    struct cbp_seed_loc *seeds, *seeds_r, *seedLoc;
    struct cbp_alignment alignment;
    char *kmer, *revcomp;
    int32_t seed_size, ext_seed, resind, mext, new_coarse_seq_id, min_progress;
    int32_t last_match, current;
    int32_t id;
    int32_t i;
    int chunks;

    int32_t max_chunk_size, max_section_size;
    int32_t overlap;

    int32_t start_of_section;
    int32_t end_of_chunk;
    int32_t end_of_section;

    bool has_end, changed, found_match;
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
    end_of_chunk = start_of_section + max_chunk_size - 1;
    end_of_section = start_of_section + max_section_size - 1;
    chunks = 0;

    found_match = false;

    matches = malloc(max_section_size*sizeof(*matches));
    matches_temp = malloc(max_section_size*sizeof(*matches_temp));
    for(i = 0; i < max_section_size; i++){
        matches[i] = true;
        matches_temp[i] = true;
    }
    for (current = 0; current < org_seq->length - seed_size - ext_seed; current++) {
if(chunks >= 42)break;
        if(current == 0 && coarse_db->seqs->size == 0){
            add_without_match(coarse_db, org_seq, 0, max_chunk_size);
            start_of_section += max_chunk_size - overlap;
            end_of_chunk = min(start_of_section + max_chunk_size,
                               org_seq->length - seed_size - ext_seed);
            end_of_section = min(start_of_section + max_section_size,
                                 org_seq->length - seed_size - ext_seed);
            current = start_of_section-1;
/********print_seeds(coarse_db->seeds);********/
            chunks++;
            continue;
        }
        kmer = org_seq->residues + current;
	revcomp = kmer_revcomp(kmer);
        seeds = cbp_seeds_lookup(coarse_db->seeds, kmer);
        seeds_r = cbp_seeds_lookup(coarse_db->seeds, revcomp);

        for (seedLoc = seeds; seedLoc != NULL; seedLoc = seedLoc->next) {
char *base = kmer;
            resind = seedLoc->residue_index;
            coarse_seq = cbp_coarse_get(coarse_db, seedLoc->coarse_seq_id);
            if (resind + seed_size + ext_seed > coarse_seq->seq->length)
                continue;
            if (0 != strncmp(
                        coarse_seq->seq->residues + seed_size,
                        org_seq->residues + seed_size,
                        ext_seed))
                continue;
            if(attempt_ext(current, -1, org_seq->residues, end_of_section - start_of_section, start_of_section+1,
                           resind, -1, coarse_seq->seq->residues, coarse_seq->seq->length, 0) +
               attempt_ext(current+seed_size-1, 1, org_seq->residues, end_of_section - start_of_section, start_of_section+1,
                           resind+seed_size-1, 1, coarse_seq->seq->residues, coarse_seq->seq->length, 0) > 50){

printf("-->\n");
                mlens_rev = extend_match(mem, coarse_seq->seq->residues, 0, coarse_seq->seq->length, resind, -1,
                                          org_seq->residues, start_of_section, end_of_section, current, -1);
                mlens_fwd = extend_match(mem, coarse_seq->seq->residues, 0, coarse_seq->seq->length, resind+seed_size-1, 1,
                                          org_seq->residues, start_of_section, end_of_section, current+seed_size-1, 1);
                /*printf(">>>%d %d\n", mlens_rev.olen, mlens_fwd.olen);*/
            }


            /*mlens = extend_match(
                mem,
                coarse_seq->seq->residues, resind, coarse_seq->seq->length,
                org_seq->residues, current, org_seq->length);

            has_end = mlens.olen + mext >= org_seq->length - current;
            if (mlens.olen < compress_flags.min_match_len && !has_end)
                continue;

            alignment = cbp_align_nw(
                mem,
                coarse_seq->seq->residues, resind, resind + mlens.rlen,
                org_seq->residues, current, current + mlens.olen);
            id = cbp_align_identity(
                alignment.ref, 0, alignment.length,
                alignment.org, 0, alignment.length);
            if (id < compress_flags.match_seq_id_threshold)
                continue;

            changed = false;
            if (has_end) {
                mlens.olen = org_seq->length - current;
                changed = true;
            }
            if (current - last_match <= mext) {
                mlens.olen += current - last_match;
                current = last_match;
                changed = true;
            }

            if (changed)
                alignment = cbp_align_nw(
                    mem,
                    coarse_seq->seq->residues, resind, resind + mlens.rlen,
                    org_seq->residues, current, current + mlens.olen);

            if (current - last_match > 0) {
                new_coarse_seq_id = add_without_match(
                    coarse_db, org_seq, last_match, current);
                cbp_compressed_seq_addlink(
                    cseq,
                    cbp_link_to_coarse_init_nodiff(
                        new_coarse_seq_id, 0, current - last_match));
            }

            cbp_compressed_seq_addlink(
                cseq,
                cbp_link_to_coarse_init(
                    coarse_seq->id, resind, resind + mlens.rlen, alignment));
            cbp_coarse_seq_addlink(
                coarse_seq,
                cbp_link_to_compressed_init(
                    org_seq->id, resind, resind + mlens.rlen));

            last_match = current + mlens.olen;
            current = last_match - 1;

            found_match = true;
            break;*/
        }
        for (seedLoc = seeds_r; seedLoc != NULL; seedLoc = seedLoc->next) {
            if(found_match)
              break;
            resind = seedLoc->residue_index;
            coarse_seq = cbp_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            if (resind + seed_size + ext_seed > coarse_seq->seq->length)
                continue;

            /*if (0 != strncmp(
                        coarse_seq->seq->residues + seed_size,
                        org_seq->residues + seed_size,
                        ext_seed))
                continue;*/

            if(attempt_ext(current, -1, org_seq->residues, end_of_section - start_of_section, start_of_section+1,
                            resind+seed_size-1, 1, coarse_seq->seq->residues, coarse_seq->seq->length, 0) +
               attempt_ext(current+seed_size-1, 1, org_seq->residues, end_of_section - start_of_section, start_of_section+1,
                            resind, -1, coarse_seq->seq->residues, coarse_seq->seq->length, 0) > 50){
printf("<--\n");
                mlens_rev = extend_match(mem, coarse_seq->seq->residues, 0, coarse_seq->seq->length, resind, -1,
                                         org_seq->residues, start_of_section, end_of_section, current+seed_size-1, 1);
                mlens_fwd = extend_match(mem, coarse_seq->seq->residues, 0, coarse_seq->seq->length, resind+seed_size-1, 1,
                                         org_seq->residues, start_of_section, end_of_section, current, -1);
                /*printf("<<<%d %d\n", mlens_rev.olen, mlens_fwd.olen);*/
            }

            /*mlens = extend_match(
                mem,
                coarse_seq->seq->residues, resind, coarse_seq->seq->length,
                org_seq->residues, current, org_seq->length);

            has_end = mlens.olen + mext >= org_seq->length - current;
            if (mlens.olen < compress_flags.min_match_len && !has_end)
                continue;

            alignment = cbp_align_nw(
                mem,
                coarse_seq->seq->residues, resind, resind + mlens.rlen,
                org_seq->residues, current, current + mlens.olen);
            id = cbp_align_identity(
                alignment.ref, 0, alignment.length,
                alignment.org, 0, alignment.length);
            if (id < compress_flags.match_seq_id_threshold)
                continue;

            changed = false;
            if (has_end) {
                mlens.olen = org_seq->length - current;
                changed = true;
            }
            if (current - last_match <= mext) {
                mlens.olen += current - last_match;
                current = last_match;
                changed = true;
            }

            if (changed)
                alignment = cbp_align_nw(
                    mem,
                    coarse_seq->seq->residues, resind, resind + mlens.rlen,
                    org_seq->residues, current, current + mlens.olen);

            if (current - last_match > 0) {
                new_coarse_seq_id = add_without_match(
                    coarse_db, org_seq, last_match, current);
                cbp_compressed_seq_addlink(
                    cseq,
                    cbp_link_to_coarse_init_nodiff(
                        new_coarse_seq_id, 0, current - last_match));
            }

            cbp_compressed_seq_addlink(
                cseq,
                cbp_link_to_coarse_init(
                    coarse_seq->id, resind, resind + mlens.rlen, alignment));
            cbp_coarse_seq_addlink(
                coarse_seq,
                cbp_link_to_compressed_init(
                    org_seq->id, resind, resind + mlens.rlen));

            last_match = current + mlens.olen;
            current = last_match - 1;

            break;*/
        }

        cbp_seed_loc_free(seeds);
        cbp_seed_loc_free(seeds_r);
if(current >= end_of_chunk){
    add_without_match(coarse_db, org_seq, start_of_section, end_of_chunk);
    start_of_section = end_of_chunk - overlap;
    end_of_chunk = min(start_of_section + max_chunk_size,
                       org_seq->length - seed_size - ext_seed);
    end_of_section = min(start_of_section + max_section_size,
                         org_seq->length - seed_size - ext_seed);
    current = start_of_section - 1;
    chunks++;
}
    }

    if (org_seq->length - last_match > 0) {
        new_coarse_seq_id = add_without_match(
            coarse_db, org_seq, last_match, org_seq->length);
        cbp_compressed_seq_addlink(
            cseq,
            cbp_link_to_coarse_init_nodiff(
                new_coarse_seq_id, 0, org_seq->length - last_match));
    }
fprintf(stderr, "Compress finished\n");
    return cseq;
}

struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind, int32_t dir1,
             char *oseq, int32_t ostart, int32_t oend, int32_t current, int32_t dir2)
{
    struct cbp_alignment alignment;
    struct extend_match mlens;
    int32_t id;
    int32_t gwsize;
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

    matches = malloc(2*compress_flags.max_chunk_size*sizeof(bool));
    matches_past_clump = malloc(2*compress_flags.max_chunk_size*sizeof(bool));
    matches_index = compress_flags.gapped_window_size;
    max_section_size = 2 * compress_flags.max_chunk_size;

    resind += dir1;
    current += dir2;

    for(i = 0; i < max_section_size; i++){
        matches[i] = true;
        matches_past_clump[i] = true;
    }

    gwsize = compress_flags.gapped_window_size;
    rlen = rend - rstart;
    olen = oend - ostart;

    mlens.rlen = 0;
    mlens.olen = 0;
    while (true) {
        int dp_len1, dp_len2, i, r_align_len, o_align_len;
        if (mlens.rlen == rlen || mlens.olen == olen)
            break;
        ungapped = cbp_align_ungapped(rseq, rstart, rend, dir1, resind,
                               oseq, ostart, oend, dir2, current,
                               matches, matches_past_clump, &matches_index);
        m = ungapped.length;
        found_bad_window = ungapped.found_bad_window;
        mlens.rlen += m;
        mlens.olen += m;
        resind += m * dir1;
        current += m * dir2;
        if(found_bad_window)
            break;
        dp_len1 = max_dp_len(resind-rstart, dir1, rend-rstart);
        dp_len2 = max_dp_len(current-ostart, dir2, oend-ostart);
printf("%d@@@%d\n", dp_len2, dp_len1);
        alignment = cbp_align_nw(mem, rseq, dp_len1, resind, dir1,
                                      oseq, dp_len2, current, dir2,
                                 matches, &matches_index);
        if(alignment.length == -1)
            break;
        printf("%s\n%s\n", alignment.org, alignment.ref);
        matches_count = 0;
        for(i = matches_index - 100; i < matches_index; i++)
            if(matches[i])
                matches_count++;
        if(matches_count < compress_flags.window_ident_thresh)
            break;

        r_align_len = cbp_align_length_nogaps(alignment.ref);
        o_align_len = cbp_align_length_nogaps(alignment.org);
        
        mlens.rlen += r_align_len;
        mlens.olen += o_align_len;
        resind += r_align_len * dir1;
        current += o_align_len * dir2;
    }
    return mlens;
}

static int32_t
add_without_match(struct cbp_coarse *coarse_db,
                  struct cbp_seq *org_seq, int32_t ostart, int32_t oend)
{
    struct cbp_coarse_seq *coarse_seq;

    coarse_seq = cbp_coarse_add(coarse_db, org_seq->residues, ostart, oend);
    cbp_coarse_seq_addlink(
        coarse_seq,
        cbp_link_to_compressed_init(org_seq->id, 0, oend - ostart));
    return coarse_seq->id;
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
        cbp_compressed_write(args->db->com_db, cseq);
        cbp_seq_free(s);
        cbp_compressed_seq_free(cseq);
    }
    cbp_align_nw_memory_free(mem);

    return NULL;
}

static int32_t
min(int32_t a, int32_t b)
{
    if (a < b)
        return a;
    return b;
}
