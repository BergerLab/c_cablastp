#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitpack.h"
#include "coarse.h"
#include "DNAutils.h"
#include "seeds.h"

struct cbp_coarse *
cbp_coarse_init(int32_t seed_size,
                FILE *file_fasta, FILE *file_seeds, FILE *file_links)
{
    struct cbp_coarse *coarse_db;
    int32_t errno;

    coarse_db = malloc(sizeof(*coarse_db));
    assert(coarse_db);

    coarse_db->seqs = ds_vector_create_capacity(10000000);
    coarse_db->seeds = cbp_seeds_init(seed_size);
    coarse_db->file_fasta = file_fasta;
    coarse_db->file_seeds = file_seeds;
    coarse_db->file_links = file_links;

    if (0 != (errno = pthread_rwlock_init(&coarse_db->lock_seq, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return coarse_db;
}

void
cbp_coarse_free(struct cbp_coarse *coarse_db)
{
    int32_t errno;
    int32_t i;

    fclose(coarse_db->file_fasta);
    fclose(coarse_db->file_seeds);
    fclose(coarse_db->file_links);

    if (0 != (errno = pthread_rwlock_destroy(&coarse_db->lock_seq))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (i = 0; i < coarse_db->seqs->size; i++)
        cbp_coarse_seq_free(
            (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i));

    ds_vector_free_no_data(coarse_db->seqs);
    cbp_seeds_free(coarse_db->seeds);
    free(coarse_db);
}

struct cbp_coarse_seq *
cbp_coarse_add(struct cbp_coarse *coarse_db,
               char *residues, int32_t start, int32_t end)
{
    struct cbp_coarse_seq *seq;
    int32_t id;

    pthread_rwlock_wrlock(&coarse_db->lock_seq);
    id = coarse_db->seqs->size;
    seq = cbp_coarse_seq_init(id, residues, start, end);
    ds_vector_append(coarse_db->seqs, (void*) seq);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    cbp_seeds_add(coarse_db->seeds, seq);

    return seq;
}

struct cbp_coarse_seq *
cbp_coarse_get(struct cbp_coarse *coarse_db, int32_t i)
{
    struct cbp_coarse_seq *seq;

    pthread_rwlock_rdlock(&coarse_db->lock_seq);
    seq = (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    return seq;
}

/*Outputs the sequences in the coarse database to a FASTA file in plain text
  and outputs the links to the compressed database in a binary format*/
void
cbp_coarse_save_binary(struct cbp_coarse *coarse_db)
{
    struct cbp_coarse_seq *seq;
    struct cbp_link_to_compressed *link;
    int64_t i;
    
    int16_t mask = (((int16_t)1)<<8)-1;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        seq = (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
        fprintf(coarse_db->file_fasta, "> %ld\n%s\n", i, seq->seq->residues);

        fprintf(coarse_db->file_links, "> %ld\n", i);
        for (link = seq->links; link != NULL; link = link->next){
            int j;
            char *id_bytes = read_int_to_bytes(link->org_seq_id, 8);
            for (j = 0; j < 8; j++)
                putc(id_bytes[j], coarse_db->file_links);
            /*Convert the start and end indices for the link to two
              characters.*/
            int16_t start = (int16_t)link->coarse_start;
            int16_t end   = (int16_t)link->coarse_end;
            char start_left  = (start >> 8) & mask;
            char start_right = start & mask;
            char end_left    = (end >> 8) & mask;
            char end_right   = end & mask;
            /*Prints the binary representations of the indices and the
              direction of the link to the links file*/
            putc(start_left, coarse_db->file_links);
            putc(start_right, coarse_db->file_links);
            putc(end_left, coarse_db->file_links);
            putc(end_right, coarse_db->file_links);
            putc((link->dir?'0':'1'), coarse_db->file_links);
            /*0 is used as a delimiter to signify that there are more links
              for this sequence*/
            if(link->next != NULL)
                putc(0, coarse_db->file_links);
        }
        /*'#' is used as a delimiter to signify the last link of the sequence*/
        if(i+1 < coarse_db->seqs->size)
            putc('#', coarse_db->file_links);
    }
}

void
cbp_coarse_save_plain(struct cbp_coarse *coarse_db)
{
    struct cbp_coarse_seq *seq;
    struct cbp_link_to_compressed *link;
    int32_t i;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        seq = (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
        fprintf(coarse_db->file_fasta, "> %d\n%s\n", i, seq->seq->residues);

        fprintf(coarse_db->file_links, "> %d\n", i);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(coarse_db->file_links,
                "original sequence id: %d, reference range: (%d, %d), direction: %c\n",
                link->org_seq_id, link->coarse_start, link->coarse_end,
                (link->dir?'0':'1'));
    }
}

void
cbp_coarse_save_seeds_plain(struct cbp_coarse *coarse_db)
{
    struct cbp_coarse_seq *seq;
    struct cbp_link_to_compressed *link;
    int32_t i, j;
    uint32_t i2;
    uint32_t mask = (uint32_t)3;
    char *kmer;

    for (i = 0; i < coarse_db->seeds->locs_length; i++) {
        struct cbp_seed_loc *loc;
        i2 = (uint32_t)0;
        for (j = 0; j < coarse_db->seeds->seed_size; j++) {
            i2 <<= 2;
            i2 |= ((i >> (2*j)) & mask);
        }
        kmer = unhash_kmer(coarse_db->seeds, i2);
        fprintf(coarse_db->file_seeds, "%s\n", kmer);
        loc = cbp_seeds_lookup(coarse_db->seeds, kmer);
        while (loc) {
            if(loc->coarse_seq_id<500)
                fprintf(coarse_db->file_seeds,"(%d, %d) > ", loc->coarse_seq_id, loc->residue_index);
            loc = loc->next;
        }
        fprintf(coarse_db->file_seeds, "\n");
        free(kmer);
    }
}

struct cbp_coarse_seq *
cbp_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end)
{
    struct cbp_coarse_seq *seq;
    int32_t errno;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->seq = cbp_seq_init_range(id, "", residues, start, end);
    seq->links = NULL;

    if (0 != (errno = pthread_rwlock_init(&seq->lock_links, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return seq;
}

void
cbp_coarse_seq_free(struct cbp_coarse_seq *seq)
{
    struct cbp_link_to_compressed *link1, *link2;
    int32_t errno;

    if (0 != (errno = pthread_rwlock_destroy(&seq->lock_links))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cbp_link_to_compressed_free(link1);
        link1 = link2;
    }

    cbp_seq_free(seq->seq);
    free(seq);
}

void
cbp_coarse_seq_addlink(struct cbp_coarse_seq *seq,
                       struct cbp_link_to_compressed *newlink)
{
/*printf("addlink %d %d %d %d\n", (newlink -> dir == 1 ? 1:-1), newlink->org_seq_id, newlink->coarse_start, newlink->coarse_end);*/
    struct cbp_link_to_compressed *link;

    assert(newlink->next == NULL);

    pthread_rwlock_wrlock(&seq->lock_links);
    if (seq->links == NULL) {
        seq->links = newlink;
        pthread_rwlock_unlock(&seq->lock_links);
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
    pthread_rwlock_unlock(&seq->lock_links);
}

struct cbp_link_to_compressed *
cbp_link_to_compressed_init(int32_t org_seq_id,
                            int16_t coarse_start, int16_t coarse_end, bool dir)
{
    struct cbp_link_to_compressed *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->org_seq_id = org_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->dir = dir;
    link->next = NULL;

    return link;
}

void
cbp_link_to_compressed_free(struct cbp_link_to_compressed *link)
{
    free(link);
}
