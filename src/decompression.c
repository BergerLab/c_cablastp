#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ds.h"
#include "stdbool.h"
#include "stdlib.h"

#include "coarse.h"
#include "compressed.h"
#include "decompression.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"

/*Takes in an entry in a compressed database and a coarse database and returns
  the decompressed sequence that the database entry came from as a pointer to
  a struct cbp_seq.*/
struct cbp_seq *cbp_decompress_seq(struct cbp_compressed_seq *cseq,
                                   struct cbp_coarse *coarsedb){
    uint64_t last_end = 0;
    int overlap;
    struct cbp_link_to_coarse *link = NULL;
    struct cbp_seq *seq = NULL;
    int decompressed_length = 0;
    int i = 0;
    int j = 0;
    int copied = 0;
    struct DSVector *decompressed_chunks = ds_vector_create();
    for (link = cseq->links; link != NULL; link = link->next) {
        int length;
        char *dec_chunk;
        struct fasta_seq *chunk = cbp_coarse_read_fasta_seq(coarsedb,
                                                           link->coarse_seq_id);

        int coarse_len = link->coarse_end - link->coarse_start;
        char *coarse_sub = malloc((coarse_len+1)*sizeof(*coarse_sub));
        memcpy(coarse_sub, chunk->seq + link->coarse_start, coarse_len);
        coarse_sub[coarse_len] = '\0';

        for (length = 0; chunk->seq[length] != '\0'; length++);

        /*overlap represents the length of the overlap of the parts of the
          decompressed sequence that has been printed and the parts of the
          decompressed sequence currently being decompressed.*/
        overlap = last_end - link->original_start;
        dec_chunk = read_edit_script(link->diff, coarse_sub, length);

        /*Print all characters of the decompressed chunk past the index
          "overlap" unless overlap is greater than the length of the
          decompressed chunk.*/
        dec_chunk += overlap;
        if ((unsigned int)overlap < link->original_end - link->original_start) {
            int chunk_length = 0;
            char *section = NULL;

            for (chunk_length=0; dec_chunk[chunk_length]!='\0'; chunk_length++);

            section = malloc((chunk_length + 1)*sizeof(*section));
            for (i = 0; i <= chunk_length; i++)
                section[i] = dec_chunk[i];
            ds_vector_append(decompressed_chunks, (void *)section);
            decompressed_length += chunk_length;
            if (link->original_end > last_end)
                last_end = link->original_end;
        }

        dec_chunk -= overlap;
        free(dec_chunk);
        free(coarse_sub);
        fasta_free_seq(chunk);
    }
    char *residues = malloc((decompressed_length+1)*sizeof(*residues));
    for (i = 0; i < decompressed_chunks->size; i++) {
        char *current_chunk = (char *)ds_vector_get(decompressed_chunks, i);
        for (j = 0; current_chunk[j] != '\0'; j++)
            residues[copied++] = current_chunk[j];
    }
    residues[copied] = '\0';
    seq = cbp_seq_init(cseq->id, cseq->name, residues);
    ds_vector_free(decompressed_chunks);
    free(residues);
    return seq;
}

int get_min(int a, int b){return a<b?a:b;}
int get_max(int a, int b){return a>b?a:b;}

struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length){
fprintf(stderr, "==========================================cbp_coarse_expand %d   %d-%d===================================================\n", id, start, end);
    FILE *links = coarsedb->file_links;
    FILE *coarse_links_index = coarsedb->file_links_index;
    FILE *fasta = coarsedb->file_fasta;
    FILE *compressed = comdb->file_compressed;

    struct DSVector *oseqs = ds_vector_create();
    struct DSVector *coarse_seq_links =
        get_coarse_sequence_links_at(links, coarse_links_index, id);
    struct fasta_seq *residues = cbp_coarse_read_fasta_seq(coarsedb, id);
    int fasta_length = strlen(residues->seq);
    int64_t *seq_lengths = cbp_compressed_get_lengths(comdb);

    int i = 0, j = 0;
    for (i = 0; i < coarse_seq_links->size; i++) {
        struct cbp_link_to_compressed *link =
            (struct cbp_link_to_compressed *)ds_vector_get(coarse_seq_links, i);
        if (link->coarse_start <= end && link->coarse_end >= start) {
            bool dir = link->dir;

            uint64_t original_start =
                get_max(0, (dir ? get_min(start + (link->original_start -
                                                   link->coarse_start),
                                          start + (link->original_end -
                                                   link->coarse_end)) :
                        get_min(link->original_start + link->coarse_end-end,
                                link->original_start - (end-link->coarse_end)))
                             - hit_pad_length);
            uint64_t original_end =
                get_min((dir ? get_max(end + (link->original_start -
                                              link->coarse_start),
                                       end + (link->original_end -
                                              link->coarse_end)) :
                     get_max(link->original_end - (start - link->coarse_start),
                             link->coarse_start + link->coarse_end - start))
                 + hit_pad_length, seq_lengths[link->org_seq_id]);

            fprintf(stderr, "%lu-%lu\n", original_start, original_end);

            struct cbp_compressed_seq *seq =
                       cbp_compressed_read_seq_at(comdb,link->org_seq_id);
            struct cbp_link_to_coarse *links_to_decompress = seq->links;
            for (; links_to_decompress;
                   links_to_decompress = links_to_decompress->next)
                if (links_to_decompress->original_end >= original_start &&
                    links_to_decompress->original_start <= original_end)
                    break;

            char *orig_str = malloc((original_end-original_start+1) *
                             sizeof(*orig_str));
            for (j = 0; j < original_end-original_start; orig_str[j++]='?');
            orig_str[original_end-original_start] = '\0';
 
            struct cbp_link_to_coarse *current = links_to_decompress;
            while (current && current->original_start <= original_end &&
                              current->original_end >= original_start) {
                pr_read_edit_script(orig_str, original_end-original_start+1,
                                    original_start, coarsedb, current);
                fprintf(stderr, "%s\n", orig_str);
                current = current -> next;
            }
            cbp_compressed_seq_free(seq);
        }
    }

    return NULL;
}
