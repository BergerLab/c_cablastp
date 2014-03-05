#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "edit_scripts.h"
#include "read_compressed.h"
#include "read_coarse.h"
#include "seq.h"
#include "util.h"

int get_min(int a, int b){return a<b?a:b;}
int get_max(int a, int b){return a>b?a:b;}

struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length){
fprintf(stderr, "==========================================cbp_coarse_expand %d   %d-%d\n===================================================", id, start, end);
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
            for (j = 0; j < original_end-original_start+1; orig_str[j++]='?');
 
            struct cbp_link_to_coarse *current = links_to_decompress;
            while (current && current->original_start <= original_end &&
                              current->original_end >= original_start) {
                pr_read_edit_script(orig_str, original_end-original_start,
                                    original_start, coarsedb, current);
                fprintf(stderr, "%s\n", orig_str);
                current = current -> next;
            }
            cbp_compressed_seq_free(seq);
        }
    }

    return NULL;
}
