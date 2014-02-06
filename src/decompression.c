#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "stdbool.h"

#include "decompression.h"

struct cbp_seq *cbp_decompress_seq(struct cbp_compressed_seq *cseq,
                                   struct cbp_coarse *coarsedb){
    uint64_t num_coarse_sequences = 0;
    uint64_t last_end;
    int overlap;
    struct cbp_link_to_coarse *link;
    last_end = 0;
    int current_chunk = 0;
    printf("%s", compressed[i]->name);
    for (link = (compressed[i])->links; link != NULL; link = link->next) {
        /*overlap represents the length of the overlap of the parts of the
          decompressed sequence that has been printed and the parts of the
          decompressed sequence currently being decompressed.*/
        overlap = last_end - link->original_start;

        struct cbp_seq *chunk =
            cbp_seq_init_range(-1, "",
                               coarse_sequences[link->coarse_seq_id]->seq,
                               link->coarse_start, link->coarse_end);
        int length;
        for (length = 0; chunk->residues[length] != '\0'; length++);

        char *decompressed = read_edit_script(link->diff, chunk->residues,
                                                                  length);

        /*Print all characters of the decompressed chunk past the index
          "overlap" unless overlap is greater than the length of the
          decompressed chunk.*/
        decompressed += overlap;
        if (overlap < link->original_end - link->original_start)
            printf("%s", decompressed);
            decompressed -= overlap;
            free(decompressed);

            current_chunk++;
            if (link->original_end > last_end)
                last_end = link->original_end;

            cbp_seq_free(chunk);
        }

    return NULL;
}
