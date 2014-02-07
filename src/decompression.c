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
#include "read_coarse.h"
#include "seq.h"

/*Takes in an entry in a compressed database and a coarse database and returns
  the decompressed sequence that the database entry came from as a pointer to
  a struct cbp_seq.*/
struct cbp_seq *cbp_decompress_seq(struct cbp_compressed_seq *cseq,
                                   struct cbp_coarse *coarsedb){
    uint64_t last_end = 0;
    int overlap;
    struct cbp_link_to_coarse *link;
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
        for (length = 0; chunk->seq[length] != '\0'; length++);/*length = strlen(chunk->seq);*/
        /*overlap represents the length of the overlap of the parts of the
          decompressed sequence that has been printed and the parts of the
          decompressed sequence currently being decompressed.*/
        overlap = last_end - link->original_start;
        dec_chunk = read_edit_script(link->diff, chunk->seq, length);
        /*Print all characters of the decompressed chunk past the index
          "overlap" unless overlap is greater than the length of the
          decompressed chunk.*/
        dec_chunk += overlap;
        if ((unsigned int)overlap < link->original_end - link->original_start) {
            int chunk_length = 0;
            char *section = NULL;

            for (i = 0; dec_chunk[chunk_length] != '\0'; chunk_length++);

            section = malloc((chunk_length + 1)*sizeof(*section));
            for (; i <= chunk_length; i++)
                section[i] = dec_chunk[i];
            ds_vector_append(decompressed_chunks, (void *)section);

            decompressed_length += chunk_length;
            dec_chunk -= overlap;
            free(dec_chunk);
            if (link->original_end > last_end)
                last_end = link->original_end;
        } else {
            dec_chunk -= overlap;
            free(dec_chunk);
        }
    }
    char *residues = malloc((decompressed_length+1)*sizeof(residues));
    for (i = 0; i < decompressed_chunks->size; i++) {
        char *current_chunk = (char *)ds_vector_get(decompressed_chunks, i);
        for (j = 0; current_chunk[j] != '\0'; j++)
            residues[copied++] = current_chunk[j];
    }
    residues[copied] = '\0';
    seq = cbp_seq_init(cseq->id, cseq->name, residues);
/*fprintf(stderr, "!!!!!\n");*/
    return seq;
}
