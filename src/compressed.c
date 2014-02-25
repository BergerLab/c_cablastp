#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "bitpack.h"
#include "compressed.h"
#include "edit_scripts.h"
#include "read_compressed.h"

struct cbp_compressed *
cbp_compressed_init(FILE *file_compressed, FILE *file_index)
{
    struct cbp_compressed *com_db;

    com_db = malloc(sizeof(*com_db));
    assert(com_db);

    com_db->file_compressed = file_compressed;
    com_db->file_index = file_index;
    com_db->seqs = ds_vector_create_capacity(100);

    return com_db;
}

void
cbp_compressed_free(struct cbp_compressed *com_db)
{
    int i;

    fclose(com_db->file_compressed);
    fclose(com_db->file_index);

    for (i = 0; i < com_db->seqs->size; i++)
        cbp_compressed_seq_free(cbp_compressed_seq_at(com_db, i));

    ds_vector_free_no_data(com_db->seqs);
    free(com_db);
}

int32_t
cbp_compressed_size(struct cbp_compressed *com_db)
{
    return com_db->seqs->size;
}

void
cbp_compressed_add(struct cbp_compressed *com_db,
                   struct cbp_compressed_seq *seq)
{
    ds_vector_append(com_db->seqs, (void*) seq);
}

/*Takes in as input the compressed database and converts its sequences to a
 *binary output format, which is printed to the file pointed to by
 *com_db->file_compressed.
 */
void
cbp_compressed_save_binary(struct cbp_compressed *com_db)
{
    int i;
    struct cbp_compressed_seq *seq;
    struct cbp_link_to_coarse *link;

    int16_t mask = (((int16_t)1)<<8)-1;
    uint64_t index = (uint64_t)0;

    for (i = 0; i < com_db->seqs->size; i++) {
        int j;
        char *id_bytes;

        output_int_to_file(index, 8, com_db->file_index);
        seq = cbp_compressed_seq_at(com_db, i);
        id_bytes = read_int_to_bytes(seq->id, 8);

        putc('>', com_db->file_compressed);
        putc(' ', com_db->file_compressed);
        for (j = 0; j < 8; j++)
            putc(id_bytes[j], com_db->file_compressed);
        putc(';', com_db->file_compressed);
        putc(' ', com_db->file_compressed);

        index += 12;

        for(j = 0; seq->name[j] != '\0'; j++)
            putc(seq->name[j], com_db->file_compressed);
        index += j;
        putc('\n', com_db->file_compressed);
        index++;

        struct cbp_link_to_coarse *find_length = seq->links;
        uint64_t original_length = 0;
        for (; find_length; find_length = find_length->next)
            original_length = find_length->original_end + 1;

        output_int_to_file(original_length, 8, com_db->file_compressed);

        for (link = seq->links; link != NULL; link = link->next){
            /*Convert the start and end indices for the link to two
              characters.*/
            uint64_t org_start = link->original_start;
            uint64_t org_end   = link->original_end;
            int16_t cor_start = (int16_t)link->coarse_start;
            int16_t cor_end   = (int16_t)link->coarse_end;
            char cor_start_left  = (cor_start >> 8) & mask;
            char cor_start_right = cor_start & mask;
            char cor_end_left    = (cor_end >> 8) & mask;
            char cor_end_right   = cor_end & mask;

            char script_left, script_right;
            char *edit_script = link->diff;
            char *script;
            int16_t script_length = (int16_t)0;
            int j;

            /*Output the ID of the current chunk as 8 characters*/
            int coarse_seq_id = link->coarse_seq_id;
            for (j = 7; j >= 0; j--) {
                char b = shift_right(coarse_seq_id, j*8) & mask;
                putc(b, com_db->file_compressed);
            }
            index += 8;

            /*Represent the length of the edit script as two characters and get
              the edit script as a sequence of half-bytes*/
            while (edit_script[script_length] != '\0')
                script_length++;

            script_left  = (script_length >> 8) & mask;
            script_right = script_length & mask;
            script = edit_script_to_half_bytes(edit_script);

            /*Output the indices of the start and end of the sequence being
             *linked to and the length of the edit script represented in 16
             *bits, and the length of the edit script.
             */
            output_int_to_file(org_start, 8, com_db->file_compressed);
            output_int_to_file(org_end, 8, com_db->file_compressed);
            putc(cor_start_left, com_db->file_compressed);
            putc(cor_start_right, com_db->file_compressed);
            putc(cor_end_left, com_db->file_compressed);
            putc(cor_end_right, com_db->file_compressed);
            putc(script_left, com_db->file_compressed);
            putc(script_right, com_db->file_compressed);
            index += 6;

            /*Output all of the characters of the edit script as half-bytes*/
            for (j = 0; j < script_length/2; j++)
                putc(script[j], com_db->file_compressed);
            index += j;


            /*If there are more links for this sequence, the character after
             *the edit script is a space.  Otherwise, the character after the
             *edit script is the > of the FASTA header.  If we are printing the
             *last link of the last sequence, print a newline.
             */
            if (link->next) {
                putc(' ', com_db->file_compressed);
                index++;
            }
            else
                if (i == com_db->seqs->size-1) {
                    putc('\n', com_db->file_compressed);
                    index++;
                }
        }
    }
}

void
cbp_compressed_save_plain(struct cbp_compressed *com_db)
{
    int i;
    struct cbp_compressed_seq *seq;
    struct cbp_link_to_coarse *link;

    for (i = 0; i < com_db->seqs->size; i++) {
        seq = cbp_compressed_seq_at(com_db, i);
        fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(com_db->file_compressed,
                "reference sequence id: %d, reference range: (%d, %d), original sequence range: (%ld %ld)\n%s\n",
                link->coarse_seq_id, link->coarse_start, link->coarse_end, link->original_start, link->original_end,
                link->diff);
    }
}

void
cbp_compressed_write(struct cbp_compressed *com_db,
                     struct cbp_compressed_seq *seq)
{
    struct cbp_link_to_coarse *link;

    fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
    for (link = seq->links; link != NULL; link = link->next) {
        fprintf(com_db->file_compressed,
            "reference sequence id: %d, reference range: (%d, %d)\n%s\n",
            link->coarse_seq_id, link->coarse_start, link->coarse_end,
            link->diff);
    }
}

/*Outputs a compressed sequence in the compressed database to the database's
  compressed file in binary format.*/
void
cbp_compressed_write_binary(struct cbp_compressed *com_db,
                            struct cbp_compressed_seq *seq)
{
    int i;
    struct cbp_link_to_coarse *link;
    int16_t mask = (((int16_t)1)<<8)-1;
    char *id_string = malloc(20*sizeof(*id_string));
    uint64_t index = ftell(com_db->file_compressed);
    output_int_to_file(index, 8, com_db->file_index);
    sprintf(id_string, "%ld", seq->id);

    /*Output the header for the sequence*/
    putc('>', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; id_string[i] != '\0'; i++)
        putc(id_string[i], com_db->file_compressed);
    putc(';', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; seq->name[i] != '\0'; i++)
        putc(seq->name[i], com_db->file_compressed);
    putc('\n', com_db->file_compressed);

    free(id_string);

    struct cbp_link_to_coarse *find_length = seq->links;
    uint64_t original_length = 0;
    for (; find_length; find_length = find_length->next)
        original_length = find_length->original_end + 1;

    output_int_to_file(original_length, 8, com_db->file_compressed);

    for (link = seq->links; link != NULL; link = link->next){
        /*Convert the start and end indices for the link to two
          characters.*/
        uint64_t org_start = link->original_start;
        uint64_t org_end   = link->original_end;
        int16_t cor_start = (int16_t)link->coarse_start;
        int16_t cor_end   = (int16_t)link->coarse_end;
        char cor_start_left  = (cor_start >> 8) & mask;
        char cor_start_right = cor_start & mask;
        char cor_end_left    = (cor_end >> 8) & mask;
        char cor_end_right   = cor_end & mask;

        char script_left, script_right;
        char *edit_script = link->diff;
        char *script = edit_script_to_half_bytes(edit_script);
        int16_t script_length = (int16_t)0;
        int odd;

        /*Output the ID of the current chunk as 8 characters*/
        uint64_t coarse_seq_id = link->coarse_seq_id;
        for (i = 7; i >= 0; i--) {
            char b = (char)(shift_right(coarse_seq_id, i*8) & mask);
            putc(b, com_db->file_compressed);
        }

        /*Represent the length of the edit script as two characters and get
          the edit script as a sequence of half-bytes*/
        while (edit_script[script_length] != '\0')
            script_length++;
        odd = script_length % 2 == 1 ? 1 : 0;
        script_left  = (script_length >> 8) & mask;
        script_right = script_length & mask;

        /*Output the indices of the start and end of the sequence being
         *linked to and the length of the edit script represented in 16
         *bits, and the length of the edit script.
         */
        output_int_to_file(org_start, 8, com_db->file_compressed);
        output_int_to_file(org_end, 8, com_db->file_compressed);
        putc(cor_start_left, com_db->file_compressed);
        putc(cor_start_right, com_db->file_compressed);
        putc(cor_end_left, com_db->file_compressed);
        putc(cor_end_right, com_db->file_compressed);
        putc(script_left, com_db->file_compressed);
        putc(script_right, com_db->file_compressed);

        /*Output all of the characters of the edit script as half-bytes*/
        for (i = 0; i < script_length/2+odd; i++)
            putc(script[i], com_db->file_compressed);

        /*If there are more links for this sequence, the character after
         *the edit script is a space.  Otherwise, the character after the
         *edit script is the > of the FASTA header.  If we are printing the
         *last link of the last sequence, print a newline.
         */
        if (link->next)
            putc(' ', com_db->file_compressed);

        free(script);
    }
    putc('\n', com_db->file_compressed);
}

struct cbp_compressed_seq *
cbp_compressed_seq_init(int32_t id, char *name)
{
    struct cbp_compressed_seq *seq;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->links = NULL;

    seq->name = malloc((1 + strlen(name)) * sizeof(*seq->name));
    assert(seq->name);
    strcpy(seq->name, name);

    return seq;
}

void
cbp_compressed_seq_free(struct cbp_compressed_seq *seq)
{
    struct cbp_link_to_coarse *link1, *link2;

    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cbp_link_to_coarse_free(link1);
        link1 = link2;
    }
    free(seq->name);
    free(seq);
}

void
cbp_compressed_seq_addlink(struct cbp_compressed_seq *seq,
                           struct cbp_link_to_coarse *newlink)
{
    struct cbp_link_to_coarse *link;

    assert(newlink->next == NULL);

    if (seq->links == NULL) {
        seq->links = newlink;
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
}


struct cbp_link_to_coarse *
cbp_link_to_coarse_init(int32_t coarse_seq_id,
                        uint64_t original_start, uint64_t original_end,
                        uint16_t coarse_start, uint16_t coarse_end,
                        struct cbp_alignment alignment, bool dir){
    struct cbp_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->coarse_seq_id = coarse_seq_id;
    link->original_start = original_start;
    link->original_end = original_end;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    link->diff = make_edit_script(alignment.org, alignment.ref, dir, alignment.length);
    assert(link->diff);

    return link;
}

struct cbp_link_to_coarse *
cbp_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               uint64_t original_start, uint64_t original_end,
                               uint16_t coarse_start, uint16_t coarse_end,
                               bool dir){
    struct cbp_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->diff = malloc(2*sizeof(char));
    link->diff[0] = dir ? '0' : '1';
    link->diff[1] = '\0';
    link->coarse_seq_id = coarse_seq_id;
    link->original_start = original_start;
    link->original_end = original_end;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    return link;
}

void
cbp_link_to_coarse_free(struct cbp_link_to_coarse *link)
{
    free(link->diff);
    free(link);
}

struct cbp_compressed_seq *
cbp_compressed_seq_at(struct cbp_compressed *com_db, int32_t i)
{
    return (struct cbp_compressed_seq *) ds_vector_get(com_db->seqs, i);
}
