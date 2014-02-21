#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "read_compressed.h"
#include "read_coarse.h"
#include "seq.h"
#include "util.h"

/*A function for getting the header for an entry in the coarse links file.
  Returns NULL if EOF is found before a newline.*/
char *get_coarse_header(FILE *f){
    int c = 0;
    char *header = malloc(30*sizeof(*header));
    int header_length = 30;
    int i = 0;
    while (c != EOF && c != '\n') {
        c = getc(f);
        if (c != EOF) {
            header[i] = (char)c;
            i++;
            if (i == header_length - 1) {
                header_length *= 2;
                header = realloc(header, header_length*sizeof(*header));
            }
        }
        else
            return NULL;
    }
    header[i] = '\0';
    header = realloc(header, (i+1)*sizeof(*header));
    return header;
}

/*Reads one link from a file with the links to the compressed database and
  converts its data to a struct cbp_link_to_compressed*/
struct cbp_link_to_compressed *read_coarse_link(FILE *f){
    struct cbp_link_to_compressed *link = malloc(sizeof(*link));
    link->org_seq_id = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->coarse_start = (uint16_t)read_int_from_file(2, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->coarse_end = (uint16_t)read_int_from_file(2, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_start = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_end = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->dir = getc(f) == '0';
    link->next = NULL;
    return link;
}

/*Takes in a pointer to the coarse.links file generated by cablastp-compress and
 *returns a vector containing all of the links in the coarse database entry for
 *the sequence the file pointer currently points to.  For this function to work
 *properly the file pointer must be pointing to the start of the header of a
 *sequence entry in the links file.
 */
struct DSVector *get_coarse_sequence_links(FILE *f){
    struct DSVector *links = ds_vector_create();
    char *h = get_coarse_header(f);
    if (h == NULL) {
        ds_vector_free(links);
        return NULL;
    }
    free(h);

    while (true) {
        char c = 1;
        struct cbp_link_to_compressed *current_link = read_coarse_link(f);
        if (current_link == NULL)
            break;
        ds_vector_append(links, (void *)current_link);
        c = getc(f);
        if (c == '#')
            break;
    }
    return links;
}

/*A wrapper function for get_coarse_sequence_links that handles fseek calls;
 *seeks to the index in the coarse.links file for the sequence at index id
 *and then calls get_coarse_sequence_links.  If fseek is successful, then
 *return the coarse sequence links for the coarse sequence at index id.
 *Otherwise, return NULL.
 */
struct DSVector *get_coarse_sequence_links_at(FILE *links, FILE *index,
                                                           int32_t id){
    int64_t offset = cbp_coarse_find_offset(index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(links, offset, SEEK_SET) == 0;
    if (!fseek_success) { 
        fprintf(stderr, "Error in seeking to offset %lu\n", offset);
        return NULL;
    }
    return get_coarse_sequence_links(links);
}

struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length){
fprintf(stderr, "cbp_coarse_expand %d   %d-%d\n", id, start, end);
    FILE *links = coarsedb->file_links;
    FILE *coarse_links_index = coarsedb->file_links_index;
    FILE *fasta = coarsedb->file_fasta;
    FILE *compressed = comdb->file_compressed;

    int i = 0;

    struct DSVector *oseqs = ds_vector_create();
    struct DSVector *coarse_seq_links =
        get_coarse_sequence_links_at(links, coarse_links_index, id);
    for (i = 0; i < coarse_seq_links->size; i++) {
        struct cbp_link_to_compressed *link =
            (struct cbp_link_to_compressed *)ds_vector_get(coarse_seq_links, i);
        if (link->coarse_start <= end && link->coarse_end >= start) {
            bool dir = link->dir;

            uint64_t original_start =
                fmax(0, (dir ? fmin(start + (link->original_start -
                                             link->coarse_start),
                                    start + (link->original_end -
                                             link->coarse_end)) :
                               fmin(link->original_start + link->coarse_end-end,
                                    link->original_start - 
                                          (end-link->coarse_end)))
                             - hit_pad_length);

            uint64_t original_end =
                (dir ? fmax(end + (link->original_start -
                                  link->coarse_start),
                           end + (link->original_end -
                                  link->coarse_end)) :
                       fmax(link->original_end - (start - link->coarse_start),
                            link->coarse_start + link->coarse_end - start))
                 + hit_pad_length;
        }
    }

    return NULL;
    /*Seek in the coarse links file to the links for the sequence being
      expanded*//*
    int64_t offset = cbp_coarse_find_offset(coarse_links_index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(links, offset, SEEK_SET) == 0;
    if (!fseek_success) { 
        fprintf(stderr, "Error in seeking to offset %lu\n", offset);
        return NULL;
    }

    struct DSHashMap *ids = ds_hashmap_create();
    struct DSVector *links_vector = get_coarse_sequence_links(links);
    struct DSVector *oseqs = ds_vector_create();

    int32_t links_count = links_vector->size;
    int32_t i = 0;
    for (; i < links_count; i++) {
        struct cbp_link_to_compressed *current_link =
            (struct cbp_link_to_compressed *)ds_vector_get(links_vector, i);
        if ((int16_t)start < current_link->coarse_start ||
            (int16_t)end > current_link->coarse_end)
            continue;
        if (ds_geti(ids, current_link->org_seq_id))
            continue;
        struct cbp_seq *oseq = cbp_compressed_read_seq(comdb, coarsedb,
                                                    current_link->org_seq_id);
        bool *t = malloc(sizeof(*t));
        *t = true;
        if (oseq != NULL)
            ds_vector_append(oseqs, oseq);
        ds_puti(ids, current_link->org_seq_id, t);
    }
    ds_hashmap_free(ids, false, true);
    ds_vector_free(links_vector);
    /*go_to_seq(id);*/
    /*return oseqs;*/
}

/*struct DSVector *
cbp_coarse_expand(struct cbp_coarse *coarsedb, struct cbp_compressed *comdb,
                  int32_t id, int32_t start, int32_t end){
    FILE *links = coarsedb->file_links;
    FILE *coarse_links_index = coarsedb->file_links_index;
    FILE *fasta = coarsedb->file_fasta;
    FILE *compressed = comdb->file_compressed;

    /*Seek in the coarse links file to the links for the sequence being
      expanded*//*
    int64_t offset = cbp_coarse_find_offset(coarse_links_index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(links, offset, SEEK_SET) == 0;
    if (!fseek_success) { 
        fprintf(stderr, "Error in seeking to offset %lu\n", offset);
        return NULL;
    }

    struct DSHashMap *ids = ds_hashmap_create();
    struct DSVector *links_vector = get_coarse_sequence_links(links);
    struct DSVector *oseqs = ds_vector_create();

    int32_t links_count = links_vector->size;
    int32_t i = 0;
    for (; i < links_count; i++) {
        struct cbp_link_to_compressed *current_link =
            (struct cbp_link_to_compressed *)ds_vector_get(links_vector, i);
        if ((int16_t)start < current_link->coarse_start ||
            (int16_t)end > current_link->coarse_end)
            continue;
        if (ds_geti(ids, current_link->org_seq_id))
            continue;
        struct cbp_seq *oseq = cbp_compressed_read_seq(comdb, coarsedb,
                                                    current_link->org_seq_id);
        bool *t = malloc(sizeof(*t));
        *t = true;
        if (oseq != NULL)
            ds_vector_append(oseqs, oseq);
        ds_puti(ids, current_link->org_seq_id, t);
    }
    ds_hashmap_free(ids, false, true);
    ds_vector_free(links_vector);
    /*go_to_seq(id);*//*
    return oseqs;
}*/

/*Takes in an index file from a coarse database and the ID number of the
  sequence in the corresponding database file that the user wants to seek to
  and returns the byte offset of the sequence in its database file.*/
int64_t cbp_coarse_find_offset(FILE *index_file, int id){
    int i;
    int try_off = id * 8;
    bool fseek_success = fseek(index_file, try_off, SEEK_SET) == 0;
    int64_t mask = make_mask(8);
    int64_t offset = (int64_t)(-1);
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d\n", try_off);
        return (int64_t)(-1);
    }
    for (i = 0; i < 8; i++) {
        int64_t current_byte=((int64_t)(getc(index_file))&mask);
        offset <<= 8;
        offset |= current_byte;
    }
    return offset;
}

/*Takes in as arguments a coarse database and the ID number of the sequence in
 *the coarse FASTA file to read in and gets a struct fasta_seq for that
 *sequence.
 */
struct fasta_seq *cbp_coarse_read_fasta_seq(struct cbp_coarse *coarsedb, int id){
    int64_t offset = cbp_coarse_find_offset(coarsedb->file_fasta_index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(coarsedb->file_fasta, offset, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d\n", offset);
        return NULL;
    }
    return fasta_read_next(coarsedb->file_fasta, "");
}
