#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "read_compressed.h"

/*A function for getting the header for an entry in the compressed links file.
  Returns NULL if EOF is found before a newline.*/
char *get_compressed_header(FILE *f){
    int c = 0;
    char *header = malloc(10000*sizeof(*header));
    int header_length = 10000;
    int i = 0;
    while (c != EOF && c != '\n'){
        c = getc(f);
        if (c != EOF){
            header[i] = (char)c;
            i++;
            if (i == header_length-1){
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

/*Reads one link from a file with the links to the coarse database and converts
  its data to a struct cbp_link_to_coarse*/
struct cbp_link_to_coarse *read_compressed_link(FILE *f){
    int i;
    unsigned int c = 0;
    struct cbp_link_to_coarse *link = malloc(sizeof(*link));
    uint64_t coarse_seq_id = (uint64_t)0;
    uint64_t original_start = (uint64_t)0;
    uint64_t original_end = (uint64_t)0;
    uint16_t coarse_start = (uint16_t)0;
    uint16_t coarse_end = (uint16_t)0;
    uint16_t script_length = (uint16_t)0;
    int chars_to_read;
    char *half_bytes;
    char *diff;
    unsigned char indices[22];

    for (i = 0; i < 8; i++){
        c = getc(f);
        if (feof(f)) {
            free(link);
            return NULL;
        }
        coarse_seq_id <<= 8;
        coarse_seq_id |= (uint64_t)c;
    }
    for (i = 0; i < 22; i++){
        c = getc(f);
        if (feof(f)) {
            free(link);
            return NULL;
        }
        indices[i] = (char)c;
    }

    for (i = 0; i < 8; i++) {
        original_start <<= 8;
        original_start |= (uint64_t)indices[i];
        original_end <<= 8;
        original_end |= (uint64_t)indices[8+i];
    }
    coarse_start |= (uint16_t)indices[16];
    coarse_start <<= 8;
    coarse_start |= (uint16_t)indices[17];
    coarse_end |= (uint16_t)indices[18];
    coarse_end <<= 8;
    coarse_end |= (uint16_t)indices[19];
    script_length |= (uint16_t)indices[20];
    script_length <<= 8;
    script_length |= (uint16_t)indices[21];

    chars_to_read = script_length / 2;
    if (script_length % 2 == 1)
        chars_to_read++;
    half_bytes = malloc(chars_to_read*sizeof(*half_bytes));
    for (i = 0; i < chars_to_read; i++) {
        c = getc(f);
        if (feof(f)) {
            free(link);
            return NULL;
        }
        half_bytes[i] = (char)c;
    }
    diff = half_bytes_to_ASCII(half_bytes, script_length);
    link->diff = diff;
    link->coarse_seq_id = coarse_seq_id;
    link->original_start = original_start;
    link->original_end = original_end;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    return link;
}

/*Takes the compressed.cbp generated by cablastp-compress and parses it to get
  an array of compressed sequences.*/
struct cbp_compressed_seq **read_compressed(FILE *f){
    int length = 0;
    struct cbp_compressed_seq **compressed_seqs =
        malloc(1000*sizeof(*compressed_seqs));
    /*Read each sequence*/
    while (true) {
        struct cbp_link_to_coarse *links = NULL;
        char *header = get_compressed_header(f);
        if (header == NULL)
            break;
        /*Read each link in the sequence*/
        while (true) {
            char c = 1;
            struct cbp_link_to_coarse *current_link = read_compressed_link(f);
            if (current_link == NULL)
                break;
            current_link->next = links;
            if (current_link->next == NULL)
                links = current_link;
            else {
                while (current_link->next->next != NULL)
                    current_link->next = current_link->next->next;
                current_link->next->next = current_link;
                current_link->next = NULL;
            }
            c = getc(f);
            /*If we find a newline at the end of the link we are at the end of
              the sequence.*/
            if (c == '\n')
                break;
        }
        compressed_seqs[length] = malloc(sizeof(*(compressed_seqs[length])));
        (compressed_seqs[length])->name = header;
        (compressed_seqs[length]->links) = links;
        length++;
    }

    compressed_seqs = realloc(compressed_seqs,
                              (length+1)*sizeof(*compressed_seqs));
    compressed_seqs[length] = NULL;
    return compressed_seqs;
}
