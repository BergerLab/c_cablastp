#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "decompression.h"

/*A function for getting the header for an entry in the compressed links file.
  Returns NULL if EOF is found before a newline.*/
char *get_header(FILE *f){
    int c = 0;
    char *header = malloc(10000*sizeof(*header));
    int header_length = 10000;
    int i = 0;
    while(c != EOF && c != '\n'){
        c = getc(f);
        if(c != EOF){
            header[i] = (char)c;
            i++;
            if(i == header_length-1){
                header_length *= 2;
                realloc(header, header_length*sizeof(*header));
            }
        }
        else
            return NULL;
    }
    header[i] = '\0';
    realloc(header, (i+1)*sizeof(*header));
    return header;
}

/*Reads one link from a file with the links to the coarse database and converts
  its data to a struct cbp_link_to_coarse*/
struct cbp_link_to_coarse *read_link(FILE *f){
    int i;
    int c = 0;
    struct cbp_link_to_coarse *link = malloc(sizeof(*link));
    uint64_t coarse_seq_id = (uint64_t)0;
    int16_t coarse_start = (int16_t)0;
    int16_t coarse_end = (int16_t)0;
    int16_t script_length;
    int chars_to_read;
    char *half_bytes;
    char indices[6];

    for(i = 0; i < 8; i++){
        c = getc(f);
        if(c == EOF){
            free(link);
            return NULL;
        }
        coarse_seq_id |= (uint64_t)c;
        coarse_seq_id <<= 8;
    }

    for(i = 0; i < 6; i++){
        c = getc(f);
        if(c == EOF){
            free(link);
            return NULL;
        }
        indices[i] = (char)c;
    }
    coarse_start |= (int16_t)indices[0];
    coarse_start <<= 8;
    coarse_start |= (int16_t)indices[1];
    coarse_end |= (int16_t)indices[2];
    coarse_end <<= 8;
    coarse_end |= (int16_t)indices[3];
    script_length |= (int16_t)indices[4];
    script_length <<= 8;
    script_length |= (int16_t)indices[5];

    chars_to_read = script_length / 2;
    if(script_length % 2 == 1)
        chars_to_read++;
    half_bytes = malloc(chars_to_read*sizeof(*half_bytes));
    for(i = 0; i < chars_to_read; i++){
        c = getc(f);
        if(c == EOF){
            free(link);
            return NULL;
        }
        half_bytes[i] = (char)c;
    }
    edit_script = half_bytes_to_ASCII(half_bytes, script_length);

    link->diff = edit_script;
    link->coarse_seq_id = coarse_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    return link;
}
