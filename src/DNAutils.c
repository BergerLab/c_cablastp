#include <stdlib.h>
#include "DNAutils.h"

char base_complement(char base){
    switch (base) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
    }
}

/*Checks if two bases match.  If one or both bases is an N, it is an
  automatic mismatch*/
int bases_match(char a, char b, int dir_prod){
    if (dir_prod > 0) /*Same-direction matches*/
        return (a == b && a != 'N') ? 1 : 0;
    else              /*Reverse-complement matches*/
        return a == base_complement(b) && a != 'N' ? 1 : 0;
}

char *kmer_revcomp(char *kmer){
    int seed_size = 10; /*compress_flags.seed_size;*/
    char *revcomp = malloc(seed_size * sizeof(char));
    int i;
    for (i = 0; i < seed_size; i++)
        revcomp[i] = base_complement(kmer[seed_size-i-1]);
    return revcomp;
}

char *string_revcomp(char *sequence, int length){
    char *revcomp;
    int i;
    /*Find the length of the sequence up to the null terminator if no length
      is given*/
    if (length < 0)
        for (length = 0; sequence[length] != '\0'; length++);
    revcomp = malloc((length+1) * sizeof(*revcomp));
    for (i = 0; i < length; i++)
        revcomp[i] = base_complement(sequence[length-i-1]);
    revcomp[length] = '\0';
    return revcomp;
}
