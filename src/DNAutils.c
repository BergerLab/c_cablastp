#include <stdlib.h>

char base_complement(char base){
  switch (base){
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      return '?';
  }
}

int bases_match(char a, char b, int dir_prod){
  if(dir_prod > 0)
    return (a == b && a != 'N') ? 1 : 0;
  else
    return a == base_complement(b) ? 1 : 0;
}

char *kmer_revcomp(char *kmer){
  int seed_size = 10; /*compress_flags.seed_size;*/
  char *kmer_revcomp = malloc(seed_size * sizeof(char));
  int i;
  for(i = 0; i < seed_size; i++)
    kmer_revcomp[i] = base_complement(kmer[seed_size-i-1]);
  return kmer_revcomp;
}
