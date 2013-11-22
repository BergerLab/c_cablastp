#ifndef __CABLASTP_DNAUTILS_H__
#define __CABLASTP_DNAUTILS_H__

int bases_match(char a, char b, int dir_prod);
char base_complement(char base);
char *revcomp(char *kmer);

#endif
