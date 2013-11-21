#include <assert.h>
#include <stdbool.h>

#include "blosum62.h"

const int8_t BLOSUM62_RESIDUE_TO_INDEX[] = {
	0,  /* 'A' */
	1,  /* 'C' */
	2,  /* 'G' */
	3,  /* 'T' */
	4,  /* 'N' */
};
