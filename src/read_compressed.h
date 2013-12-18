#ifndef __CABLASTP_READ_COMPRESSED_H__
#define __CABLASTP_READ_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "align.h"
#include "bitpack.h"
#include "compressed.h"
#include "coarse.h"
#include "database.h"
#include "DNAutils.h"
#include "edit_scripts.h"

#include "stdbool.h"

struct links_to_coarse{
    char *header;
    struct cbp_link_to_coarse *links;
};

char *get_header(FILE *f);

struct cbp_link_to_coarse *read_link(FILE *f);

#endif
