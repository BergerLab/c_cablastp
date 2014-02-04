#ifndef __CABLASTP_DATABASE_H__
#define __CABLASTP_DATABASE_H__

#include <stdbool.h>

#include "coarse.h"
#include "compressed.h"

#define CABLASTP_COARSE_FASTA "coarse.fasta"
#define CABLASTP_COARSE_LINKS "coarse.links"
#define CABLASTP_COARSE_SEEDS "coarse.seeds"
#define CABLASTP_COARSE_LINKS_INDEX "coarse.links.index"
#define CABLASTP_COARSE_FASTA_INDEX "coarse.fasta.index"
#define CABLASTP_COMPRESSED "compressed.cbp"
#define CABLASTP_COMPRESSED_INDEX "compressed.cbp.index"

struct cbp_database {
    char *name;
    struct cbp_coarse *coarse_db;
    struct cbp_compressed *com_db;
};

struct cbp_database *
cbp_database_init(char *dir, int32_t seed_size, bool add);

struct cbp_database *
cbp_database_read(char *dir, int32_t seed_size);

void cbp_database_populate(struct cbp_database *db, const char *pfasta,
                           const char *plinks);

void
cbp_database_free(struct cbp_database *db);

#endif
