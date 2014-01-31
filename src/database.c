#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "fasta.h"
#include "read_coarse.h"

static FILE * open_db_file(char *path, char *fopen_mode);

static char * path_join(char *a, char *b);

static char * basename(char *path);

struct cbp_database *
cbp_database_init(char *dir, int32_t seed_size, bool add)
{
    struct cbp_database *db;
    struct stat buf;
    FILE *ffasta, *fseeds, *flinks, *fcompressed, *findex;
    char *pfasta, *pseeds, *plinks, *pcompressed, *pindex;

    pfasta = path_join(dir, CABLASTP_COARSE_FASTA);
    pseeds = path_join(dir, CABLASTP_COARSE_SEEDS);
    plinks = path_join(dir, CABLASTP_COARSE_LINKS);
    pcompressed = path_join(dir, CABLASTP_COMPRESSED);
    pindex = path_join(dir, CABLASTP_INDEX);

    /* If we're not adding to a database, make sure `dir` does not exist. */
    if (!add && 0 == stat(dir, &buf)) {
        /* fprintf(stderr, */
            /* "The directory '%s' already exists. A new compressed " */
            /* "database cannot be created in the same directory as an " */
            /* "existing database. If you want to append to an existing " */
            /* "database, use the '--append' flag.\n", dir); */
        /* exit(1); */

        /* Just for testing purposes. */
        unlink(pfasta);
        unlink(pseeds);
        unlink(plinks);
        unlink(pcompressed);
        unlink(pindex);
        rmdir(dir);
    }
    /* Otherwise, check to make sure it *does* exist. */
    if (add && 0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for appending.", dir);
        exit(1);
    }

    if (0 != mkdir(dir, 0777)) {
        fprintf(stderr, "cbp_database_init: 'mkdir %s' failed because: %s\n",
            dir, strerror(errno));
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);
    db->name = basename(dir);

    ffasta = open_db_file(pfasta, "r+");
    fseeds = open_db_file(pseeds, "r+");
    flinks = open_db_file(plinks, "r+");
    fcompressed = open_db_file(pcompressed, "r+");
    findex = open_db_file(pindex, "r+");

    db->coarse_db = cbp_coarse_init(seed_size, ffasta, fseeds, flinks);
    db->com_db = cbp_compressed_init(fcompressed, findex);

    free(pfasta);
    free(pseeds);
    free(plinks);
    free(pcompressed);
    free(pindex);

    return db;
}

struct cbp_database *
cbp_database_read(char *dir, int32_t seed_size)
{
    struct cbp_database *db;
    struct stat buf;
    FILE *ffasta, *fseeds, *flinks, *fcompressed, *findex;
    char *pfasta  = path_join(dir, CABLASTP_COARSE_FASTA);
    char *pseeds = path_join(dir, CABLASTP_COARSE_SEEDS);
    char *plinks = path_join(dir, CABLASTP_COARSE_LINKS);
    char *pcompressed = path_join(dir, CABLASTP_COMPRESSED);
    char *pindex = path_join(dir, CABLASTP_INDEX);

    /* Make sure the database directory exists. */
    if (0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for reading.", dir);
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);
    db->name = basename(dir);

    ffasta = open_db_file(pfasta, "r");
    fseeds = open_db_file(pseeds, "r");
    flinks = open_db_file(plinks, "r");
    fcompressed = open_db_file(pcompressed, "r");
    findex = open_db_file(pindex, "r");

    db->coarse_db = cbp_coarse_init(seed_size, ffasta, fseeds, flinks);
    db->com_db = cbp_compressed_init(fcompressed, findex);
    cbp_database_populate(db, pfasta, plinks);

    return db;
}

void cbp_database_populate(struct cbp_database *db, const char *pfasta,
                           const char *plinks){
    struct fasta_seq_gen *fsg = fasta_generator_start(pfasta, "", 100);
    FILE *flinks = fopen(plinks, "r");
    read_coarse(db->coarse_db, flinks, fsg);
    fclose(flinks);
}

void
cbp_database_free(struct cbp_database *db)
{
    /* All files opened in cbp_database_init are close in subsequent frees. */
    cbp_coarse_free(db->coarse_db);
    cbp_compressed_free(db->com_db);
    free(db->name);
    free(db);
}

static FILE *
open_db_file(char *path, char *fopen_mode)
{
    struct stat buf;
    FILE *fp;

    if (0 != stat(path, &buf))
        if (-1 == creat(path, 0666)) {
            fprintf(stderr, "open_db_file: 'creat %s' failed: %s\n",
                path, strerror(errno));
            exit(1);
        }
    if (NULL == (fp = fopen(path, fopen_mode))) {
        fprintf(stderr, "open_db_file: 'fopen %s' failed: %s\n",
            path, strerror(errno));
        exit(1);
    }

    return fp;
}

static char *
path_join(char *a, char *b)
{
    char *joined;

    joined = malloc((1 + strlen(a) + 1 + strlen(b)) * sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

static char *
basename(char *path)
{
    char *base;
    int i;
    int len;

    len = strlen(path);
    for (i = len; i >= 0 && path[i] != '/'; i--);
    if (i > 0)
        i++;

    base = malloc((1 + len - i) * sizeof(*base));
    assert(base);

    strncpy(base, path + i, len - i);

    return base;
}
