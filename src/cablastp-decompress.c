#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "opt.h"

#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "flags.h"
#include "fasta.h"
#include "read_compressed.h"
#include "seq.h"

static char *path_join(char *a, char *b)
{
    char *joined;

    joined = malloc((1 + strlen(a) + 1 + strlen(b)) * sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

int
main(int argc, char **argv)
{
    struct cbp_database *db;
    struct opt_config *conf;
    struct opt_args *args;
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct cbp_seq *org_seq;
    int i, org_seq_id;
    char *dir;
    struct timeval start, current;
    long double elapsed;
    conf = load_compress_args();
    args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir output-fasta-file\n",
            argv[0]);
        exit(1);
    }

    db = cbp_database_read(args->args[0], compress_flags.map_seed_size);

    org_seq_id = 0;
    gettimeofday(&start, NULL);
    char *compressed_filename = path_join(args->args[0], CABLASTP_COMPRESSED);
    FILE *compressed_file = fopen(compressed_filename, "r");
    struct cbp_compressed_seq **compressed = read_compressed(compressed_file);
    struct fasta_seq **coarse_sequences = malloc(10000*sizeof(*coarse_sequences));
    uint64_t num_coarse_sequences = 0;
    struct cbp_link_to_coarse *link;
    /*for (i = 0; compressed[i] != NULL; i++)
        for (link = (compressed[i])->links; link != NULL; link = link->next)
            printf("%ld %d %d\n%s\n\n\n", link->coarse_seq_id, link->coarse_start, link->coarse_end, link->diff);*/
    fsg = fasta_generator_start(path_join(args->args[0],CABLASTP_COARSE_FASTA),
                                FASTA_EXCLUDE_NCBI_BLOSUM62, 100);
    while (NULL != (seq = fasta_generator_next(fsg)))
        coarse_sequences[num_coarse_sequences++] = seq;

    fasta_generator_free(fsg);

    /*cbp_database_free(db);*/
    opt_config_free(conf);
    opt_args_free(args);

    return 0;
}
