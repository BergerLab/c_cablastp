#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "opt.h"

#include "blosum62.h"
#include "compression.h"
#include "flags.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"

/*This function will only be here temporarily*/
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
    fsg = fasta_generator_start(path_join(args->args[0],CABLASTP_COARSE_FASTA),
                                FASTA_EXCLUDE_NCBI_BLOSUM62, 100);
    while (NULL != (seq = fasta_generator_next(fsg))) {
        printf("%s\n%s\n", seq->name, seq->seq);
        /*org_seq = cbp_seq_init(org_seq_id, seq->name, seq->seq);
        cbp_compress_send_job(workers, org_seq);

        fasta_free_seq(seq);

        org_seq_id++;
        if (org_seq_id % 1000 == 0) {
            gettimeofday(&current, NULL);
            elapsed = (long double)(current.tv_sec - start.tv_sec);
            printf("%d sequences compressed (%0.4Lf seqs/sec)\n",
                   org_seq_id, ((long double) org_seq_id) / elapsed);
        }*/
    }

    fasta_generator_free(fsg);

    /*cbp_database_free(db);*/
    opt_config_free(conf);
    opt_args_free(args);

    return 0;
}
