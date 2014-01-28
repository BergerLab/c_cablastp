#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "opt.h"

#include "blosum62.h"
#include "coarse.h"
#include "compression.h"
#include "database.h"
#include "flags.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"

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
    struct cbp_compress_workers *workers;
    struct opt_config *conf;
    struct opt_args *args;
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct cbp_seq *org_seq;
    int i, org_seq_id;
    struct timeval start, current;
    long double elapsed;
    conf = load_compress_args();
    args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir fasta-file [ fasta-file ... ]\n",
            argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    db = cbp_database_init(args->args[0], compress_flags.map_seed_size, false);
    workers = cbp_compress_start_workers(db, compress_flags.procs);

    org_seq_id = 0;
    gettimeofday(&start, NULL);
    for (i = 1; i < args->nargs; i++) {
        fsg = fasta_generator_start(
            args->args[i], FASTA_EXCLUDE_NCBI_BLOSUM62, 100);

        while (NULL != (seq = fasta_generator_next(fsg))) {
            org_seq = cbp_seq_init(org_seq_id, seq->name, seq->seq);
            cbp_compress_send_job(workers, org_seq);

            fasta_free_seq(seq);

            org_seq_id++;
            if (org_seq_id % 1000 == 0) {
                gettimeofday(&current, NULL);
                elapsed = (long double)(current.tv_sec - start.tv_sec);
                printf("%d sequences compressed (%0.4Lf seqs/sec)\n",
                    org_seq_id, ((long double) org_seq_id) / elapsed);
            }
        }

        fasta_generator_free(fsg);
    }

    cbp_compress_join_workers(workers);
    cbp_coarse_save_plain(db->coarse_db);
    /*cbp_coarse_save_binary(db->coarse_db);*/
    cbp_coarse_save_seeds_plain(db->coarse_db);
    cbp_compressed_save_plain(db->com_db);
    /*cbp_compressed_save_binary(db->com_db);*/

    char *coarse_filename = path_join(args->args[0], "coarse.fasta");
    int len_filename = strlen(coarse_filename);
    int len_command = strlen("makeblastdb -dbtype nucl -in  -out") + len_filename + 1;
    char *makeblastdb = malloc(len_command * sizeof(makeblastdb));
    strcpy(makeblastdb, "makeblastdb -dbtype nucl -in ");
    strcat(makeblastdb, coarse_filename);
    strcat(makeblastdb, " -out ");
    strcat(makeblastdb, coarse_filename);

    system(makeblastdb);

    free(makeblastdb);
    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);
    cbp_compress_free_workers(workers);

    return 0;
}
