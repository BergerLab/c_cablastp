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
#include "DNAutils.h"
#include "edit_scripts.h"
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
    fsg = fasta_generator_start(path_join(args->args[0],CABLASTP_COARSE_FASTA),
                                FASTA_EXCLUDE_NCBI_BLOSUM62, 100);
    while (NULL != (seq = fasta_generator_next(fsg)))
        coarse_sequences[num_coarse_sequences++] = seq;
    for (i = 0; compressed[i] != NULL; i++) {
        printf("%s", compressed[i]->name);
        int current_chunk = 0;
        bool prev_match = false;
        for (link = (compressed[i])->links; link != NULL; link = link->next) {
            /*If link -> diff[1] is a null terminator, this means this link is
              to a chunk added without any matches*/
            int start = link->diff[1] == '\0' ? 100 : 0;
            if (!prev_match && start == 100)
                start--;
            struct cbp_seq *chunk =
                cbp_seq_init_range(-1, "",
                                   coarse_sequences[link->coarse_seq_id]->seq,
                                   link->coarse_start, link->coarse_end);
            if(link->diff[0] == '1'){
                /*char *revcomp = string_revcomp(chunk->residues, -1);
                free(chunk->residues);
                chunk->residues = revcomp;*/
fprintf(stderr, "\n\n\n%s\n\n\n", chunk->residues);
            }
            int length;
            for (length = 0; chunk->residues[length] != '\0'; length++);
            if (start == 0 || current_chunk == 0)
                printf("%s", read_edit_script(link->diff, chunk->residues, length));
            else
                printf("%s", read_edit_script(link->diff, chunk->residues+start, length-start));
            prev_match = start == 0;
            current_chunk++;
if(current_chunk == 3)break;
        }
    }
    putc('\n', stdout);

    fasta_generator_free(fsg);

    /*cbp_database_free(db);*/
    opt_config_free(conf);
    opt_args_free(args);

fprintf(stderr, "%s\n", read_edit_script("1s0--GAs6As2-s2TGi3Ts1Gs4Ci2Ts2-s2-i2Gs1Gi2Ci2Gs1Ts4Ti4Cs1Gi2CCi3Gi2CGi3TTs3Gi3As2Ai2Gs1Cs3Gs2Cs4-i2Ci2Gs1Ti3TCi5Gs2Ai3Gs2Ci2Gs3-i2Ci3CCGGi6Cs1Ci6GGi4Ts2-Cs3CGi3Ts3-s3Cs4-s4ATs4-s2TTs3Gs2Ci4TTs2Ti5As4-Cs5Ts3Ai7Ti2Gi6GCi4CCs3-Cs3Cs2Ci2CCs7-Gi3Ai2GCi3Gs5-i2Ai3Gs6--s4-s6ATCi5Cs1Ts2Ts3ATi4As3CCTCCs6Ts2Gs3-s2Ts3-GCs4Ai3Ci2Cs3-i2TTTTs5-s2As3CCs4--i4Cs3Ci6Cs4--s11-s3--s3-s3AAs3Gs3CCs4Gi4Ts1Cs5--i5Cs1Gi2Cs1Ts3As2As2--TTs11-s2-s5TTTs4Cs4Cs2Ci2As4-s2-s3ACs5---i5Ti2Ci2CTs3Ts2Ci4Ai2GGs3Ai10Gs1Ci2Ai2As2-i3ACi6Cs2A","CTTCAAGGAACAGCTTCGCGATATCGACCTTCTGGTGATCGACGACATGCAGTTCCTGCAGGGCAAGTCGATCCAGCACGAGTTCTGCCATCTTCTGAATACGCTGCTCGACAGCGCAAAACAGGTCGTGGTTGCCGCCGACCGCGCCCCTTCGGAACTGGAATCGCTCGACGTTCGCGTCCGTTCCCGCCTTCAGGGCGGTGTCGCTCTGGAAGTGGCCGCACCCGATTATGAAATGCGCCTGGAAATGCTGCGTCGCCGTCTGGCCAGCGCACAGTGCGAGGATGCAAGCCTCGATATCGGTGAAGAGATTCTCGCCCATGTCGCGCGCACGGTCACGGGTTCCGGCCGTGAACTGGAAGGCGCATTCAACCAGCTTCTGTTCCGCCAGTCTTTCGAGCCGAATATTTCCATCGACCGCGTGGACGAGTTGCTCGGCC",440));

    return 0;
}
