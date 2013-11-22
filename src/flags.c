#include <stdint.h>

#include "flags.h"
#include "util.h"

struct opt_config *
load_compress_args()
{
    struct opt_config *conf;
    int32_t cpus;

    conf = opt_config_init();
    cpus = num_cpus();

    opt_flag_int(conf, 
        &compress_flags.gapped_window_size, "gapped-window-size", 75,
        "The size of the gapped match window.");
    opt_flag_int(conf,
        &compress_flags.ungapped_window_size, "ungapped-window-size", 10,
        "The size of the ungapped match window.");
    opt_flag_int(conf,
        &compress_flags.match_extend, "match-extend", 30,
        "The maximum number of residues to blindly extend a match without\n"
        "\tregard to sequence identity. This is to avoid small sequences\n"
        "\tin the coarse database.");
    opt_flag_int(conf,
        &compress_flags.match_kmer_size, "match-kmer-size", 10,
        "The size of the K-mer fragments to match in ungapped extension.");
    opt_flag_int(conf,
        &compress_flags.match_seq_id_threshold, "match-seq-id-threshold", 70,
        "The sequence identity threshold of an entire match.");
    opt_flag_int(conf,
        &compress_flags.min_match_len, "min-match-len", 50,
        "The minimum length of a match.");
    opt_flag_int(conf,
        &compress_flags.procs, "procs", cpus,
        "The number of total CPUs to use to divide work.");
    opt_flag_int(conf,
        &compress_flags.map_seed_size, "map-seed-size", 6,
        "The size of a seed in the K-mer map. This size combined with "
        "'ext-seed-size' forms the total seed size.");
    opt_flag_int(conf,
        &compress_flags.ext_seed_size, "ext-seed-size", 4,
        "The additional residues to require for each seed match.");
    opt_flag_int(conf,
        &compress_flags.ext_seq_id_threshold, "ext-seq-id-threshold", 50,
        "The sequence identity threshold of [un]gapped extension.");

    opt_flag_int(conf,
        &compress_flags.min_match_len, "min-match-len", 300,
        "The minimum size of an extension needed for a match to be recorded.");
    opt_flag_int(conf,
        &compress_flags.max_kmer_freq, "max-kmer-freq", 500,
        "The maximum number of entries for a k-mer in the seeds table.");
    opt_flag_int(conf,
        &compress_flags.overlap, "overlap", 100,
        "The maximum number of entries for a k-mer in the seeds table.");
    opt_flag_int(conf,
        &compress_flags.max_chunk_size, "max-chunk-size", 10000,
        "The maximum number of bases that are checked before adding a chunk without a match.");


    return conf;
}

