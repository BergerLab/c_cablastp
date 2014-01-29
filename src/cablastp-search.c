#include <assert.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "ds.h"
#include "opt.h"

#include "blosum62.h"
#include "coarse.h"
#include "compression.h"
#include "database.h"
#include "flags.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"
#include "xml.h"

static char *path_join(char *a, char *b)
{
    char *joined;

    joined = malloc((1 + strlen(a) + 1 + strlen(b)) * sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}


/*Takes in the filename for a FASTA file and returns the file's sequences in a
  struct fasta_file*/
/*struct fasta_file *get_input_fasta(const char *filename){
    struct fasta_file *input_fasta_query = fasta_read_all(filename, "");
    return input_fasta_query;
}*/

/*Runs BLAST on the coarse database and stores the results in a temporary XML file*/
void blast_coarse(char *input_dir, char *query){
    char *input_path = path_join(input_dir, CABLASTP_COARSE_FASTA);
    char *blastn_command = "blastn -db  -outfmt 5 -query %s > CaBLAST_temp_blast_results.xml";
    int command_length = strlen(blastn_command) + strlen(input_path) + strlen(query) + 1;
    char *blastn = malloc(command_length * sizeof(*blastn));
    sprintf(blastn,
           "blastn -db %s -outfmt 5 -query %s > CaBLAST_temp_blast_results.xml",
           input_path, query);
    fprintf(stderr, "%s\n", blastn);
    system(blastn);
}
/*A function for traversing a parsed XML tree.  Takes in the root node, a
 *void * accumulator, and a function that takes in an xmlNode and a void *
 *accumulator and traverses the tree, applying the function on each node.
 */
void traverse_blast_xml(xmlNode *root, void (*f)(xmlNode *, void *), void *acc){
    for (; root; root = root->next) {
        f(root, acc);
        traverse_blast_xml(root->children, f, acc);
    }
}

/*Takes in the xmlNode representing a BLAST hit and populates a struct hit
  with its data.*/
struct hit *populate_blast_hit(xmlNode *node){
    struct hit *h = malloc(sizeof(*h));
    for (; node; node = node->next) {
        if (!strcmp((char *)node->name, "Hit_id")) {
            int name_length = strlen((char *)node->name) + 1;
            h->xml_name = malloc(name_length * sizeof(*(h->xml_name)));
            strcpy(h->xml_name, (char *)node->children->content);
        }
        if (!strcmp((char *)node->name, "Hit_num"))
            h->num = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hit_accession")) 
            h->accession = atoi((char *)node->children->content);
        /*if (!strcmp((char *)node->name, "Hit_accession")) 
            h->hsps = populate_blast_hsp(node->children);*/
    }
    return h;
}

/*Takes in an xmlNode and a linked list of xmlNodes converted to a void *
  and adds the node to the list if the node represents a BLAST hit.*/
void add_blast_hit(xmlNode *node, void *hits){
    if (!strcmp((char *)(node->name), "Hit"))
        ds_list_append((struct DSLinkedList *)hits,
                       (void *)populate_blast_hit(node->children));
}

/*Takes in an xmlNode and a linked list of xmlNodes converted to a void *
  and adds the node to the list if the node represents a BLAST hsp.*/
void add_blast_hsp(xmlNode *node, void *hsps){
    if (!strcmp((char *)(node->name), "Hsp"))
        ds_list_append((struct DSLinkedList *)hsps, (void *)node);
}

/*Takes in the xmlNode for a hit in a parsed BLAST XML tree and returns a
  linked list of all of the nodes in the tree that represent BLAST hits.*/
struct DSLinkedList *get_hit_hsps(xmlNode *hit){ 
    struct DSLinkedList *hsps = ds_list_create();
    traverse_blast_xml(hit, add_blast_hsp, hsps);
    return hsps;
}

/*Takes in the xmlNode for the root of a parsed BLAST XML tree and returns
  a linked list of all of the nodes in the tree that represent BLAST hits.*/
struct DSLinkedList *get_blast_hits(xmlNode *node){
    struct DSLinkedList *hits = ds_list_create();
    traverse_blast_xml(node, add_blast_hit, hits);
    return hits;
}

void expand_blast_hits(){
    xmlDoc *doc = NULL;
    xmlNode *root = NULL;
    doc = xmlReadFile("CaBLAST_temp_blast_results.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Could not parse CaBLAST_temp_blast_results.xml\n");
        return;
    }
    root = xmlDocGetRootElement(doc);
    struct DSLinkedList *hits = get_blast_hits(root);
    xmlFreeDoc(doc);
    xmlCleanupParser();
}

int
main(int argc, char **argv)
{ 
    /*struct cbp_database *db;
    struct cbp_compress_workers *workers;*/
    struct opt_config *conf;
    struct opt_args *args;
    /*struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct cbp_seq *org_seq;
    int i, org_seq_id;
    struct timeval start, current;
    long double elapsed;*/
    conf = load_search_args();
    args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir fasta-file [ --blast_args BLASTN_ARGUMENTS ]\n",
            argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }
    /*struct fasta_file *input_fasta_query = get_input_fasta(path_join(argv[1], CABLASTP_COARSE_FASTA));
    int i = 0;
    for(i = 0; i < input_fasta_query->length; i++){
        fprintf(stderr, "%s\n", input_fasta_query->seqs[i]->seq);
    }
    struct cbp_database *db = cbp_database_read(argv[1], search_flags.map_seed_size);*/
    blast_coarse(args->args[0], args->args[1]);
    expand_blast_hits();
/*
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
    /*cbp_coarse_save_binary(db->coarse_db);*//*
    cbp_coarse_save_seeds_plain(db->coarse_db);
    cbp_compressed_save_plain(db->com_db);
    /*cbp_compressed_save_binary(db->com_db);*//*

    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);
    cbp_compress_free_workers(workers);*/

    return 0;
}
