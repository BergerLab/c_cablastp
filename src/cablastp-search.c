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
#include "compressed.h"
#include "compression.h"
#include "database.h"
#include "fasta.h"
#include "flags.h"
#include "read_coarse.h"
#include "read_compressed.h"
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

/*Runs BLAST on the coarse database and stores the results in a temporary XML file*/
void blast_coarse(char *input_dir, char *query){
    char *input_path = path_join(input_dir, CABLASTP_COARSE_FASTA);
    char *blastn_command =
           "blastn -db  -outfmt 5 -query %s > CaBLAST_temp_blast_results.xml";
    int command_length = strlen(blastn_command) + strlen(input_path) +
                                                       strlen(query) + 1;
    char *blastn = malloc(command_length * sizeof(*blastn));
    sprintf(blastn,
           "blastn -db %s -outfmt 5 -query %s > CaBLAST_temp_blast_results.xml",
           input_path, query);
    fprintf(stderr, "%s\n", blastn);
    system(blastn);
    free(blastn);
}

/*Runs BLAST on the fine FASTA file*/
void blast_fine(char *subject, char *query){
    char *blastn_command =
           "blastn -subject  -query  > CaBLAST_results.txt";
    int command_length = strlen(blastn_command) + strlen(subject) +
                                                       strlen(query) + 1;
    char *blastn = malloc(command_length * sizeof(*blastn));
    sprintf(blastn,
            "blastn -subject %s -query %s > CaBLAST_results.txt",
            subject, query);
    fprintf(stderr, "%s\n", blastn);
    system(blastn);
    free(blastn);
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

/*Takes in the xmlNode representing a BLAST hsp and populates a struct hsp
  with its data.*/
struct hsp *populate_blast_hsp(xmlNode *node){
    struct hsp *h = malloc(sizeof(*h));
    h->xml_name = "Hsp";
    for (; node; node = node->next){
        if (!strcmp((char *)node->name, "Hsp_num"))
            h->num = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_evalue"))
            h->evalue = atof((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_query-from"))
            h->query_from = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_query-to"))
            h->query_to = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_hit-from"))
            h->hit_from = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hsp_hit-to"))
            h->hit_to = atoi((char *)node->children->content);
    }
    return h;
}

/*Takes in the xmlNode representing a BLAST hit and populates a struct hit
  with its data.*/
struct hit *populate_blast_hit(xmlNode *node){
    struct hit *h = malloc(sizeof(*h));
    h->xml_name = "Hit";
    for (; node; node = node->next) {
        if (!strcmp((char *)node->name, "Hit_num"))
            h->num = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hit_accession")) 
            h->accession = atoi((char *)node->children->content);
        if (!strcmp((char *)node->name, "Hit_hsps")){
            struct DSVector *hsps = ds_vector_create();
            xmlNode *hsp_node = node->children;
            for(; hsp_node; hsp_node = hsp_node->next)
                if(!strcmp((char *)hsp_node->name, "Hsp"))
                    ds_vector_append(hsps,
                           (void *)populate_blast_hsp(hsp_node->children));
            h->hsps = hsps;
        }
    }
    return h;
}

/*Takes in an xmlNode and a vector of xmlNodes converted to a void *
  and adds the node to the vector if the node represents a BLAST hit.*/
void add_blast_hit(xmlNode *node, void *hits){
    if (!strcmp((char *)(node->name), "Hit"))
        ds_vector_append((struct DSVector *)hits,
                         (void *)populate_blast_hit(node->children));
}

/*Takes in the xmlNode for the root of a parsed BLAST XML tree and returns
  a vector of all of the nodes in the tree that represent BLAST hits.*/
struct DSVector *get_blast_hits(xmlNode *node){
    struct DSVector *hits = ds_vector_create();
    traverse_blast_xml(node, add_blast_hit, hits);
    return hits;
}

struct DSVector *expand_blast_hits(struct cbp_database *db){
    int32_t i, j, k = 0;
    xmlDoc *doc = NULL;
    xmlNode *root = NULL;
    
    struct DSHashMap *used = ds_hashmap_create();
    struct DSVector *oseqs = ds_vector_create();

    doc = xmlReadFile("CaBLAST_temp_blast_results.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Could not parse CaBLAST_temp_blast_results.xml\n");
        ds_hashmap_free(used, true, true);
        ds_vector_free(oseqs);
        return NULL;
    }
    root = xmlDocGetRootElement(doc);
    struct DSVector *hits = get_blast_hits(root);
    for (i = 0; i < hits->size; i++) {
        struct hit *current_hit = (struct hit *)ds_vector_get(hits, i); 
        struct DSVector *hsps = current_hit->hsps;
        for (j = 0; j < hsps->size; j++) {
            struct hsp *current_hsp = (struct hsp *)ds_vector_get(hsps, j);
            if (current_hsp->evalue > search_flags.coarse_evalue)
                continue;
            struct DSVector *seqs = cbp_coarse_expand(db->coarse_db,db->com_db,
                                                      current_hit->accession,
                                                      current_hsp->hit_from,
                                                      current_hsp->hit_to);
            for (k = 0; k < seqs->size; k++) {
                struct cbp_seq *s = (struct cbp_seq *)ds_vector_get(seqs, k);
                if (ds_geti(used,s->id))
                    cbp_seq_free(s);
                else {
                    ds_vector_append(oseqs, (void *)s);
                    bool *t = malloc(sizeof(*t));
                    ds_puti(used, s->id, (void *)t);
                }
            }
            ds_vector_free_no_data(seqs);
        }
        ds_vector_free(hsps);
    }
    ds_vector_free(hits);
    
    xmlFreeDoc(doc);
    xmlCleanupParser();
    ds_hashmap_free(used, true, true);
    return oseqs;
}

void write_fine_fasta(struct DSVector *oseqs){
    int i;
    FILE *temp = fopen("CaBLAST_fine.fasta", "w");
    if (!temp) {
        fprintf(stderr, "Could not open CaBLAST_fine.fasta for writing\n");
        return;
    }
    for (i = 0; i < oseqs->size; i++) {
        struct cbp_seq *current_seq = (struct cbp_seq *)ds_vector_get(oseqs, i);
        fprintf(temp, "> %s\n%s\n", current_seq->name, current_seq->residues);
    }
    fclose(temp);
}

int
main(int argc, char **argv)
{
    int i = 0;
    struct cbp_database *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;
    conf = load_search_args();
    args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir fasta-file [ --blast_args BLASTN_ARGUMENTS ]\n",
            argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }
    db = cbp_database_read(args->args[0], search_flags.map_seed_size);

    blast_coarse(args->args[0], args->args[1]);
    struct DSVector *expanded_hits = expand_blast_hits(db);

    write_fine_fasta(expanded_hits);
    blast_fine("CaBLAST_fine.fasta", args->args[1]);

    for(i = 0; i < expanded_hits->size; i++)
        cbp_seq_free((struct cbp_seq *)ds_vector_get(expanded_hits, i));

    ds_vector_free_no_data(expanded_hits);
    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);
    /*system("rm CaBLAST_temp_blast_results.xml");
    system("rm CaBLAST_fine.fasta");*/

    return 0;
}
