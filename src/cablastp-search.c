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


/*Takes in the filename for a FASTA file and returns the file's sequences in a
  struct fasta_file*/
/*struct fasta_file *get_input_fasta(const char *filename){
    struct fasta_file *input_fasta_query = fasta_read_all(filename, "");
    return input_fasta_query;
}*/

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
/*fprintf(stderr, "populate_blast_hsp\n");*/
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
            struct DSVector *seqs = cbp_coarse_expand(db->coarse_db,db->com_db,
                                        current_hit->accession,
                                        current_hsp->hit_from,
                                        current_hsp->hit_to);
            for (k = 0; k < seqs->size; k++) {
                struct cbp_seq *s = (struct cbp_seq *)ds_vector_get(seqs, k);
fprintf(stderr, "%d\n", s->id);
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

int
main(int argc, char **argv)
{ 
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
fprintf(stderr,"%d!!!!!!!!\n",expanded_hits->size);
    ds_vector_free(expanded_hits);
    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);

    return 0;
}
