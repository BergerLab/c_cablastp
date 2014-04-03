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

#include "bitpack.h"
#include "coarse.h"
#include "compressed.h"
#include "compression.h"
#include "database.h"
#include "decompression.h"
#include "DNAalphabet.h"
#include "fasta.h"
#include "flags.h"
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

/*Runs BLAST on the coarse database and stores the results in a temporary
  XML file.*/
void blast_coarse(struct opt_args *args, uint64_t dbsize){
    char *input_path = path_join(args->args[0], CABLAST_COARSE_FASTA);
    char *blastn_command =
           "blastn -db  -outfmt 5 -query  -dbsize  -task blastn -evalue "
           "%s > CaBLAST_temp_blast_results.xml";
    int command_length = strlen(blastn_command) + strlen(input_path) +
                                                  strlen(args->args[1]) + 31;
    char *blastn = malloc(command_length * sizeof(*blastn));
    sprintf(blastn,"blastn -db %s -outfmt 5 -query %s -dbsize %lu -task blastn"
                   " -evalue %s > CaBLAST_temp_blast_results.xml",
           input_path, args->args[1], dbsize, search_flags.coarse_evalue);
    fprintf(stderr, "%s\n", blastn);
    system(blastn);
    free(blastn);
}

/*Runs BLAST on the fine FASTA file*/
void blast_fine(char *subject, struct opt_args *args, uint64_t dbsize){
    char *blastn_command =
           "blastn -subject  -query -outfmt 5 -dbsize  -task blastn "
           "> CaBLAST_results.txt";
    int command_length = strlen(blastn_command) + strlen(subject) +
                                                  strlen(args->args[1]) + 31;
    char *blastn = malloc(command_length * sizeof(*blastn));
    sprintf(blastn,
            "blastn -subject %s -query %s -outfmt 5 -dbsize %lu -task blastn "
            "> CaBLAST_results.txt",
            subject, args->args[1], dbsize);
    fprintf(stderr, "%s\n", blastn);
    system(blastn);
    free(blastn);
}


/*A function for traversing a parsed XML tree.  Takes in the root node, a
 *void * accumulator, and a function that takes in an xmlNode and a void *
 *accumulator and traverses the tree, applying the function on each node.
 */
void traverse_blast_xml(xmlNode *root, void (*f)(xmlNode *,void *), void *acc){
    for (; root; root = root->next) {
        f(root, acc);
        traverse_blast_xml(root->children, f, acc);
    }
}

/*A wrapper function for traverse_blast_xml that treats the xmlNode passed into
  r as the root of the XML tree, skipping any sibling nodes r has*/
void traverse_blast_xml_r(xmlNode *r, void (*f)(xmlNode *, void *), void *acc){
    f(r, acc);
    traverse_blast_xml(r->children, f, acc);
}

/*Takes in the xmlNode representing a BLAST hsp and populates a struct hsp
  with its data.*/
struct hsp *populate_blast_hsp(xmlNode *node){
    struct hsp *h = malloc(sizeof(*h));
    bool hit_frame_fwd = true;
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
        if (!strcmp((char *)node->name, "Hsp_hit-frame"))
            hit_frame_fwd = atoi((char *)node->children->content) == 1;
    }
    if (!hit_frame_fwd) {
        int hit_to = h->hit_from;
        h->hit_from = h->hit_to;
        h->hit_to = hit_to;
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
    traverse_blast_xml_r(node, add_blast_hit, hits);
    return hits;
}

/*Takes in an xmlNode and a vector of xmlNodes converted to a void *
  and adds the node to the vector if the node represents an iteration.*/
void add_blast_iteration(xmlNode *node, void *iterations){
    if (!strcmp((char *)(node->name), "Iteration"))
        ds_vector_append((struct DSVector *)iterations, (void *)node);
}


/*Takes in the xmlNode for the root of a parsed BLAST XML tree and returns
  a vector of all of the nodes in the tree that represent iterations.*/
struct DSVector *get_blast_iterations(xmlNode *node){
    struct DSVector *iterations = ds_vector_create();
    traverse_blast_xml_r(node, add_blast_iteration, iterations);
    return iterations;
}

void write_fine_fasta(struct DSVector *oseqs){
    int i;
    FILE *temp = fopen("CaBLAST_fine.fasta", "w");
    if (!temp) {
        fprintf(stderr, "Could not open CaBLAST_fine.fasta for writing\n");
        return;
    }
    for (i = 0; i < oseqs->size; i++) {
        struct cb_seq *current_seq = (struct cb_seq *)ds_vector_get(oseqs, i);
        fprintf(temp, "> %s\n%s\n", current_seq->name, current_seq->residues);
    }
    fclose(temp);
}

/*Takes in the arguments for the program and returns a string of all of the
 *arguments for fine BLAST (the arguments after --blast-args) concatenated and
 *separated by a space.
 */
char *get_blast_args(struct opt_args *args){
    int i = 0;
    int index = -1;
    int length = 1;
    char *blast_args = NULL;
    for (i = 0; i < args->nargs; i++)
        if (index >= 0)
            length += (strlen(args->args[i]) + 1);
        else
            index = strcmp(args->args[i], "--blast-args") == 0 ? i : -1;
    blast_args = malloc(length*sizeof(*args));
    if (index == -1)
        *blast_args = '\0';
    else
        for (i = index + 1; i < args->nargs; i++) {
            blast_args = strcat(blast_args, args->args[i]);
            if (i < args->nargs - 1)
                blast_args = strcat(blast_args, " ");
        }
    return blast_args;
}

/*Takes in as input the vector of iterations from coarse BLAST, an index into
 *the vector (representing the index of which query we want to get the results
 *from), and the database we are using for search and returns a vector of
 *every original sequence section re-created from the calls to
 *cb_coarse_expand for the hits in the iteration we are expanding.
 */
struct DSVector *expand_blast_hits(struct DSVector *iterations, int index,
                                   struct cb_database *db){
    struct DSVector *expanded_hits = ds_vector_create();
    int i = 0, j = 0, k = 0;
    struct DSVector *hits = get_blast_hits((xmlNode *)
                                ds_vector_get(iterations, index));
    fprintf(stderr, "iteration: %d/%d\n", index+1, iterations->size);
    for (i = 0; i < hits->size; i++) {
        struct hit *current_hit = (struct hit *)ds_vector_get(hits, i);
        struct DSVector *hsps = current_hit->hsps;
        fprintf(stderr, "        hit: %d/%d\n", i+1, hits->size);
        for (j = 0; j < hsps->size; j++) {
            struct hsp *h = (struct hsp *)ds_vector_get(hsps, j);
            int16_t coarse_start = h->hit_from-1;
            int16_t coarse_end = h->hit_to-1;
            int32_t coarse_seq_id = current_hit->accession;
            struct DSVector *oseqs =
                cb_coarse_expand(db->coarse_db, db->com_db, coarse_seq_id,
                                  coarse_start, coarse_end, 50);
            fprintf(stderr, "                Hsp: %d/%d\n",
                                           j+1, hsps->size);
            for (k = 0; k < oseqs->size; k++)
                ds_vector_append(expanded_hits, ds_vector_get(oseqs, k));
            ds_vector_free_no_data(oseqs);
        }
    }
    return expanded_hits;
}


int
main(int argc, char **argv)
{
    int i = 0, j = 0;
    uint64_t dbsize = 0;
    struct cb_database *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;
    xmlDoc *doc = NULL;
    xmlNode *root = NULL;
    struct DSVector *iterations = NULL, *expanded_hits = NULL;

    conf = load_search_args();
    args = opt_config_parse(conf, argc, argv);

    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir fasta-file "
            "[ --blast_args BLASTN_ARGUMENTS ]\n",
            argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    db = cb_database_read(args->args[0], search_flags.map_seed_size);
    dbsize = read_int_from_file(8, db->coarse_db->file_params);
    blast_coarse(args, dbsize);

    doc = xmlReadFile("CaBLAST_temp_blast_results.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Could not parse CaBLAST_temp_blast_results.xml\n");
        return 0;
    }
    root = xmlDocGetRootElement(doc);
    iterations = get_blast_iterations(root);
    for (i = 0; i < iterations->size; i++) {
        expanded_hits = expand_blast_hits(iterations, i, db);
        ds_vector_free(expanded_hits);
    }

    /*Free the XML data and expanded hits*/
    for (i = 0; i < iterations->size; i++) {
        struct DSVector *iteration =
            (struct DSVector *)ds_vector_get(iterations, i);
        for (j = 0; j < iteration->size; j++) {
            struct hit *h = (struct hit *)ds_vector_get(iteration, j);
            ds_vector_free(h->hsps);
        }
    }
    cb_database_free(db);
    xmlFreeDoc(doc);
    return 0;
}
