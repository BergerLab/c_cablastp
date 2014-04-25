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
#include "range_tree.h"
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
           "%s -out CaBLAST_temp_blast_results.xml";
    int command_length = strlen(blastn_command) + strlen(input_path) +
                                                  strlen(args->args[1]) + 31;
    char *blastn = malloc(command_length * sizeof(*blastn));
    assert(blastn);

    sprintf(blastn,"blastn -db %s -outfmt 5 -query %s -dbsize %lu -task blastn"
                   " -evalue %s -out CaBLAST_temp_blast_results.xml",
           input_path, args->args[1], dbsize, search_flags.coarse_evalue);
    if (!search_flags.hide_messages)
        fprintf(stderr, "Running coarse BLAST:\n%s\n\n", blastn);
    system(blastn); /*Run coarse BLAST*/
    free(blastn);
}

/*Runs BLAST on the fine FASTA file*/
void blast_fine(char *query, uint64_t dbsize,
                char *blast_args, bool has_evalue){
    char *default_evalue = "-evalue 1e-30";
    int evalue_len = has_evalue ? 0 : strlen(default_evalue);

    char *blastn_command =
           "blastn -db CaBLAST_fine.fasta -query  "
           "-dbsize  -task blastn ";
    int command_length = strlen(blastn_command) + strlen(query) + 31 +
                         evalue_len + strlen(blast_args);
    char *blastn = malloc(command_length * sizeof(*blastn));
    assert(blastn);

    sprintf(blastn,
            "blastn -db CaBLAST_fine.fasta -query %s "
            "-dbsize %lu -task blastn %s %s",
            query, dbsize, (has_evalue ? "" : default_evalue), blast_args);
    if (!search_flags.hide_messages)
        fprintf(stderr, "Running fine BLAST:\n%s\n\n", blastn);
    system(blastn); /*Run fine BLAST*/
    free(blastn);

    /*Delete the FASTA file full of hit expansions and its database unless the
      no-cleanup flag is activated.*/
    if (!search_flags.no_cleanup)
        system("rm CaBLAST_fine.fasta CaBLAST_fine.fasta.*");
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
    struct hsp *h;
    bool hit_frame_fwd = true;

    h = malloc(sizeof(*h));
    assert(h);

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
    assert(h);

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

/*Takes in a vector of expansion structs for the expanded BLAST hits from a
  query and outputs them to the FASTA file CaBLAST_fine.fasta.*/
void write_fine_db(struct DSVector *oseqs){
    int i;
    FILE *fine_fasta = fopen("CaBLAST_fine.fasta", "a");
    if (!fine_fasta) {
        if (!search_flags.hide_messages)
            fprintf(stderr,
                    "Could not open CaBLAST_fine.fasta for appending\n");
        return;
    }
    for (i = 0; i < oseqs->size; i++) {
        struct cb_seq *current_seq =
            ((struct cb_hit_expansion *)ds_vector_get(oseqs, i))->seq;
        fprintf(fine_fasta, "> %s\n%s\n", current_seq->name,
                                          current_seq->residues);
    }
    fclose(fine_fasta);
    system("makeblastdb -dbtype nucl -in CaBLAST_fine.fasta "
                                    "-out CaBLAST_fine.fasta");
}

void write_fine_db_from_trees(struct DSHashMap *range_trees){
    FILE *tree_expansions = fopen("tree.fasta", "w");
    int i = 0;

    ds_hashmap_sort_keys(range_trees);
    for (i = 0; i < range_trees->keys->size; i++){
        int32_t key =
            ((struct DSHashKey *)ds_vector_get(range_trees->keys, i))->key.i;
        struct cb_range_tree *tree =
            (struct cb_range_tree *)ds_geti(range_trees, key);
        cb_range_tree_output_fasta(tree, tree_expansions);
    }
    fclose(tree_expansions);
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
            index = (strcmp(args->args[i],
                            "--blast-args") == 0 && index == -1) ? i : index;
    blast_args = malloc(length*sizeof(*args));
    assert(blast_args);
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

/*Checks the string of additional BLAST arguments to CaBLAST and returns true
 *if the arguments includes the specified argument or false if the additional
 *BLAST arguments does not include the argument or if no additional BLAST
 *arguments were given.
 */
bool has_blast_arg(struct opt_args *args, char *arg){
    int i = 0;
    int index = -1;
    for (i = 0; i < args->nargs; i++)
        index = (strcmp(args->args[i], 
                        "--blast-args") == 0 && index == -1) ? i : index;
    if (index == -1)
        return false;
    else
        for (i = index + 1; i < args->nargs; i++)
            if (strcmp(args->args[i], arg) == 0)
                return true;
    return false;
}

/*Takes in as input the vector of iterations from coarse BLAST, an index into
 *the vector (representing the index of which query we want to get the results
 *from), and the database we are using for search and returns a vector of
 *every original sequence section re-created from the calls to
 *cb_coarse_expand for the hits in the iteration we are expanding.
 */
struct DSVector *expand_blast_hits(struct DSVector *iterations, int index,
                                   struct cb_database *db,
                                   struct DSHashMap *range_trees){
    struct DSVector *expanded_hits = ds_vector_create();
    int i = 0, j = 0, k = 0;
    struct DSVector *hits = get_blast_hits((xmlNode *)
                                ds_vector_get(iterations, index));
    for (i = 0; i < hits->size; i++) {
        struct hit *current_hit = (struct hit *)ds_vector_get(hits, i);
        struct DSVector *hsps = current_hit->hsps;
        for (j = 0; j < hsps->size; j++) {
            struct hsp *h = (struct hsp *)ds_vector_get(hsps, j);
            int32_t coarse_start = h->hit_from-1;
            int32_t coarse_end = h->hit_to-1;
            int32_t coarse_seq_id = current_hit->accession;
            struct DSVector *oseqs =
                cb_coarse_expand(db->coarse_db, db->com_db,
                                 range_trees, coarse_seq_id,
                                 coarse_start, coarse_end, 50);
            for (k = 0; k < oseqs->size; k++)
                ds_vector_append(expanded_hits, ds_vector_get(oseqs, k));
            ds_vector_free_no_data(oseqs);
        }
        ds_vector_free(current_hit->hsps);
        free(current_hit);
    }
    ds_vector_free_no_data(hits);
    return expanded_hits;
}

char *progress_bar(int current, int iterations){
    char *bar;
    int bars = (int)(((float)(current+1)/iterations)*50);
    int b = 0;

    bar = malloc(53*sizeof(*bar));
    assert(bar);

    bar[0] = '[';
    for(b = 0; b < 50; b++)
        bar[b+1] = b < bars ? '|' : ' ';
    bar[51] = ']';
    bar[52] = '\0';
    return bar;
}

int32_t range_start_cmp(void *a, void *b){
    return ((struct cb_hit_expansion *)a)->offset -
           ((struct cb_hit_expansion *)b)->offset;
}

int
main(int argc, char **argv)
{
    int i = 0, j = 0;
    uint64_t dbsize = 0, expanded_dbsize = 0;
    struct cb_database *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;
    xmlDoc *doc = NULL;
    xmlNode *root = NULL;
    struct DSVector *iterations = NULL, *expanded_hits = NULL,
                    *queries = NULL, *oseqs = ds_vector_create();
    struct DSHashMap *range_trees = ds_hashmap_create();
    struct fasta_seq *query = NULL;
    FILE *query_file = NULL;
    char *blast_args = NULL;
    bool has_evalue = false;

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
    blast_args = get_blast_args(args);
    has_evalue = has_blast_arg(args, "-evalue");

    db = cb_database_read(args->args[0], search_flags.map_seed_size);
    dbsize = read_int_from_file(8, db->coarse_db->file_params);
    blast_coarse(args, dbsize);
    query_file = fopen(args->args[1], "r");

    queries = ds_vector_create();
    query = fasta_read_next(query_file, "");
    while (query) {
        ds_vector_append(queries, (void *)query);
        query = fasta_read_next(query_file, "");
    }

    fclose(query_file);

    system("rm CaBLAST_fine.fasta");

    /*Parse the XML file generated from coarse BLAST and get its iterations.*/
    doc = xmlReadFile("CaBLAST_temp_blast_results.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Could not parse CaBLAST_temp_blast_results.xml\n");
        return 0;
    }
    root = xmlDocGetRootElement(doc);
    iterations = get_blast_iterations(root);

    if (!search_flags.hide_messages)
        fprintf(stderr, "Expanding coarse BLAST hits\n");

    /*Expand any BLAST hits we got from the current query sequence during
      coarse BLAST.*/
    for (i = 0; i < iterations->size; i++) {
        /*If the show-progress flag is activated, display the progress bar*/
        if (!search_flags.hide_messages) {
            char *bar = progress_bar(i, iterations->size);
            fprintf(stderr, "\r");
            fprintf(stderr, "iteration: %d/%d", i+1, iterations->size);
            fprintf(stderr, " %s", bar);
            free(bar);
        }

        expanded_hits = expand_blast_hits(iterations, i, db, range_trees);
        query = (struct fasta_seq *)ds_vector_get(queries, i);
        for (j = 0; j < expanded_hits->size; j++) {
            struct cb_hit_expansion *current_expansion =
                (struct cb_hit_expansion *)ds_vector_get(expanded_hits, j);
            expanded_dbsize += current_expansion->seq->length;
            ds_vector_append(oseqs, (void *)current_expansion);
        }
        ds_vector_free_no_data(expanded_hits);
    }

    /*if (!search_flags.hide_messages)
        fprintf(stderr, "\n\nWriting database for fine BLAST\n\n");
    write_fine_db(oseqs);*/

    ds_vector_sort(oseqs, range_start_cmp);

    for (i = 0; i < oseqs->size; i++) {
        struct cb_hit_expansion *current =
            (struct cb_hit_expansion *)ds_vector_get(oseqs, i);
        int current_start = current->offset;
        int current_end = current_start + current->seq->length;
        int org_seq_id = current->seq->id;
        char *current_seq = current->seq->residues;
        struct cb_range_tree *tree = (struct cb_range_tree *)
                                         ds_geti(range_trees, org_seq_id);
        struct cb_range_node *current_node =
            cb_range_tree_find(tree, current_start, current_end);
        if ((!is_complete_overlap(current_node->seq, current_seq))){
            fprintf(stderr, "Error in tree from sequence #%d, %d-%d\n"
                    "%s\ndoes not overlap\n%d-%d\n%s\n\n", org_seq_id,
                    current_start, current_end, current_seq,
                    current_node->start, current_node->end, current_node->seq);
            printf("Error in tree from sequence #%d, %d-%d\n"
                    "%s\ndoes not overlap\n%d-%d\n%s\n\n", org_seq_id,
                    current_start, current_end, current_seq,
                    current_node->start, current_node->end, current_node->seq);
        }
    }

printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

    for (i = 0; i < 10; i++) {
        printf("\n\nPrinting sequences from original sequence #%d\n\n", i);
        for (j = 0; j < oseqs->size; j++) {
            struct cb_hit_expansion *current =
                (struct cb_hit_expansion *)ds_vector_get(oseqs, j);
            int current_start = current->offset;
            int current_end = current_start + current->seq->length;
            int org_seq_id = current->seq->id;
            if (org_seq_id == i)
                printf("%d-%d\n", current_start, current_end);
        }
        printf("___________________________________________________\n");
    }

printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");

printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    for (i = 0; i < 10; i++) {
        printf("\n\nPrinting tree #%d\n\n", i);
        struct cb_range_tree *tree = (struct cb_range_tree *)
                                      ds_geti(range_trees, i);
        cb_range_tree_output_ranges(tree, stdout);
    }
printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");


/*  write_fine_db_from_trees(range_trees);
    blast_fine(args->args[1], expanded_dbsize, blast_args, has_evalue);*/

    ds_vector_free_no_data(queries);

    /*Free the XML data and expanded hits*/
    for (i = 0; i < iterations->size; i++) {
        struct DSVector *iteration =
            (struct DSVector *)ds_vector_get(iterations, i);
        for (j = 0; j < iteration->size; j++) {
            struct hit *h = (struct hit *)ds_vector_get(iteration, j);
            ds_vector_free(h->hsps);
            free(h);
        }
    }

    for (i = 0; i < oseqs->size; i++)
        cb_hit_expansion_free(
            (struct cb_hit_expansion *)ds_vector_get(oseqs, i));

    ds_vector_free_no_data(oseqs);
/*    ds_hashmap_free(range_trees, true, true);*/

    cb_database_free(db);
    xmlFreeDoc(doc);

    /*Free the coarse BLAST results file if the --no-cleanup flag is not being
      used.*/
    if(!search_flags.no_cleanup)
        system("rm CaBLAST_temp_blast_results.xml");
    return 0;
}
