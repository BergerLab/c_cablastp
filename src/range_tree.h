#include <stdio.h>
#include <stdlib.h>

#ifndef __CABLAST_RANGE_TREE_H__
#define __CABLAST_RANGE_TREE_H__

struct cb_range_node {
    char *seq;
    int start;
    int end;
    struct cb_range_node *left;
    struct cb_range_node *right;
};

struct cb_range_tree {
    char *seq_name;
    struct cb_range_node *root;
};

/*A struct for getting the data of a range tree node without its pointers.*/
struct cb_range_node_data{
    int start;
    int end;
    char *seq;
};

struct cb_range_node *cb_range_node_create(char *sequence, int start, int end);
struct cb_range_tree *cb_range_tree_create(char *name);
struct cb_range_node_data *cb_range_node_data_create(int start, int end,
                                                     char *sequence);

void cb_range_node_free(struct cb_range_node *node);
void cb_range_tree_free(struct cb_range_tree *tree);
void cb_range_node_data_free(struct cb_range_node_data *d);

void cb_range_node_update(struct cb_range_node *node, char *sequence,
                          int start, int end);
void cb_range_node_data_update(struct cb_range_node_data *d,
                               int start, int end, char *sequence);

void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end);
struct cb_range_node_data *remove_overlap(struct cb_range_node *cur,
                                          int start, int end, int dir,
                                          struct cb_range_node_data *data);
struct cb_range_node *cb_range_node_insert(struct cb_range_node *current,
                                           char *sequence, int start, int end);
char *cb_range_merge(char *left_seq, int left_start, int left_end,
                     char *right_seq, int right_start, int right_end);


void cb_range_node_traverse(struct cb_range_node *node,
                            struct cb_range_tree *tree,
                            void (*f)(struct cb_range_node *,
                                      struct cb_range_tree *,void *),
                            void *acc);
void cb_range_tree_traverse(struct cb_range_tree *tree,
                            void (*f)(struct cb_range_node *,
                                      struct cb_range_tree *,void *),
                            void *acc);
void cb_range_node_output(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out);
void cb_range_tree_output(struct cb_range_tree *tree, FILE *out);

void cb_range_node_output_data(struct cb_range_node *node,
                               struct cb_range_tree *tree, void *out);

void cb_range_tree_output_data(struct cb_range_tree *tree, FILE *out);

#endif
