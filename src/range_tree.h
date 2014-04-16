#include <stdio.h>
#include <stdlib.h>

#ifndef __CABLAST_RANGE_TREE_H__
#define __CABLAST_RANGE_TREE_H__

struct cb_range{
    int start;
    int end;
};

struct cb_range_node {
    struct cb_range *range;
    char *sequence;
    struct cb_range_node *left;
    struct cb_range_node *right;
};

struct cb_range_tree {
    char *seq_name;
    struct cb_range_node *root;
};

struct cb_range_node *cb_range_node_create(char *sequence, int start, int end);
struct cb_range_tree *cb_range_tree_create(char *name);
void cb_range_node_free(struct cb_range_node *node);
void cb_range_tree_free(struct cb_range_tree *tree);
void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end);
struct cb_range_node *find_last_in_range(struct cb_range_node *current, int dir,
                                         int start, int end);
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

void cb_range_node_output_range(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out);

void cb_range_tree_output_range(struct cb_range_tree *tree, FILE *out);

#endif
