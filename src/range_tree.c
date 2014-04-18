#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "range_tree.h"

/*Creates a new cb_range_node with the range passed in start must be greater
 *than end.  The sequence passed into "sequence" must not be a string literal
 *so it can be freed.
 */
struct cb_range_node *cb_range_node_create(char *sequence, int start, int end){
    struct cb_range_node *node = NULL;
    assert(end > start);

    node = malloc(sizeof(*node));
    node->start = start;
    node->end = end;

    node->seq = malloc((end-start+2)*sizeof(*(node->seq)));
    strncpy(node->seq, sequence, (end-start+1));
    node->seq[end-start+1] = '\0';

    node->left = NULL;
    node->right = NULL;
    return node;
}

/*Creates a new empty cb_range_tree.*/
struct cb_range_tree *cb_range_tree_create(char *seq_name){
    struct cb_range_tree *tree = malloc(sizeof(*tree));
    if (seq_name != NULL) {
        tree->seq_name = malloc((strlen(seq_name)+1)*sizeof(*(tree->seq_name)));
        strncpy(tree->seq_name, seq_name, strlen(seq_name));
        tree->seq_name[strlen(seq_name)] = '\0';
    }
    else
        tree->seq_name = NULL;
    tree->root = NULL;
    return tree;
}

/*Frees a cb_range_node and all nodes in its subtrees.*/
void cb_range_node_free(struct cb_range_node *node){
    if (node != NULL) {
        if (node->seq != NULL)
            free(node->seq);
        cb_range_node_free(node->left);
        cb_range_node_free(node->right);
        free(node);
    }
}

/*Frees a cb_range_tree and all of its nodes*/
void cb_range_tree_free(struct cb_range_tree *tree){
    if (tree != NULL) {
        cb_range_node_free(tree->root);
        free(tree);
    }
}

/*Takes in a range tree node, a new DNA sequence, and starting and the new
 *sequence's starting and ending indices and sets the node's sequence and
 *indices to the sequence and indices passed in.
 */
void cb_range_node_update(struct cb_range_node *node, char *sequence,
                          int start, int end){
    int len = -1;

    assert(end > start);
    assert(sequence != NULL);


    len = strlen(sequence);
    if (node->seq != NULL)
        free(node->seq);
    node->seq = malloc((len+1)*sizeof(*(node->seq)));
    strncpy(node->seq, sequence, len);
    node->seq[len] = '\0';

    node->start = start;
    node->end = end;
}

/*Takes in a tree, an expanded DNA sequence, and its starting and ending
 *indices and inserts a node created from those data into the tree and updates
 *the tree so that no ranges overlap.
 */
void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end){
    assert(end > start);
    assert(sequence != NULL);
    tree->root = cb_range_node_insert(tree->root, sequence, start, end);
}

/*Create a range node data struct*/
struct cb_range_node_data *cb_range_node_data_create(int start, int end,
                                                     char *sequence){
    struct cb_range_node_data *d = malloc(sizeof(*d));
    int len = -1;

    assert(end > start);
    assert(sequence != NULL);


    d->start = start;
    d->end = end;

    len = strlen(sequence);
    d->seq = malloc((len+1)*sizeof(*(d->seq)));
    strncpy(d->seq, sequence, len);
    d->seq[len] = '\0';
    return d;
}

/*Update a range node data struct*/
void cb_range_node_data_update(struct cb_range_node_data *d,
                               int start, int end, char *sequence){
    int len = -1;

    assert(end > start);
    assert(sequence != NULL);

    d->start = start;
    d->end = end;

    len = strlen(sequence);
    if (d->seq != NULL)
        free(d->seq);
    d->seq = malloc((len+1)*sizeof(*(d->seq)));
    strncpy(d->seq, sequence, len);
    d->seq[len] = '\0';   
}

/*Free a range node data struct*/
void cb_range_node_data_free(struct cb_range_node_data *d){
    assert(d != NULL);
    if (d->seq != NULL)
        free(d->seq);
    free(d);
}

/*Removes nodes in the tree that overlap with range being passed in with
 *"start" and "end" in the direction specified by "dir" (-1 = left, 1 = right)
 *and returns the data of the node overlapping the range that is farthest in
 *the specified direction.
 */
struct cb_range_node_data *remove_overlap(struct cb_range_node *cur,
                                          int start, int end, int dir,
                                          struct cb_range_node_data *data){
    struct cb_range_node *next, *new_next;
    bool left, overlap, next_overlap;
    
    assert(end > start); /*Make sure we are using a valid range*/
    assert(cur!= NULL); /*remove_overlap is never called from a NULL pointer.*/

    left = dir < 0;
    overlap = start < cur->end && end > cur->start;
    next = overlap ? (left ? cur->left : cur->right) :
                     (left ? cur->right : cur->left);
    /*If there is no next node, we already found the farthest node in the
      specified direction to overlap the range, so return its data.*/
    if (!next)
        return data;
    next_overlap = start < next->end && end > next->start;

    /*If the next node overlaps the range, the next node needs to be deleted,
     *so see if the data we are returning needs to be updated and then free the
     *the next node and its subtree in the direction opposite of dir, making it
     *so the next node is now the next node after next.  Since we don't know if
     *new_next is NULL or overlaps the range, call the next call to
     *remove_overlap from the current node.
     */
    if (next_overlap) {
        if (left && next->start < data->start)
            cb_range_node_data_update(data, next->start, next->end, next->seq);
        new_next = left ? next->left : next->right;
        if (left) {
            cur->left = new_next;
            next->left = NULL;
        }
        else {
            cur->right = new_next;
            next->right = NULL;
        }
        /*Free the next node and the subtree new_next is not in.*/
        cb_range_node_free(next);

        return remove_overlap(cur, start, end, dir, data);
    }
    /*If the next node does not overlap the range, call the next call to
      remove_overlap from the next node.*/
    else
        return remove_overlap(next, start, end, dir, data);
}

struct cb_range_node *cb_range_node_insert(struct cb_range_node *node,
                                           char *sequence, int start, int end){
    bool overlap;
    struct cb_range_node_data *last;
    if (!node)
        return cb_range_node_create(sequence, start, end);
    overlap = start < node->end && end > node->start;
    if (!overlap) {
        if (end < node->start)
            node->left=cb_range_node_insert(node->left, sequence, start, end);
        else
            node->right=cb_range_node_insert(node->right, sequence, start, end);
    }
    else {
        if (start < node->start) {
            cb_range_node_update(node,
                                 cb_range_merge(sequence, start, end, node->seq,
                                                node->start, node->end),
                                 start, node->end);
            last = cb_range_node_data_create(start, end, sequence);
            remove_overlap(node, start, end, -1, last);
        }
        if (start < node->end) {
            ;
        }
    }
    return node;
}

/*Takes in a range tree node, starting and ending indices for the range being
 *added to the tree, the direction that the range overlaps the node this
 *function was originally called from, and a pointer to a range tree node and
 *finds the node in the tree that is farthest in that direction to overlap the
 *range being added.
 */
/*struct cb_range_node *find_last_overlap(struct cb_range_node *node,
                                        int start, int end, int dir,
                                        struct cb_range_node *last){
    assert(end > start);
    /*If we are at the end of the tree, return the last node we found that
      overlaps the range.*//*
    if (!node) 
        return last;
    else {
        /*If the node overlaps the range passed in, go in the direction
         *specified by dir to look for nodes further in that direction that
         *overlap the range.  Otherwise, we are past the range so look in the
         *opposite direction.
         *//*
        bool overlap = start < node->end && end > node->start;
        struct cb_range_node *next = overlap ? (dir?node->left:node->right) :
                                               (dir?node->right:node->left);
        return find_last_overlap(node, start, end, dir, (overlap?node:last));
    }
}*/

/*Takes in the sequences and starting and ending indices of two sections of an
  original sequence and merges the sequences together.*/
char *cb_range_merge(char *left_seq, int left_start, int left_end,
                     char *right_seq, int right_start, int right_end){
fprintf(stderr, "Merging %d-%d and %d-%d\n", left_start, left_end, right_start, right_end);
    int index = 0, i = 0;
    char *merged = NULL;

    assert(left_end > left_start && right_end > right_start);
    assert(left_end > right_start && right_end > left_start);
    assert(left_seq != NULL && right_seq != NULL);

    merged = malloc((right_end-left_start)*sizeof(*merged));
    for (i = 0; i < right_start - left_start; i++)
        merged[index++] = left_seq[i];
    for (i = 0; i < right_end - right_start; i++)
        merged[index++] = right_seq[i];
    merged[right_end-left_start] = '\0';

    return merged;
}

/*A pre-order traversal function for range trees.  Takes in a range tree node,
 *its tree, a function that takes in a range tree node, its tree, and a void *,
 *and an accumulator and calls the function for the node and each node in its
 *subtrees.
 */
void cb_range_node_traverse(struct cb_range_node *node,
                            struct cb_range_tree *tree,
                            void (*f)(struct cb_range_node *,
                                      struct cb_range_tree *,void *),
                            void *acc){
    if (node->left != NULL)
        cb_range_node_traverse(node->left, tree, f, acc);
    f(node, tree, acc);
    if (node->right != NULL)
        cb_range_node_traverse(node->right, tree, f, acc);
}

/*Takes in a range tree, a function taking in a range tree node, its tree,
 *and a void *, and an accumulator and carries out a pre-order traversal on
 *the tree.
 */
void cb_range_tree_traverse(struct cb_range_tree *tree,
                            void (*f)(struct cb_range_node *,
                                      struct cb_range_tree *,void *),
                            void *acc){
    if (tree != NULL && tree->root != NULL)
        cb_range_node_traverse(tree->root, tree, f, acc);
}

/*A function to be passed into cb_range_tree_traverse that takes in a range
 *tree node, its tree, and a FILE pointer and outputs the node's sequence to
 *the file.
 */
void cb_range_node_output(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out){
    fprintf((FILE *)out, "> %s\n%s\n", tree->seq_name, node->seq);
}

/*Uses cb_range_tree_traverse and cb_range_node_output to output the sequences
  in the range tree passed in to the file pointer passed in.*/
void cb_range_tree_output(struct cb_range_tree *tree, FILE *out){
    cb_range_tree_traverse(tree, cb_range_node_output, (void *)out);
}

/*A function to be passed into cb_range_tree_traverse that takes in a range
 *tree node, its tree, and a FILE pointer and outputs the node's range to
 *and its sequence to the specified file pointer.
 */
void cb_range_node_output_data(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out){
    fprintf((FILE *)out, "> %d-%d... %s\n", node->start, node->end, node->seq);
}

/*Uses cb_range_tree_traverse and cb_range_node_output to output the sequences
  and ranges in the range tree passed in to the file pointer passed in.*/
void cb_range_tree_output_data(struct cb_range_tree *tree, FILE *out){
    fprintf(stderr, "Outputting range tree\n");
    cb_range_tree_traverse(tree, cb_range_node_output_data, (void *)out);
    fprintf(stderr, "Finished outputting range tree\n");
}

