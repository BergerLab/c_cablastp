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
    if (end <= start) {
        fprintf(stderr, "Invalid cb_range_node start (%d) must be greater than "
                        "end (%d)\n", start, end);
        return NULL;
    }
    node = malloc(sizeof(*node));
    node->range = malloc(sizeof(*(node->range)));
    node->range->start = start;
    node->range->end = end;
    node->sequence = malloc((end-start+2)*sizeof(*(node->sequence)));
    strncpy(node->sequence, sequence, (end-start+1));
    node->sequence[end-start+1] = '\0';
    node->left = NULL;
    node->right = NULL;
    return node;
}

/*Created a new empty cb_range_tree.*/
struct cb_range_tree *cb_range_tree_create(char *seq_name){
    struct cb_range_tree *tree = malloc(sizeof(*tree));
    tree->seq_name = malloc((strlen(seq_name)+1)*sizeof(*(tree->seq_name)));
    strncpy(tree->seq_name, seq_name, strlen(seq_name));
    tree->seq_name[strlen(seq_name)] = '\0';
    tree->root = NULL;
    return tree;
}

/*Frees a cb_range_node and all nodes in its subtrees.*/
void cb_range_node_free(struct cb_range_node *node){
    if (node != NULL) {
        free(node->range);
        free(node->sequence);
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
    free(node->sequence);
    node->sequence = sequence;
    node->range->start = start;
    node->range->end = end;
}

/*Takes in a tree, an expanded DNA sequence, and its starting and ending
 *indices and inserts a node created from those data into the tree and updates
 *the tree so that no ranges overlap.
 */
void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end){
    if (end <= start)
       fprintf(stderr, "Invalid range start (%d) must be less than end (%d)",
                       start, end);
    else
        tree->root = cb_range_node_insert(tree->root, sequence, start, end);
}

/*Find the last node in a node's left or right subtree that the range specified
  by start and end would overlap*/
struct cb_range_node *find_last_in_range(struct cb_range_node *current, int dir,
                                         int start, int end){
   if (end <= start) {
       fprintf(stderr, "Invalid range start (%d) must be less than end (%d)",
                       start, end);
       return NULL;
   }
   if (dir < 0)
       return (!current->left || current->left->range->end < start) ?
              current : find_last_in_range(current->left, dir, start, end);
   else
       return (!current->right || current->right->range->start > end) ?
              current : find_last_in_range(current->right, dir, start, end);
}

/*Called from cb_range_tree_insert to add an expanded DNA sequence and its
  indices to a cb_range_tree.*/
struct cb_range_node *cb_range_node_insert(struct cb_range_node *current,
                                           char *sequence, int start, int end){
    struct cb_range_node *temp = NULL, *parent = NULL;

    /*If the node we are at is a null pointer, then create a new node with the
      sequence and range passed in and return the node to add it to the tree.*/
    if (!current)
        return cb_range_node_create(sequence, start, end);

    /*If the new range overlaps the current node's range and goes past the
     *current node's range, update the current node's range and delete any
     *nodes in the current node's subtrees that the current node's new range
     *overlaps.
     */
    if (start <= current->range->end && end >= current->range->start) {
        /*If the start goes past the current node's start, update the current
         *node's start and delete any nodes in the left subtree that would
         *overlap the new range, setting the current node's start to the start
         *of the last node the current node overlaps and updating the node's
         *sequence to include any sequences the range overlaps to the left.
         */
        if (start < current->range->start) {
            char *merged = cb_range_merge(sequence, start, end,
                                          current->sequence,
                                          current->range->start,
                                          current->range->end);
            free(current->sequence);
            current->sequence = merged;
            current->range->start = start;

            temp = current->left;
            parent = find_last_in_range(current, -1, start, end);

            if (current != parent) {
fprintf(stderr, "<=start = %d...     ", start);
                if (parent->range->start < start) {
                    char *merged =
                        cb_range_merge(parent->sequence, parent->range->start,
                                       parent->range->end,
                                       current->sequence, current->range->start,
                                       current->range->end);
                    free(current->sequence);
                    current->sequence = merged;
fprintf(stderr, "                                 %s\n", merged);
                    current->range->start = parent->range->start;
                }
                current->left = parent->left;
                parent->left = NULL;
fprintf(stderr, "Now the current node's range is %d-%d\n", current->range->start, current->range->end);
                cb_range_node_free(temp);
            }
        }
        /*If the end goes past the current node's end, update the current node's
         *end and delete any nodes in the right subtree that would overlap the
         *new range, setting the current node's end to the end of the last node
         *the current node overlaps and updating the node's sequence to include
         *any sequences the range overlaps to the right.
         */
        if (end > current->range->end) {
            temp = current->right;
            parent = find_last_in_range(current, 1, start, end);
            
            if (current != parent) {
fprintf(stderr, "=>");
                char *merged =
                    cb_range_merge(current->sequence, current->range->start,
                                   current->range->end,
                                   parent->sequence, parent->range->start,
                                   parent->range->end);
                free(current->sequence);
                current->sequence = merged;
                current->range->start = parent->range->start;
                current->right = parent->right;
fprintf(stderr, "Now the current node's range is %d-%d\n", current->range->start, current->range->end);
                parent->right = NULL;
                cb_range_node_free(temp);
            }
        }
    }
    /*If there is no overlap between the current node's range and the range
      passed in, call cb_range_node_insert on whichever subtree the range would
      be in to make any necessary updates to that subtree.*/
    else
        if (end < current->range->start)
            current->left = cb_range_node_insert(current->left, sequence,
                                                 start, end);
        else
            current->right = cb_range_node_insert(current->right, sequence,
                                                  start, end);
    return current;
}

/*Takes in the sequences and starting and ending indices of two sections of an
  original sequence and merges the sequences together.*/
char *cb_range_merge(char *left_seq, int left_start, int left_end,
                     char *right_seq, int right_start, int right_end){
fprintf(stderr, "Merging %d-%d (%s) and %d-%d (%s) to ", left_start, left_end, left_seq, right_start, right_end, right_seq);
    char *merged = malloc((right_end-left_start)*sizeof(*merged));
    int index = 0, i = 0;
    for (i = 0; i < right_start - left_start; i++)
        merged[index++] = left_seq[i];
    for (i = 0; i < right_start - left_start; i++)
        fprintf(stderr, "%d\n", left_seq[i]);
    for (i = 0; i < right_end - right_start; i++)
        merged[index++] = right_seq[i];
    merged[right_end-left_start] = '\0';
fprintf(stderr, "%s      index = %d\n", merged, index);
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
    /*fprintf(stderr, "%s\n%s\n", tree->seq_name, node->sequence);*/
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
    fprintf((FILE *)out, "> %s\n%s\n", tree->seq_name, node->sequence);
}

/*Uses cb_range_tree_traverse and cb_range_node_output to output the sequences
  in the range tree passed in to the specified file pointer passed in.*/
void cb_range_tree_output(struct cb_range_tree *tree, FILE *out){
    cb_range_tree_traverse(tree, cb_range_node_output, (void *)out);
}


void cb_range_node_output_range(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out){
    fprintf((FILE *)out, "> %d-%d... %s\n", node->range->start,
                                            node->range->end, node->sequence);
}

void cb_range_tree_output_range(struct cb_range_tree *tree, FILE *out){
    fprintf(stderr, "Outputting range tree\n");
    cb_range_tree_traverse(tree, cb_range_node_output_range, (void *)out);
    fprintf(stderr, "Finished outputting range tree\n");
}

