#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "range_tree.h"

/*Creates a new cb_range_node with the range passed in start must be greater
 *than end.  The sequence passed into "sequence" must not be a string literal
 *so it can be freed.
 */
struct cb_range_node *cb_range_node_create(char *sequence, int start, int end){
fprintf(stderr, "    Started cb_node_create\n");
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
    node->sequence = malloc((end-start+1)*sizeof(*(node->sequence)));
    strcpy(node->sequence, sequence);
    node->left = NULL;
    node->right = NULL;
fprintf(stderr, "    Completed cb_node_create\n");
    return node;
}

/*Created a new empty cb_range_tree.*/
struct cb_range_tree *cb_range_tree_create(){
    struct cb_range_tree *tree = malloc(sizeof(*tree));
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

/*Takes in a tree, an expanded DNA sequence, and its starting and ending
 *indices and inserts a node created from those data into the tree and updates
 *the tree so that no ranges overlap.
 */
void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end){
fprintf(stderr, "----------Starting cb_range_tree_insert\n");
    if (end <= start)
       fprintf(stderr, "Invalid range start (%d) must be less than end (%d)",
                       start, end);
    else
        tree->root = cb_range_node_insert(tree->root, sequence, start, end);
fprintf(stderr, "----------Finished cb_range_tree_insert\n");
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
       return (!current->left || current->left->range->end > start) ?
              current : find_last_in_range(current->left, dir, start, end);
   else
       return (!current->right || current->right->range->start < end) ?
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
    if (start < current->range->end && end > current->range->start) {
fprintf(stderr, "Starting overlap\n");
        /*If the start goes past the current node's start, update the current
         *node's start and delete any nodes in the left subtree that would
         *overlap the new range, setting the current node's start to the start
         *of the last node the current node overlaps and updating the node's
         *sequence to include any sequences the range overlaps to the left.
         */
        if (start < current->range->start) {
            current->range->start = start;
            temp = current->left;
fprintf(stderr, "Called find_last_in_range\n");
            parent = find_last_in_range(current, -1, start, end);
fprintf(stderr, "Finished find_last_in_range\n");

            if (current != parent) {
                char *merged =
                    cb_range_merge(current->sequence, current->range->start,
                                   current->range->end,
                                   parent->sequence, parent->range->start,
                                   parent->range->end);
                free(current->sequence);
                current->sequence = merged;
                current->range->start = parent->range->start;
                current->left = parent->left;
                parent->left = NULL;
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
            current->range->end = end;
            temp = current->right;
fprintf(stderr, "Called find_last_in_range\n");
            parent = find_last_in_range(current, 1, start, end);
fprintf(stderr, "Finished find_last_in_range\n");
            
            if (current != parent) {
                char *merged =
                    cb_range_merge(parent->sequence, parent->range->start,
                                   parent->range->end,
                                   current->sequence, current->range->start,
                                   current->range->end);
                free(current->sequence);
                current->sequence = merged;
                current->range->start = parent->range->start;
                current->right = parent->right;
                parent->right = NULL;
                cb_range_node_free(temp);
            }
        }
fprintf(stderr, "Finished overlap\n");
    }
    /*If there is no overlap between the current node's range and the range
      passed in, call cb_range_node_insert on whichever subtree the range would
      be in to make any necessary updates to that subtree.*/
    else{fprintf(stderr, "Entered else branch\n");
        if (end < current->range->start)
            current->left = cb_range_node_insert(current->left, sequence,
                                                 start, end);
        else
            current->right = cb_range_node_insert(current->right, sequence,
                                                  start, end);
fprintf(stderr, "Finished else branch\n");}
    return current;
}

/*Takes in the sequences and starting and ending indices of two sections of an
  original sequence and merges the sequences together.*/
char *cb_range_merge(char *left_seq, int left_start, int left_end,
                     char *right_seq, int right_start, int right_end){
fprintf(stderr, "cb_range_merge\n%s\n%s\n\n", left_seq, right_seq);
    char *merged = malloc((right_end-left_start)*sizeof(*merged));
    int index = 0, i = 0;
    for (i = 0; i < right_start-left_start; i++)
        merged[index++] = left_seq[i];
    for (i = 0; i < right_end-right_start; i++)
        merged[index++] = right_seq[i];
fprintf(stderr, "Finished cb_range_merge\n");
    return merged;
}