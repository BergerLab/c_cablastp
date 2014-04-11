#include <stdio.h>
#include <stdlib.h>
#include "range_tree.h"

/*Creates a new cb_range_node with the range passed in start must be greater
  than end.*/
struct cb_range_node *cb_range_node_create(int start, int end){
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
    node->left = NULL;
    node->right = NULL;
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

/*Insert a node with the specified start and end into a cb_range_tree struct
  and update the tree so no ranges overlap.*/
void cb_range_tree_insert(struct cb_range_tree *tree, int start, int end){
    tree->root = cb_range_node_insert(tree->root, start, end);
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

/*Called from cb_range_tree_insert to add the specified range to a
  cb_range_tree.*/
struct cb_range_node *cb_range_node_insert(struct cb_range_node *current,
                                           int start, int end){
    struct cb_range_node *temp = NULL, *parent = NULL;

    /*If the node we are at is a null pointer, then create a new node with the
      range passed in and return it to add it to the tree.*/
    if (!current)
        return cb_range_node_create(start, end);

    /*If the new range overlaps the current node's range and goes past the
     *current node's range, update the current node's range and delete any
     *nodes in the current node's subtrees that the current node's new range
     *overlaps.
     */
    if (start < current->range->end && end > current->range->start) {
        /*If the start goes past the current node's start, update the current
          node's start and delete any nodes in the left subtree that would
          overlap the new range, setting the current node's start to the start
          of the last node the current node overlaps*/
        if (start < current->range->start) {
            current->range->start = start;
            temp = current->left;
            parent = find_last_in_range(current, -1, start, end);

            if (current != parent) {
                current->range->start = parent->range->start;
                current->left = parent->left;
                parent->left = NULL;
                cb_range_node_free(temp);
            }
        }
        /*If the end goes past the current node's end, update the current node's
          end and delete any nodes in the right subtree that would overlap the
          new range, setting the current node's end to the end of the last node
          the current node overlaps*/
        if (end > current->range->end) {
            current->range->end = end;
            temp = current->right;
            parent = find_last_in_range(current, 1, start, end);
            
            if (current != parent) {
                current->range->start = parent->range->start;
                current->right = parent->right;
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
            current->left = cb_range_node_insert(current->left, start, end);
        else
            current->right = cb_range_node_insert(current->right, start, end);

    return current;
}
