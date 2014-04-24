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
    assert(node);

    node->seq = malloc((end-start+2)*sizeof(*(node->seq)));
    assert(node->seq);

    strncpy(node->seq, sequence, (end-start+1));
    node->seq[end-start+1] = '\0';
    node->start = start;
    node->end = end;

    node->left = NULL;
    node->right = NULL;
    return node;
}

/*Creates a new empty cb_range_tree.*/
struct cb_range_tree *cb_range_tree_create(char *seq_name){
    struct cb_range_tree *tree = malloc(sizeof(*tree));
    assert(tree);

    if (seq_name != NULL) {
        tree->seq_name = malloc((strlen(seq_name)+1)*sizeof(*(tree->seq_name)));
        assert(tree->seq_name);

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
    assert(node->seq);

    strncpy(node->seq, sequence, len);
    node->seq[len] = '\0';
    node->start = start;
    node->end = end;
}

/*Create a range node data struct*/
struct cb_range_node_data *cb_range_node_data_create(int start, int end,
                                                     char *sequence){
    int len = -1;
    struct cb_range_node_data *d = malloc(sizeof(*d));
    assert(d);

    assert(end > start);
    assert(sequence != NULL);

    len = end - start;
    d->seq = malloc((len+1)*sizeof(*(d->seq)));
    assert(d->seq);

    strncpy(d->seq, sequence, len);
    d->seq[len] = '\0';
    d->start = start;
    d->end = end;
    return d;
}

/*Update a range node data struct*/
void cb_range_node_data_update(struct cb_range_node_data *d,
                               int start, int end, char *sequence){
    int len = -1;

    assert(end > start);
    assert(sequence != NULL);

    len = end - start;
    if (d->seq != NULL)
        free(d->seq);
    d->seq = malloc((len+1)*sizeof(*(d->seq)));
    assert(d->seq);

    strncpy(d->seq, sequence, len);
    d->seq[len] = '\0';   
    d->start = start;
    d->end = end;
}

/*Free a range node data struct*/
void cb_range_node_data_free(struct cb_range_node_data *d){
    assert(d != NULL);
    if (d->seq != NULL)
        free(d->seq);
    free(d);
}

/*Removes nodes in the tree with ranges adjacent to or overlapping the range
 *being passed in with "start" and "end" in the direction specified by "dir"
 *(-1 = left, 1 = right) and returns the data of the node adjacent to the range
 *passed in that is farthest in the specified direction.
 */
struct cb_range_node_data *remove_adjacent(struct cb_range_node *cur,
                                           int start, int end, int dir,
                                           struct cb_range_node_data *data){
    struct cb_range_node *next, *new_next;
    bool left, adjacent, next_adjacent, next_is_left;
    
    assert(end > start); /*Make sure we are using a valid range*/
    assert(cur != NULL);/*remove_adjacent is never called from a NULL pointer.*/

    left = dir < 0;
    adjacent = start <= cur->end && end >= cur->start;
    next_is_left = adjacent ? left : !left;
    next = next_is_left ? cur->left : cur->right;

    /*If there is no next node, we already found the farthest node in the
      specified direction that is adjacent to the range, so return its data.*/
    if (!next)
        return data;

    next_adjacent = start <= next->end && end >= next->start;

    /*If the next node's range is adjacent to the range, the next node needs to
     *be deleted, so see if the data we are returning needs to be updated and
     *then free the the next node and its subtree in the direction opposite of
     *dir, making it so the next node is now the next node after "next".  Since
     *we don't know if new_next is NULL or has a range adjacent to the range
     *passed in, call the next call to remove_adjacent from the current node.
     */
    if (next_adjacent) {
        if (left && next->start < data->start || !left && next->end > data->end)
            cb_range_node_data_update(data, next->start, next->end, next->seq);

        new_next = left ? next->left : next->right;

        if (next_is_left)
            cur->left = new_next;
        else
            cur->right = new_next;

        if (left)
            next->left = NULL;
        else
            next->right = NULL;

        /*Free the next node and the subtree new_next is not in.*/
        cb_range_node_free(next);

        return remove_adjacent(cur, start, end, dir, data);
    }
    /*If the next node's range is not adjacent to the range, call the next call
      to remove_adjacent from the next node.*/
    else
        return remove_adjacent(next, start, end, dir, data);
}

/*Takes in a tree, an expanded DNA sequence, and its starting and ending
 *indices in the sequence it is from and inserts a node created from those
 *data into the tree and updates the tree so that no two ranges overlap or are
 *adjacent, merging ranges that do overlap or are adjacent.
 */
void cb_range_tree_insert(struct cb_range_tree *tree,
                          char *sequence, int start, int end){
    assert(end > start);
    assert(sequence != NULL);
    tree->root = cb_range_node_insert(tree->root, sequence, start, end);
}


/*Carries out the recursive calls for cb_range_tree_insert.  Takes in the
 *current node and the sequence and range for the node to be inserted and
 *searches the tree to find a location for adding the new range.
 *
 *If cb_range_node_insert does not find a node that overlaps or is adjacent to
 *the range passed in, the range and sequence are used to create a new range
 *node.  Otherwise, the range is merged with the highest node in the tree whose
 *range is adjacent to the range passed in, which is in turn merged with any
 *other nodes that with ranges adjacent to the new range, such that when any
 *insertion is complete no two ondes in the tree have overlapping or adjacent
 *ranges.
 */
struct cb_range_node *cb_range_node_insert(struct cb_range_node *cur,
                                           char *sequence, int start, int end){
    bool adjacent;
    struct cb_range_node_data *last;
    assert(end > start);

    if (!cur)
        return cb_range_node_create(sequence, start, end);
    adjacent = start <= cur->end && end >= cur->start;
    /*If the range being added does not overlap the range of the current node,
     *and is not adjacent to the range of the current node, call
     *cb_range_node_insert on whichever subtree of the current node would
     *contain nodes that would overlap the range.
     */
    if (!adjacent) {
        if (end < cur->start)
            cur->left=cb_range_node_insert(cur->left, sequence, start, end);
        else
            cur->right=cb_range_node_insert(cur->right, sequence, start, end);
    }
    else {
        if (start < cur->start) {
fprintf(stderr, "Starting left merge\n");
            /*If the new range overlaps the current node's range to the left,
             *update the current node's range and sequence to merge the current
             *node's range and sequence with the new range and sequence.
             */
            int left_merge_end = end < cur->end ? end : cur->end;
            cb_range_node_update(cur,
                             cb_range_merge(sequence, start, left_merge_end,
                                            cur->seq, cur->start, cur->end),
                             start, cur->end);
            /*Search the tree for the leftmost node that overlaps the new range
              to the left, deleting any node overlapping the new range.*/
            last = cb_range_node_data_create(start, end, sequence);
            remove_adjacent(cur, start, end, -1, last);

            /*If remove_adjacent found a node that overlaps the new range, merge
              it with the current node.*/
            if (last->start < start)
                cb_range_node_update(cur,
                                 cb_range_merge(last->seq,last->start,last->end,
                                                cur->seq, cur->start, cur->end),
                                 last->start, cur->end);
        }
        if (end > cur->end) {
fprintf(stderr, "Starting right merge\n");
            /*If the new range overlaps the current node's range to the right,
             *update the current node's range and sequence to merge the current
             *node's range and sequence with the new range and sequence.
             */
            cb_range_node_update(cur,
                             cb_range_merge(cur->seq, cur->start, cur->end,
                                            sequence, start, end),
                             cur->start, end);
            /*Search the tree for the rightmost node that overlaps the new range
              to the right, deleting any node overlapping the new range.*/
            last = cb_range_node_data_create(start, end, sequence);
            remove_adjacent(cur, start, end, 1, last);
fprintf(stderr, "last->end = %d\n", last->end);

            /*If remove_adjacent found a node that overlaps the new range, merge
              it with the current node.*/
            if (last->end > end)
                cb_range_node_update(cur,
                                 cb_range_merge(cur->seq, cur->start, cur->end,
                                            last->seq, last->start, last->end),
                                 cur->start, last->end);
        }
    }
    return cur;
}

/*Takes in the sequences and starting and ending indices of two sections of an
  original sequence and merges the sequences together.*/
char *cb_range_merge(char *left_seq, int left_start, int left_end,
                     char *right_seq, int right_start, int right_end){
fprintf(stderr, "Merging %d-%d (%s) and %d-%d (%s)\n", left_start, left_end, left_seq,
                                                       right_start, right_end, right_seq);
    int index = 0, i = 0;
    char *merged = NULL;

    /*Make sure we are merging valid ranges*/
    assert(left_end > left_start && right_end > right_start);
    assert(left_start <= right_end && right_end >= left_end);
    assert(left_seq != NULL && right_seq != NULL);

    merged = malloc((right_end-left_start+1)*sizeof(*merged));
    assert(merged);

    for (i = 0; i < right_start - left_start; i++)
        merged[index++] = left_seq[i];
    for (i = 0; i < right_end - right_start; i++)
        merged[index++] = right_seq[i];
    merged[right_end-left_start] = '\0';
fprintf(stderr, "Result of merge is %s\n", merged);
    return merged;
}

/*Takes in a range tree and the starting and ending indices of a range in a DNA
 *sequence and returns a node in the tree with a range that overlaps the range
 *if there is one or a null pointer if there is no node in the tree overlapping
 *the range.
 */
struct cb_range_node *cb_range_tree_find(struct cb_range_tree *tree,
                                                int start, int end){
    assert(tree);
    assert(end > start);
    return cb_range_node_find(tree->root, start, end);
}

/*The recursive function called from cb_range_tree_find.  Takes in the current
 *node and the range passed in, making another recursive call if the node's
 *range does not overlap the range passed in, returning the current node if the
 *current node's range does overlap the range passed in, or returning a null
 *pointer if the current node is a null pointer.
 */
struct cb_range_node *cb_range_node_find(struct cb_range_node *node,
                                                int start, int end){
    if(node == NULL)
        return NULL;
    else if(start <= node->end && end >= node->start)
        return node;
    else
        return cb_range_node_find((start < node->start?node->left:node->right),
                                                                   start, end);
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
 *the file in FASTA format.
 */
void cb_range_node_output_fasta(struct cb_range_node *node,
                                struct cb_range_tree *tree, void *out){
    fprintf((FILE *)out, "> %s\n%s\n", tree->seq_name, node->seq);
}

/*Uses cb_range_tree_traverse and cb_range_node_output to output the sequences
  in the range tree passed in to the file pointer passed in in FASTA format.*/
void cb_range_tree_output_fasta(struct cb_range_tree *tree, FILE *out){
    cb_range_tree_traverse(tree, cb_range_node_output_fasta, (void *)out);
}

/*A function to be passed into cb_range_tree_traverse that takes in a range
 *tree node, its tree, and a FILE pointer and outputs the node's range to
 *and its sequence to the specified file pointer.
 */
void cb_range_node_output_data(struct cb_range_node *node,
                          struct cb_range_tree *tree, void *out){
    (void) tree;
    fprintf((FILE *)out, "> %d-%d... %s\n", node->start, node->end, node->seq);
}

/*Uses cb_range_tree_traverse and cb_range_node_output to output the sequences
  and ranges in the range tree passed in to the file pointer passed in.*/
void cb_range_tree_output_data(struct cb_range_tree *tree, FILE *out){
    fprintf(stderr, "Outputting range tree\n");
    cb_range_tree_traverse(tree, cb_range_node_output_data, (void *)out);
    fprintf(stderr, "Finished outputting range tree\n");
}

