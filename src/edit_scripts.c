#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "edit_scripts.h"

/* Converts an integer to its octal representation */
char *to_octal_str(int i) {
  char *buf = malloc(16*sizeof(char));
  sprintf(buf, "%o", i);
  return buf;
}

/* Takes in as input two strings, a character representing whether or not they
 * are in the same direction, and the length of the strings and returns an edit
 * script that can convert the reference string to the original string.
 */
char *make_edit_script(char *str, char *ref, char dir, int length){
    char direction = dir == '+' ? '0' : '1';
    bool insert_open = false, subdel_open = false;
    int last_edit = 0;
    char *edit_script = malloc(3*length*sizeof(char));
    int current = 1;
    int i, j;
    char *octal;
    edit_script[0] = direction;
    for(i = 0; i < length; i++){
        if(str[i] == ref[i]){
            insert_open = false;
            subdel_open = false;
        }
        else{ /* mismatch */
            if(ref[i] == '-'){ /* insertion in str relative to ref (i.e., gap in ref) */
                subdel_open = false;
	        if(!insert_open){ /* indicate start of insertion */
	            insert_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 'i';
                    for(j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
	            last_edit = i;
                }
	        edit_script[current++] = str[i];
            }
           /* substitution or deletion in str (represented in script by '-') relative
              to ref */
            else{
                insert_open = false;
                if(!subdel_open){ /* indicate start of subdel */
	            subdel_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 's';
                    for(j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
                    last_edit = i;
                }
	        edit_script[current++] = str[i];
            }
        }
    }
    edit_script = realloc(edit_script, current+1*sizeof(char));
    edit_script[current] = '\0';
    return edit_script;
}

/* reads one edit worth of info (or fails if current decoded char is digit)
   moves pos accordingly */
bool next_edit(char *edit_script, int *pos, struct edit_info *edit){
    int editLength = 0;
    int i = 0;
    if(isdigit(edit_script[(*pos)]) || edit_script[(*pos)] == '\0')
        return false;
    edit->is_subdel = edit_script[(*pos)++] == 's';
    edit->last_dist = 0;
    edit->str = "";
    while(isdigit(edit_script[(*pos)])){
        edit->last_dist *= 8; /* octal encoding */
        edit->last_dist += edit_script[(*pos)++] - '0';
    }
    while(isupper(edit_script[(*pos)+editLength]) ||
                  edit_script[(*pos)+editLength] == '-')
        editLength++;
    edit->str = malloc((editLength+1)*sizeof(char));
    edit->str_length = editLength;
    while (isupper(edit_script[(*pos)]) || edit_script[(*pos)] == '-')
        edit->str[i++] = edit_script[(*pos)++];
    return true;
}

char *read_edit_script(char *edit_script, char *orig, int length){
    char *str = malloc((2*length+1)*sizeof(char));
    int i;
    struct edit_info edit;
    int orig_pos = 0, last_edit_str_len = 0; /* length of last edit str */
    int current = 0;
    int script_pos = 1;

    while(next_edit(edit_script, &script_pos, &edit)){
        /* chunk after previous edit */
        for(i = 0; i < edit.last_dist - last_edit_str_len; i++)
            str[current++] = orig[orig_pos+i];

        /* update position in original string */
        orig_pos += edit.last_dist - last_edit_str_len;

        /* append replacement string in edit script; get rid of dashes */
        for(i = 0; i < edit.str_length; i++)
            if(edit.str[i] != '-')
	        str[current++] = edit.str[i];

        /* skip subdel along original string */
        if (edit.is_subdel) orig_pos += edit.str_length;

        last_edit_str_len = edit.str_length;
    }
    while(orig_pos < length)
        str[current++] = orig[orig_pos++];
    str = realloc(str, current+1*sizeof(char));
    str[current] = '\0';
    printf("Exited read_edit_script\n");
    return str;
}
