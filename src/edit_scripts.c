#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "edit_scripts.h"

char *to_octal_str(int i) {
  char *buf = malloc(16*sizeof(char));
  sprintf(buf, "%o", i);
  return buf;
}

char *make_edit_script(char **axt, char dir, int length){ 
    char direction = dir == '+' ? '0' : '1';
    char *str = axt[0]; /* strings are aligned and include gaps */
    char *ref = axt[2];
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
