#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "DNAutils.h"
#include "edit_scripts.h"

/* Converts an integer to its octal representation */
char *to_octal_str(int i) {
    char *buf = malloc(16*sizeof(char));
    sprintf(buf, "%o", i);
    return buf;
}

/* Converts a character in an edit script to a half-byte */
char to_half_byte(char c){
    switch (c) {
        case '0': return (char)0;
        case '1': return (char)1;
        case '2': return (char)2;
        case '3': return (char)3;
        case '4': return (char)4;
        case '5': return (char)5;
        case '6': return (char)6;
        case '7': return (char)7;
        case 'A': return (char)8;
        case 'C': return (char)9;
        case 'G': return (char)10;
        case 'T': return (char)11;
        case '-': return (char)12;
        case 'i': return (char)13;
        case 's': return (char)14;
        /*For direction bytes with a 1 as the first bit to indicate that
          the script is from a match, move the 1 to the leftmost bit of the
          half-byte*/
        case '0' | (char)0x80: return (char)8;
        case '9' | (char)0x80: return (char)9;
        /*'N' is represented as the half byte 1111*/
        default:  return (char)15;
    }
}

/* Converts a half-byte to an edit script character */
char half_byte_to_char(char h){
    if (h < 8) /*h is an octal digit*/
        return '0' + h;
    switch (h) {
        case (char)8:  return 'A';
        case (char)9:  return 'C';
        case (char)10: return 'G';
        case (char)11: return 'T';
        case (char)12: return '-';
        case (char)13: return 'i';
        case (char)14: return 's';
        default:       return 'N';
    }
}

/* Converts the string for an edit script to half-byte format */
char *edit_script_to_half_bytes(char *edit_script){
    int i = 0;
    int length = 0;
    int odd;
    char *half_bytes;
    while (edit_script[length] != '\0')
        length++;
    odd = length % 2;
    length = length / 2 + odd;
    half_bytes = malloc(length*sizeof(*half_bytes));
 
    while (edit_script[i] != '\0') {
        if (i%2 == 0)
            half_bytes[i/2] = (char)0;
        else
            half_bytes[i/2] <<= 4;
        /*Insert the current half byte into the current byte*/
        half_bytes[i/2] |= to_half_byte(edit_script[i]);
        i++;
    }
    /*If the length of the edit script is odd, to signify the end of the edit
     *script add a 1 to signify the end of the script since this is incorrect
     *edit script syntax.
     */
    if (odd) {
        half_bytes[i/2] <<= 4;
        half_bytes[i/2] |= (char)1;
    }
    return half_bytes;
}

/* Converts an edit script in half-byte format to ASCII */
char *half_bytes_to_ASCII(char *half_bytes, int length){
    int i = 0;
    char *edit_script = malloc((length+1)*sizeof(char));
    for (i = 0; i < length; i++) {
        if (i % 2 == 0) { /* Copy the left half-byte of the current byte */
            char left = half_bytes[i/2] & (((char)15) << 4);
            left >>= 4;
            left &= (char)15;
            edit_script[i] = half_byte_to_char(left);
            /*Handle the direction byte so that the match bit of the
              half-byte is added to the direction byte separately from the
              direction bit*/
            if (i == 0) {
                left &= (char)1;
                edit_script[i] = half_byte_to_char(left);
                if ((half_bytes[0] & (char)0x80) != 0)
                    edit_script[i] |= (char)0x80;
            }
        }
        else /* Copy the right half-byte of the current byte */
            edit_script[i] = half_byte_to_char(half_bytes[i/2] & (char)15);
    }
    edit_script[length] = '\0';
    return edit_script;
}

/* Takes in as input two strings, a bool representing whether or not they are
 * in the same direction, and the length of the strings and returns an edit
 * script that can convert the reference string to the original string.
 */
char *make_edit_script(char *str, char *ref, bool dir, int length){
    /*direction has its first bit set to 1 to indicate that the edit script
      was made from a match*/
    char direction = (dir ? '0' : '1');
    direction |= ((char)0x80);
    bool insert_open = false, subdel_open = false;
    int last_edit = 0;
    char *edit_script = malloc(3*length*sizeof(char));
    int current = 1;
    int i, j;
    char *octal;
    edit_script[0] = direction;
    for (i = 0; i < length; i++) {
        if (str[i] == ref[i]) {
            insert_open = false;
            subdel_open = false;
        }
        else { /* mismatch */
            if (ref[i] == '-') { /* insertion in str relative to ref (i.e., gap in ref) */
                subdel_open = false;
	        if (!insert_open) { /* indicate start of insertion */
	            insert_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 'i';
                    for (j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
	            last_edit = i;
                }
                edit_script[current++] = str[i];
            }
           /* substitution or deletion in str (represented in script by '-') relative
              to ref */
            else {
                insert_open = false;
                if (!subdel_open) { /* indicate start of subdel */
	            subdel_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 's';
                    for (j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
                    last_edit = i;
                }
                edit_script[current++] = str[i];
            }
        }
    }
    edit_script = realloc(edit_script, (current+1)*sizeof(char));
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

/*Takes in as input an edit script in ASCII format, a sequence to read, and
  the length of the sequence and applies the edit script to the sequence to
  produce a new sequence.*/
char *read_edit_script(char *edit_script, char *orig, int length){
    char *str = malloc((2*length+1)*sizeof(char));
    int i;
    struct edit_info edit;
    int orig_pos = 0, last_edit_str_len = 0; /* length of last edit str */
    int current = 0;
    int script_pos = 1;
    char *s = str;
    /*char *original = orig;*//*(edit_script[0] == '0' ? orig : string_revcomp(orig, -1));*/

    while (next_edit(edit_script, &script_pos, &edit)) {
        /* chunk after previous edit */
        for (i = 0; i < edit.last_dist - last_edit_str_len; i++)
            str[current++] = orig[orig_pos+i];

        /* update position in original string */
        orig_pos += edit.last_dist - last_edit_str_len;

        /* append replacement string in edit script; get rid of dashes */
        for (i = 0; i < edit.str_length; i++)
            if (edit.str[i] != '-')
	        str[current++] = edit.str[i];

        /* skip subdel along original string */
        if (edit.is_subdel) orig_pos += edit.str_length;

        last_edit_str_len = edit.str_length;
    }
    while (orig_pos < length)
        str[current++] = orig[orig_pos++];
    str = realloc(str, current+1*sizeof(char));
    str[current] = '\0';
    if ((edit_script[0] & ((char)0x7f)) == '1')
        str = string_revcomp(s, -1);
    return str;
}

char *no_dashes(char *sequence){
    int length;
    int bases = 0;
    int i, j;
    j = 0;
    for (length = 0; sequence[length] != '\0'; length++)
        if (sequence[length] != '-')
            bases++;
    char *n = malloc((bases+1)*sizeof(char));
    n[bases] = '\0';
    for (i = 0; i < length; i++)
        if (sequence[i] != '-')
            n[j++] = sequence[i];
    return n;
}
