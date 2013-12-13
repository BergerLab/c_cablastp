#ifndef __CABLASTP_EDITSCRIPTS_H__
#define __CABLASTP_EDITSCRIPTS_H__

struct edit_info{
    bool is_subdel;
    int last_dist;
    char *str;
    int str_length;
};

char to_half_byte(char c);
char half_byte_to_char(char h);
char *to_octal_str(int i);
char *make_edit_script(char *str, char *ref, char dir, int length);
char *read_edit_script(char *edit_script, char *orig, int length);
bool next_edit(char *edit_script, int *pos, struct edit_info *edit);
#endif
