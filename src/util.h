#ifndef __CABLAST_UTIL_H__
#define __CABLAST_UTIL_H__

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

char *
trim(char *s, const char *totrim);

char *
trim_space(char *s);

int32_t
readline(FILE *f, char **line);

int32_t
num_cpus();

char *
str_slice(char *str, int32_t start, int32_t end);

bool
is_complete_overlap(char *s1, char *s2);

#endif
