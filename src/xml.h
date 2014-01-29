#include <stdlib.h>
#include "ds.h"

#ifndef __CABLASTP_SEARCH_XML_H__
#define __CABLASTP_SEARCH_XML_H__
struct blast{
    char *xml_name;
    struct DSLinkedList *hits;
};

struct hit{
    char *xml_name;
    int32_t num;
    int32_t accession;
    struct DSLinkedList *hsps; 
};

struct hsp{
    char *xml_name;
    int32_t num;
    double evalue;
    int32_t query_from;
    int32_t query_to;
    int32_t hit_from;
    int32_t hit_to;
};
#endif