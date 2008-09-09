#ifndef _xml_util_h
#define _xml_util_h
	
#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_fatal.h"
	
extern int find_elem(char *name,int nr,const char *names[]);

extern void add_xml_int(xmlNodePtr ptr,const char *name,int val);

extern void add_xml_double(xmlNodePtr ptr,const char *name,double val);

extern void add_xml_char(xmlNodePtr ptr,const char *name,char *val);

extern xmlNodePtr add_xml_child(xmlNodePtr parent,const char *type);

extern xmlNodePtr add_xml_comment(xmlDocPtr doc,
				  xmlNodePtr prev,char *comment);

#endif
