#ifndef _gmx_xml_h
#define _gmx_xml_h
	
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>

extern void add_xml_int(xmlNodePtr ptr,xmlChar *name,int val);

extern void add_xml_double(xmlNodePtr ptr,xmlChar *name,double val);

extern void add_xml_char(xmlNodePtr ptr,xmlChar *name,char *val);

extern xmlNodePtr add_xml_child(xmlNodePtr parent,char *type);

extern xmlNodePtr add_xml_comment(xmlDocPtr doc,
				  xmlNodePtr prev,xmlChar *comment);
				  
#endif

#endif
