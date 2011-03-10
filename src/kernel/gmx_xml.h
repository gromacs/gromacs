#ifndef _gmx_xml_h
#define _gmx_xml_h
	
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>

extern void add_xml_int(xmlNodePtr ptr,const char *name,int val);

extern void add_xml_double(xmlNodePtr ptr,const char *name,double val);

extern void add_xml_char(xmlNodePtr ptr,const char *name,const char *val);

extern xmlNodePtr add_xml_child(xmlNodePtr parent,const char *type);

extern xmlNodePtr add_xml_child_val(xmlNodePtr parent,
				    const char *type,const char *value);

extern xmlNodePtr add_xml_comment(xmlDocPtr doc,
				  xmlNodePtr prev,const char *comment);
				  
#endif

#endif
