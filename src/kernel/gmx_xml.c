#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "gmx_fatal.h"


#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_xml.h"
	
void add_xml_int(xmlNodePtr ptr,xmlChar *name,int val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%d",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

void add_xml_double(xmlNodePtr ptr,xmlChar *name,double val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%g",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

void add_xml_char(xmlNodePtr ptr,xmlChar *name,char *val)
{
  if (xmlSetProp(ptr,name,(xmlChar *)val) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

xmlNodePtr add_xml_child(xmlNodePtr parent,char *type)
{
  xmlNodePtr child;
  
  if ((child = xmlNewChild(parent,NULL,(xmlChar *)type,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating element %s",(char *)type);
  
  return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr doc,
			   xmlNodePtr prev,xmlChar *comment)
{
  xmlNodePtr comm,ptr;
  
  if ((comm = xmlNewComment(comment)) == NULL)
    gmx_fatal(FARGS,"Creating doc comment element");
  ptr = prev;
  while (ptr->next != NULL)
    ptr=ptr->next;
  ptr->next    = comm;
  comm->prev   = ptr;
  comm->doc    = doc;
  
  return comm;
}

#endif
