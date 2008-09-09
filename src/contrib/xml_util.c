#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_fatal.h"
#include "xml_util.h"
	
int find_elem(char *name,int nr,const char *names[])
{
  int i;
  
  for(i=0; (i<nr); i++)
    if (strcmp(name,names[i]) == 0) 
      break;
  if (i == nr)
    gmx_fatal(FARGS,"Unknown element name %s",name);
    
  return i;
}

void add_xml_int(xmlNodePtr ptr,const char *name,int val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%d",val);
  if (xmlSetProp(ptr,(xmlChar *)name,buf) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

void add_xml_double(xmlNodePtr ptr,const char *name,double val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%g",val);
  if (xmlSetProp(ptr,(xmlChar *)name,buf) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

void add_xml_char(xmlNodePtr ptr,const char *name,char *val)
{
  if (xmlSetProp(ptr,(xmlChar *)name,(xmlChar *)val) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

xmlNodePtr add_xml_child(xmlNodePtr parent,const char *type)
{
  xmlNodePtr child;
  
  if ((child = xmlNewChild(parent,NULL,(xmlChar *)type,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating element",(char *)type);
  
  return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr doc,
			   xmlNodePtr prev,char *comment)
{
  xmlNodePtr comm,ptr;
  
  if ((comm = xmlNewComment((xmlChar *)comment)) == NULL)
    gmx_fatal(FARGS,"Creating doc comment element","");
  ptr = prev;
  while (ptr->next != NULL)
    ptr=ptr->next;
  ptr->next    = comm;
  comm->prev   = ptr;
  comm->doc    = doc;
  
  return comm;
}

