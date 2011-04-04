/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <string.h>
#include "gmx_fatal.h"
#include "smalloc.h"
#include "macros.h"
#include "futil.h"
#include "gmx_qhop_parm.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
	
extern int xmlDoValidityCheckingDefaultValue;
	
#define NN(x) (NULL != (x))

static const char *xmltypes[] = { 
  NULL, 
  "XML_ELEMENT_NODE",
  "XML_ATTRIBUTE_NODE",
  "XML_TEXT_NODE",
  "XML_CDATA_SECTION_NODE",
  "XML_ENTITY_REF_NODE",
  "XML_ENTITY_NODE",
  "XML_PI_NODE",
  "XML_COMMENT_NODE",
  "XML_DOCUMENT_NODE",
  "XML_DOCUMENT_TYPE_NODE",
  "XML_DOCUMENT_FRAG_NODE",
  "XML_NOTATION_NODE",
  "XML_HTML_DOCUMENT_NODE",
  "XML_DTD_NODE",
  "XML_ELEMENT_DECL",
  "XML_ATTRIBUTE_DECL",
  "XML_ENTITY_DECL",
  "XML_NAMESPACE_DECL",
  "XML_XINCLUDE_START",
  "XML_XINCLUDE_END"
};
#define NXMLTYPES asize(xmltypes)
	
enum { 
  exmlQHOPS,
  exmlQHOP, exmlDONOR, exmlACCEPTOR,
  exmlPARAM, exmlNAME, exmlVALUE, 
  exmlUNIT,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "qhops",
  "qhop", "donor", "acceptor", "parameter", 
  "name", "value", "unit" 
};

typedef struct {
  int        nqh;
  gmx_qhop_t *gqh;
} t_xmlrec;

static int find_elem(char *name,int nr,const char *names[])
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

static char *sp(int n, char buf[], int maxindent)
{
  int i;
  if(n>=maxindent)
    n=maxindent-1;
  
  /* Don't indent more than maxindent characters */
  for(i=0; (i<n); i++)
    buf[i] = ' ';
  buf[i] = '\0';
  
  return buf;
}

static void qhop_process_attr(FILE *fp,xmlAttrPtr attr,int parent,
			      int elem,int indent,gmx_qhop_t qht)
{
  char *attrname,*attrval;
  char buf[100];
  int  i,kkk,eprop;
  char *xbuf[exmlNR];
  
  for(i=0; (i<exmlNR); i++)
    xbuf[i] = NULL;
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((gmx_strcasecmp(attrname,s) == 0) && (attrval != NULL))
    kkk = find_elem(attrname,exmlNR,exml_names);
    if (attrval != NULL)
      xbuf[kkk] = strdup(attrval);
      
    if (fp)
      fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent,buf,99),
	      attrname,attrval);
    attr = attr->next;
#undef atest
  }
  
  switch (elem) {
  case exmlQHOP:
    if (NN(xbuf[exmlDONOR]) && NN(xbuf[exmlACCEPTOR])) {
      gmx_qhop_set_donor(qht,xbuf[exmlDONOR]);
      gmx_qhop_set_acceptor(qht,xbuf[exmlACCEPTOR]);
    }
    break;
  case exmlPARAM:
    if (NN(xbuf[exmlNAME]) && NN(xbuf[exmlUNIT]) && NN(xbuf[exmlVALUE])) {
      gmx_qhop_add_param(qht,xbuf[exmlNAME],xbuf[exmlVALUE],
			 xbuf[exmlUNIT]);
    }
    break;
  default:
    break;
  }
  for(i=0; (i<exmlNR); i++)
    if (NN(xbuf[i]))
      sfree(xbuf[i]);
}

static void qhop_process_element(FILE *fp,xmlNodePtr tree,int parent,
				 int indent,t_xmlrec *xml)
{
  int elem;
  char buf[100];
  
  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  if (elem == exmlQHOP) {
    xml->nqh++;
    srenew(xml->gqh,xml->nqh);
    xml->gqh[xml->nqh-1] = gmx_qhop_init();
  }
  if (elem != exmlQHOPS)
    qhop_process_attr(fp,tree->properties,parent,
		      elem,indent+2,xml->gqh[xml->nqh-1]);
}
 
static void qhop_process_tree(FILE *fp,xmlNodePtr tree,int parent,
			      int indent,t_xmlrec *xml)
{
  char buf[100];
  int  elem;
  
  while (tree != NULL) {
    if (fp) {
      if ((tree->type > 0) && (tree->type < NXMLTYPES))
	fprintf(fp,"Node type %s encountered with name %s\n",
		xmltypes[tree->type],(char *)tree->name);
      else
	fprintf(fp,"Node type %d encountered\n",tree->type);
    }
    
    switch (tree->type) {
    case XML_ELEMENT_NODE:
      qhop_process_element(fp,tree,parent,indent+2,xml);
      
      if (tree->children) {
	elem = find_elem((char *)tree->name,exmlNR,exml_names);
	qhop_process_tree(fp,tree->children,elem,indent+2,xml);
      }
      break;
    default:
      break;
    }
    tree = tree->next;
  }
}

gmx_qhop_t *
gmx_qhops_read(char *fn,int *nqhop)
{
  xmlDocPtr     doc;
  int           i,npd;
  t_xmlrec      *xml;
  const char *db="qhops.dat";
  gmx_bool fna=FALSE;
  
  xmlDoValidityCheckingDefaultValue = 0;
  if (NULL == fn) 
  {
    fn = (char *)gmxlibfn(db);
    fna=TRUE;
  }
  if ((doc = xmlParseFile(fn)) == NULL) {
		fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.",
				fn);
    exit(1);
  }

  snew(xml,1);
  qhop_process_tree(NULL,doc->children,0,0,xml);
  
  xmlFreeDoc(doc);
  if (fna)
      sfree(fn);
  
  *nqhop = xml->nqh;
  
  return xml->gqh;
}

static void add_xml_qhop(xmlNodePtr parent,gmx_qhop_t qht)
{
  xmlNodePtr ptr,child,grandchild,comp;
  char   *name,*type,*value,*unit;
  
  ptr = add_xml_child(parent,exml_names[exmlQHOP]);
  add_xml_char(ptr,exml_names[exmlDONOR],gmx_qhop_get_donor(qht));
  add_xml_char(ptr,exml_names[exmlACCEPTOR],gmx_qhop_get_acceptor(qht));
  
  while (gmx_qhop_get_param(qht,&name,&value,&unit) == 1) {
    child = add_xml_child(ptr,exml_names[exmlPARAM]);
    add_xml_char(child,exml_names[exmlNAME],name);
    add_xml_char(child,exml_names[exmlVALUE],value);
    add_xml_char(child,exml_names[exmlUNIT],unit);
    sfree(name);
    sfree(value);
    sfree(unit);
  }
}

void gmx_qhops_write(char *fn,int nqhop,gmx_qhop_t qht[])
{
  xmlDocPtr  doc;
  xmlDtdPtr  dtd;
  xmlNodePtr myroot;
  int        i,nmt;
  xmlChar    *libdtdname,*dtdname,*gmx;
  
  gmx        = (xmlChar *) "qhops";
  dtdname    = (xmlChar *) "qhops.dtd";
  libdtdname = dtdname;
  
  if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
    gmx_fatal(FARGS,"Creating XML document","");
    
  if ((dtd = xmlCreateIntSubset(doc,dtdname,libdtdname,dtdname)) == NULL)
    gmx_fatal(FARGS,"Creating XML DTD","");
    
  if ((myroot = xmlNewDocNode(doc,NULL,gmx,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating root element","");
  dtd->next = myroot;
  myroot->prev = (xmlNodePtr) dtd;
    
  /* Add molecule definitions */
  for(i=0; (i<nqhop); i++)
    add_xml_qhop(myroot,qht[i]);

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving file",fn);
  xmlFreeDoc(doc);
}



#else

gmx_qhop_t gmx_qhops_read(char *fn,int *nqhop)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function gmx_qhops_read to work");
  return NULL;
}

void gmx_qhops_write(char *fn,int nqhop,gmx_qhop_t qht)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function gmx_qhops_write to work");
}

#endif
