/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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

#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_fatal.h"
#include "macros.h"
#include "grompp.h"
#include "smalloc.h"
#include "molprop.h"
#include "xml_util.h"

extern int xmlDoValidityCheckingDefaultValue;
	
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
  exmlMOLECULES,
  exmlMOLECULE, exmlFORMULA, exmlMOLNAME, exmlWEIGHT,
  exmlCATEGORY, exmlCATNAME,
  exmlPROPERTY, exmlPNAME,   exmlPVALUE, exmlPERROR, 
  exmlPMETHOD, exmlPREFERENCE,
  exmlCOMPOSITION, exmlCOMPNAME, exmlCATOM, exmlC_NAME, exmlC_NUMBER,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "molecules",
  "molecule", "formula", "molname", "weight", 
  "category", "catname",
  "property", "pname", "pvalue", "perror", "pmethod", "preference",
  "composition", "compname", "catom", "cname", "cnumber"
};

typedef struct {
  int nmolprop;
  gmx_molprop_t *mpt;
} t_xmlrec;

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

static char *process_attr(FILE *fp,xmlAttrPtr attr,int elem,
			  int indent,gmx_molprop_t mpt)
{
  char *attrname,*attrval;
  char buf[100];
  char *pname=NULL,*preference=NULL,*pmethod=NULL,*pvalue=NULL,*perr=NULL;
  char *cname=NULL,*cnumber=NULL;
  
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
    switch (elem) {
    case exmlMOLECULE: 
      if (atest("formula")) 
	gmx_molprop_set_formula(mpt,attrval);
      else if (atest("molname")) 
	gmx_molprop_set_molname(mpt,attrval);
      else if (atest("weight")) 
	gmx_molprop_set_weight(mpt,atof(attrval));
      else
	gmx_fatal(FARGS,"Unknown attribute %s for molecule",attrname);
      break;
      
    case exmlCATEGORY:
      if (atest("catname")) 
	gmx_molprop_add_category(mpt,attrval);
      else
	gmx_fatal(FARGS,"Unknown attribute %s for category",attrname);
      break;
      
    case exmlPROPERTY: {
      if (atest("pname")) 
	pname = strdup(attrval);
      else if (atest("pmethod")) 
	pmethod = strdup(attrval);
      else if (atest("preference")) 
	preference = strdup(attrval);
      else if (atest("pvalue")) 
	pvalue = strdup(attrval);
      else if (atest("perror")) 
	perr = strdup(attrval);
      else
	gmx_fatal(FARGS,"Unknown attribute %s for property",attrname);
      break;
    }
    case exmlCOMPOSITION: {
      if (atest("compname")) 
	gmx_molprop_add_composition(mpt,attrval);
      else
	gmx_fatal(FARGS,"Unknown attribute %s for catom",attrname);
      break;
    }
    case exmlCATOM:
      if (atest("cname")) 
	cname = strdup(attrval);
      else if (atest("cnumber")) 
	cnumber = strdup(attrval);
      else
	gmx_fatal(FARGS,"Unknown attribute %s for catom",attrname);
      break;
      
    default:
      break;
    }
    if (fp)
      fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent,buf,99),
	      attrname,attrval);
    attr = attr->next;
#undef atest
  }
  /* Done processing attributes for this element. Let's see if we still need
     to interpret them.
  */
  if (pname && preference && pmethod && pvalue)
    gmx_molprop_add_property(mpt,0,pname,atof(pvalue),perr ? atof(perr) : 0.0,
			     pmethod,preference);
  if (cname && cnumber)
    gmx_molprop_add_composition_atom(mpt,NULL,cname,atoi(cnumber));
}

static void process_element(FILE *fp,xmlNodePtr tree,int indent,t_xmlrec *xml)
{
  int elem;
  char buf[100];
  
  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  switch (elem) {
  case exmlMOLECULES:
    break;
  case exmlMOLECULE:
    xml->nmolprop++;
    srenew(xml->mpt,xml->nmolprop);
    xml->mpt[xml->nmolprop-1] = gmx_molprop_init();
    break;
    /* The items below are handled when treating attributes */
  case exmlCATEGORY:
    break;
  case exmlPROPERTY: 
    break;
  case exmlCOMPOSITION:
    break;
  case exmlCATOM: 
    break;
  default:
    break;
  }
  if (elem != exmlMOLECULES)
    process_attr(fp,tree->properties,elem,indent+2,xml->mpt[xml->nmolprop-1]);
}
 
static void process_tree(FILE *fp,xmlNodePtr tree,int indent,t_xmlrec *xml)
{
  char          buf[100];
  
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
      process_element(fp,tree,indent+2,xml);
      
      if (tree->children)
	process_tree(fp,tree->children,indent+2,xml);
      break;
    default:
      break;
    }
    tree = tree->next;
  }
}

gmx_molprop_t *gmx_molprops_read(char *fn,int *nmolprop)
{
  xmlDocPtr     doc;
  int           i,npd;
  t_xmlrec      xml;
    
  xmlDoValidityCheckingDefaultValue = 0;
  if ((doc = xmlParseFile(fn)) == NULL) {
    fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.",
	    fn);
    exit(1);
  }

  xml.nmolprop = 0;
  xml.mpt = NULL;
  process_tree(NULL,doc->children,0,&xml);

  xmlFreeDoc(doc);
  
  *nmolprop = xml.nmolprop;
  
  return xml.mpt;
}

static void add_xml_molprop(xmlNodePtr parent,gmx_molprop_t mpt)
{
  xmlNodePtr ptr,child,grandchild,comp;
  int    i,eMP,cnumber;
  char   *tmp,*prop_name,*prop_method,*prop_reference;
  char   *catom;
  double value,error;
  
  
  ptr = add_xml_child(parent,exml_names[exmlMOLECULE]);
  add_xml_char(ptr,exml_names[exmlMOLNAME],gmx_molprop_get_molname(mpt));
  add_xml_char(ptr,exml_names[exmlFORMULA],gmx_molprop_get_formula(mpt));
  add_xml_double(ptr,exml_names[exmlWEIGHT],gmx_molprop_get_weight(mpt));
  
  while ((tmp = gmx_molprop_get_category(mpt)) != NULL) {
    child = add_xml_child(ptr,exml_names[exmlCATEGORY]);
    add_xml_char(child,exml_names[exmlCATNAME],tmp);
  }
	
  while (gmx_molprop_get_property(mpt,&eMP,&prop_name,&value,&error,
				  &prop_method,&prop_reference) == 1) {
    child = add_xml_child(ptr,exml_names[exmlPROPERTY]);
    add_xml_char(child,exml_names[exmlPNAME],prop_name);
    add_xml_double(child,exml_names[exmlPVALUE],value);
    add_xml_double(child,exml_names[exmlPERROR],error);
    add_xml_char(child,exml_names[exmlPMETHOD],prop_method);
    add_xml_char(child,exml_names[exmlPREFERENCE],prop_reference);
    sfree(prop_name);
    sfree(prop_method);
    sfree(prop_reference);
  }
  
  while ((tmp = gmx_molprop_get_composition(mpt)) != NULL) {
    child = add_xml_child(ptr,exml_names[exmlCOMPOSITION]);
    add_xml_char(child,exml_names[exmlCOMPNAME],tmp);
    while ((gmx_molprop_get_composition_atom(mpt,&catom,&cnumber)) == 1) {
      grandchild = add_xml_child(child,exml_names[exmlCATOM]);
      add_xml_char(grandchild,exml_names[exmlC_NAME],catom);
      add_xml_int(grandchild,exml_names[exmlC_NUMBER],cnumber);
      sfree(catom);
    }
  }
}

void gmx_molprops_write(char *fn,int nmolprop,gmx_molprop_t mpt[])
{
  xmlDocPtr  doc;
  xmlDtdPtr  dtd;
  xmlNodePtr myroot;
  int        i,nmt;
  xmlChar    *libdtdname,*dtdname,*gmx;
  
  gmx        = (xmlChar *) "molecules";
  dtdname    = (xmlChar *) "molprops.dtd";
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
  for(i=0; (i<nmolprop); i++)
    add_xml_molprop(myroot,mpt[i]);

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving file",fn);
  xmlFreeDoc(doc);
}



