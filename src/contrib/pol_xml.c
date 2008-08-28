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
#include "gmx_fatal.h"
#include "macros.h"
#include "grompp.h"
#include "tune_pol.h"
#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>

extern int xmlDoValidityCheckingDefaultValue;
	
typedef struct {
  char *miller;
  char *spoel;
} t_miller2spoel;

#define NM2S 4

static t_miller2spoel miller2spoel[NM2S] = {
  { "CBR",  "C20" },
  { "OPI2", "O20" },
  { "SPI2", "S20" },
  { "NPI2", "N20" }
};

typedef struct {
  char   *pname;
  double pvalue,perror;
  char   *psource;
  char   *preference;
} t_property;

typedef struct {
  char *cname;
  int  cnumber;
} t_catom;

typedef struct {
  char    *compname;
  int     ncatom;
  t_catom *catom;
} t_composition;

typedef struct {
  char          *formula,*molname,*reference;
  double        weight;
  int           nproperty;
  t_property    *property;
  int           ncomposition;
  t_composition *composition;
  int           ncategory;
  char          **category;
} t_molprop2;

typedef struct {
  int        npd;
  t_molprop  **pd;
  t_molprop2 **pd2;
} t_xmlrec;
	
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
  exmlPSOURCE, exmlPREFERENCE,
  exmlCOMPOSITION, exmlCOMPNAME, exmlCATOM, exmlC_NAME, exmlC_NUMBER,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "molecules",
  "molecule", "formula", "molname", "weight", 
  "category", "catname",
  "property", "pname", "pvalue", "perror", "psource", "preference",
  "composition", "compname", "catom", "cname", "cnumber"
};

static int find_elem(char *name,int nr,const char *names[])
{
  int i;
  
  for(i=0; (i<nr); i++)
    if (strcmp(name,names[i]) == 0) 
      break;
  if (i == nr)
    fatal_error("Unknown element name %s",name);
    
  return i;
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

static void set_miller(t_molprop *mp,char *milnm,int number,int nh)
{
  int i;
  
  mp->emil_comp[emlH] += nh*number;
  for(i=0; (i<emlNR); i++) {
    if (strcmp(miller[i].name,milnm) == 0) {
      mp->emil_comp[i] += number;
      return;
    }
  }
  if (i == emlNR) 
    fatal_error("No such Miller type",milnm);
}

static void set_bosque(t_molprop *mp,char *bosnm,int number,int nh)
{
  int i;

  mp->elem_comp[eelemH] += nh*number;
  for(i=0; (i<eelemNR); i++) {
    if (strcmp(bosque[i].name,bosnm) == 0) {
      mp->elem_comp[i] += number;
      return;
    }
  } 
  if (i == eelemNR) 
    fatal_error("No such Bosque type",bosnm);
}

static void copy_molprop(t_molprop *dst,t_molprop2 *src,int update_bm)
{
  int i,j,k,n,nh,imiller;
  char *ptr;
    
  if (src->formula) 
    dst->formula    = strdup(src->formula);
  if (src->molname)
    dst->molname    = strdup(src->molname);
  dst->weight     = src->weight;
  
  for(i=0; (i<src->ncategory); i++) {
    for(j=0; (j<ecNR); j++) {
      if (strcasecmp(src->category[i],ec_name[j]) == 0) {
	dst->category[j] = 1;
      }
    }
  }

  for(i=0;(i<src->nproperty);i++){
    if (strcasecmp(src->property[i].psource,"experiment") == 0) {
      dst->nexperiment++;
      srenew(dst->experiment,dst->nexperiment);
      srenew(dst->reference,dst->nexperiment);
      srenew(dst->pname,dst->nexperiment);
      dst->experiment[dst->nexperiment-1] = src->property[i].pvalue;
      dst->reference[dst->nexperiment-1]  = strdup(src->property[i].preference);
      dst->pname[dst->nexperiment-1] = strdup(src->property[i].pname);
    }
    else {
      for(j=0; (j<eqmNR); j++) 
	if (strcasecmp(lbasis[j],src->property[i].psource) == 0)
	  break;
      if (j == eqmNR)
	fatal_error("No such method known",src->property[i].psource);
      dst->qm[j] = src->property[i].pvalue;
    }
  }
  
  for(i=0; (i<src->ncomposition); i++) {
    if (strcasecmp(src->composition[i].compname,"spoel") == 0) {
      for(k=0; (k<src->composition[i].ncatom); k++) {
	ptr = src->composition[i].catom[k].cname;	
	n   = src->composition[i].catom[k].cnumber;
	
	for(imiller=0; (imiller<NM2S); imiller++) {
	  if (strcasecmp(miller2spoel[imiller].miller,ptr) == 0) {
	    ptr = miller2spoel[imiller].spoel;
	    if (update_bm)
	      set_miller(dst,miller2spoel[imiller].miller,n,0);
	    printf("Converting %s to %s\n",miller2spoel[imiller].miller,
		   miller2spoel[imiller].spoel);
	    break;
	  }
	}
	
	for(j=0; (j<eatNR+eatExtra); j++)
	  if (strcasecmp(spoel[j].name,ptr) == 0) {
	    dst->frag_comp[j] += n;
	    if (update_bm) {
	      nh = spoel[j].nh;
	      if (imiller == NM2S)
		set_miller(dst,spoel[j].miller,n,nh);
	      set_bosque(dst,spoel[j].bosque,n,nh);
	    }
	    break;
	  }
      }
    }
    else if ((strcasecmp(src->composition[i].compname,"bosque") == 0) && !update_bm) {
      for(k=0; (k<src->composition[i].ncatom); k++) {
	for(j=0; (j<eelemNR); j++)
	  if (strcasecmp(bosque[j].name,src->composition[i].catom[k].cname) == 0) {
	    dst->elem_comp[j] += src->composition[i].catom[k].cnumber;
	    break;
	  }
      }
    }
    else if ((strcasecmp(src->composition[i].compname,"miller") == 0) && !update_bm) {
      for(k=0; (k<src->composition[i].ncatom); k++) {
	for(j=0; (j<emlNR); j++)
	  if (strcasecmp(miller[j].name,src->composition[i].catom[k].cname) == 0) {
	    dst->emil_comp[j] += src->composition[i].catom[k].cnumber;
	    break;
	  }
      }
    }
  }    
}

static void process_attr(FILE *fp,xmlAttrPtr attr,int elem,
			 int indent,t_xmlrec *xml)
  {
  char *attrname,*attrval;
  char buf[100];
  t_molprop2 *xptr=NULL;
  t_catom *ca;
  
  if (xml->npd > 0)
    xptr = xml->pd2[xml->npd-1];
  
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
    switch (elem) {
    case exmlMOLECULE: 
      if (atest("formula")) 
	xptr->formula = strdup(attrval);
      else if (atest("molname")) 
	xptr->molname = strdup(attrval);
      else if (atest("weight")) 
	xptr->weight = atof(attrval);
      else
	fatal_error("Unknown attribute for molecule",attrname);
      break;
      
    case exmlCATEGORY:
      if (atest("catname")) 
	xptr->category[xptr->ncategory-1] = strdup(attrval);
      else
	fatal_error("Unknown attribute for category",attrname);
      break;
      
    case exmlPROPERTY: {
      t_property *xp = &(xptr->property[xptr->nproperty-1]);
      if (atest("pname")) 
	xp->pname = strdup(attrval);
      else if (atest("psource")) 
	xp->psource = strdup(attrval);
      else if (atest("preference")) 
	xp->preference = strdup(attrval);
      else if (atest("pvalue")) 
	xp->pvalue = atof(attrval);
      else if (atest("perror")) 
	xp->perror = atof(attrval);
      else
	fatal_error("Unknown attribute for property",attrname);
      break;
    }
    case exmlCOMPOSITION: {
      if (xptr->ncomposition <= 0)
	fatal_error("Composition error","");
      t_composition *co = &(xptr->composition[xptr->ncomposition-1]);
      if (atest("compname")) 
	co->compname = strdup(attrval);
      else
	fatal_error("Unknown attribute for catom",attrname);
      break;
    }
    case exmlCATOM:
      if (xptr->ncomposition <= 0)
	fatal_error("Composition error","");
      ca = &(xptr->composition[xptr->ncomposition-1].catom[xptr->composition[xptr->ncomposition-1].ncatom-1]);
      if (atest("cname")) 
	ca->cname = strdup(attrval);
      else if (atest("cnumber")) 
	ca->cnumber = atoi(attrval);
      else
	fatal_error("Unknown attribute for catom",attrname);
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
}

static void process_element(FILE *fp,xmlNodePtr tree,int indent,t_xmlrec *xml)
{
  int elem;
  char buf[100];
  t_molprop2 *xptr=NULL;
  
  if (xml->npd > 0)
    xptr = xml->pd2[xml->npd-1];
  
  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  switch (elem) {
  case exmlMOLECULES:
    break;
  case exmlMOLECULE:
    xml->npd++;
    srenew(xml->pd2,xml->npd);
    snew((xml->pd2[xml->npd-1]),1);
    break;
  case exmlCATEGORY:
    xptr->ncategory++;
    srenew(xptr->category,xptr->ncategory);
    nullify(xptr->category[xptr->ncategory-1]);
    break;
  case exmlPROPERTY: 
    xptr->nproperty++;
    srenew(xptr->property,xptr->nproperty);
    nullify(xptr->property[xptr->nproperty-1]);
    break;
  case exmlCOMPOSITION:
    xptr->ncomposition++;
    srenew(xptr->composition,xptr->ncomposition);
    nullify(xptr->composition[xptr->ncomposition-1]);
    break;
  case exmlCATOM: {
    t_composition *co = &(xptr->composition[xptr->ncomposition-1]);
    
    co->ncatom++;
    srenew(co->catom,co->ncatom);
    nullify(co->catom[co->ncatom-1]);
    break;
  }
  default:
    break;
  }
  process_attr(fp,tree->properties,elem,indent+2,xml);
}
 
static void process_tree(FILE *fp,xmlNodePtr tree,int indent,t_xmlrec *xml)
{
  char buf[100];
  
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

static t_xmlrec *init_xml()
{
  t_xmlrec *xml;
  
  snew(xml,1);
  
  return xml;
}

static void done_xml(t_xmlrec **xml)
{
  free((*xml)->pd);
  free((*xml)->pd2);
  free(*xml);
  *xml = NULL;
}
	
int read_molprops(char *fn,t_molprop **pd,int update_bm)
{
  xmlDocPtr  doc;
  t_xmlrec   *xml;
  int        i,npd;
  
  xmlDoValidityCheckingDefaultValue = 0;
  if ((doc = xmlParseFile(fn)) == NULL) {
    fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.",
	    fn);
    exit(1);
  }

  xml = init_xml();
  process_tree(NULL,doc->children,0,xml);

  /* Copy xml to pd */
  npd = xml->npd;
  snew(*pd,npd);
  for(i=0; (i<npd); i++)
    copy_molprop(&((*pd)[i]),xml->pd2[i],update_bm);
      
  xmlFreeDoc(doc);
  done_xml(&xml);
  
  return npd;
}

static void add_xml_int(xmlNodePtr ptr,xmlChar *name,int val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%d",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    fatal_error("Setting",(char *)name);
}

static void add_xml_double(xmlNodePtr ptr,xmlChar *name,double val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%g",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    fatal_error("Setting",(char *)name);
}

static void add_xml_char(xmlNodePtr ptr,xmlChar *name,char *val)
{
  if (xmlSetProp(ptr,name,(xmlChar *)val) == 0)
    fatal_error("Setting",(char *)name);
}

static xmlNodePtr add_xml_child(xmlNodePtr parent,int type)
{
  xmlNodePtr child;
  
  if ((child = xmlNewChild(parent,NULL,(xmlChar *)exml_names[type],NULL)) == NULL)
    fatal_error("Creating element",(char *)exml_names[type]);
  
  return child;
}

static xmlNodePtr add_xml_comment(xmlDocPtr doc,
				  xmlNodePtr prev,xmlChar *comment)
{
  xmlNodePtr comm,ptr;
  
  if ((comm = xmlNewComment(comment)) == NULL)
    fatal_error("Creating doc comment element","");
  ptr = prev;
  while (ptr->next != NULL)
    ptr=ptr->next;
  ptr->next    = comm;
  comm->prev   = ptr;
  comm->doc    = doc;
  
  return comm;
}

static void add_xml_molprop(xmlNodePtr parent,t_molprop *mp)
{
  xmlNodePtr ptr,child,comp;
  int i;
  
  ptr = add_xml_child(parent,exmlMOLECULE);
  add_xml_char(ptr,(xmlChar *)exml_names[exmlMOLNAME],mp->molname);
  add_xml_char(ptr,(xmlChar *)exml_names[exmlFORMULA],mp->formula);
  add_xml_double(ptr,(xmlChar *)exml_names[exmlWEIGHT],mp->weight);
  
  for(i=0; (i<ecNR); i++)
    if (mp->category[i] > 0) {
      child = add_xml_child(ptr,exmlCATEGORY);
      add_xml_char(child,(xmlChar *)exml_names[exmlCATNAME],ec_name[i]);
    }
	
  for(i=0; (i<mp->nexperiment); i++) {
    child = add_xml_child(ptr,exmlPROPERTY);
    add_xml_char(child,(xmlChar *)exml_names[exmlPNAME],mp->pname[i]);
    add_xml_double(child,(xmlChar *)exml_names[exmlPVALUE],mp->experiment[i]);
    add_xml_char(child,(xmlChar *)exml_names[exmlPSOURCE],"Experiment");
    add_xml_char(child,(xmlChar *)exml_names[exmlPREFERENCE],mp->reference[i]);
  }
  for(i=0; (i<eqmNR); i++) {
    if (mp->qm[i] > 0) {
      child = add_xml_child(ptr,exmlPROPERTY);
      add_xml_char(child,(xmlChar *)exml_names[exmlPNAME],"Polarizability");
      add_xml_double(child,(xmlChar *)exml_names[exmlPVALUE],mp->qm[i]);
      add_xml_char(child,(xmlChar *)exml_names[exmlPSOURCE],lbasis[i]);
      add_xml_char(child,(xmlChar *)exml_names[exmlPREFERENCE],"Spoel2007a");
    }
  }
  child = add_xml_child(ptr,exmlCOMPOSITION);
  add_xml_char(child,(xmlChar *)exml_names[exmlCOMPNAME],"Miller");
  for(i=0; (i<emlNR); i++) 
    if (mp->emil_comp[i] > 0) {
      comp = add_xml_child(child,exmlCATOM);
      add_xml_char(comp,(xmlChar *)exml_names[exmlC_NAME],miller[i].name);
      add_xml_int(comp,(xmlChar *)exml_names[exmlC_NUMBER],mp->emil_comp[i]);
    }
  child = add_xml_child(ptr,exmlCOMPOSITION);
  add_xml_char(child,(xmlChar *)exml_names[exmlCOMPNAME],"Bosque");
  for(i=0; (i<eelemNR); i++) 
    if (mp->elem_comp[i] > 0) {
      comp = add_xml_child(child,exmlCATOM);
      add_xml_char(comp,(xmlChar *)exml_names[exmlC_NAME],bosque[i].name);
      add_xml_int(comp,(xmlChar *)exml_names[exmlC_NUMBER],mp->elem_comp[i]);
    }
  child = add_xml_child(ptr,exmlCOMPOSITION);
  add_xml_char(child,(xmlChar *)exml_names[exmlCOMPNAME],"Spoel");
  for(i=0; (i<eatNR+eatExtra); i++) 
    if (mp->frag_comp[i] > 0) {
      comp = add_xml_child(child,exmlCATOM);
      add_xml_char(comp,(xmlChar *)exml_names[exmlC_NAME],spoel[i].name);
      add_xml_int(comp,(xmlChar *)exml_names[exmlC_NUMBER],mp->frag_comp[i]);
    }
}

void write_molprops(char *fn,int npd,t_molprop pd[])
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
    fatal_error("Creating XML document","");
    
  if ((dtd = xmlCreateIntSubset(doc,dtdname,libdtdname,dtdname)) == NULL)
    fatal_error("Creating XML DTD","");
    
  if ((myroot = xmlNewDocNode(doc,NULL,gmx,NULL)) == NULL)
    fatal_error("Creating root element","");
  dtd->next = myroot;
  myroot->prev = (xmlNodePtr) dtd;
    
  /* Add molecule definitions */
  for(i=0; (i<npd); i++)
    add_xml_molprop(myroot,&(pd[i]));

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    fatal_error("Saving file",fn);
  xmlFreeDoc(doc);
}

void write_atoms_molprops(char *fn,char *name,char *formula,t_atoms*atoms,t_params *bonds)
{
  t_molprop *mp;
  int i,tp;
  
  snew(mp,1);
  mp->nr = atoms->nr;
  mp->formula = formula;
  mp->molname = name;
  mp->weight = 0;
  mp->nexperiment = 0;
  mp->experiment = NULL;
  mp->error = 0;
  mp->pname = NULL;
  mp->reference = NULL;
  for(i=0; (i<atoms->nr); i++) {
    tp = atoms->atom[i].type;
    range_check(tp,0,eatNR);
    mp->frag_comp[tp]++;
  }
  write_molprop(fn,1,mp);  
  sfree(mp);
}

#else

int read_molprops(char *fn,t_molprop **pd,int update_bm)
{
  gmx_fatal(FARGS,"You must compile this program with the libxml2 library");
}

void write_molprops(char *fn,int npd,t_molprop pd[])
{
  gmx_fatal(FARGS,"You must compile this program with the libxml2 library");
}

void write_atoms_molprops(char *fn,t_atoms*atoms,t_params *bonds)
{
  gmx_fatal(FARGS,"You must compile this program with the libxml2 library");
}

#endif

void fatal_error(char *str,char *val)
{
  fprintf(stderr,"Fatal: %s - %s\n",str,val);
  exit(1);
}

