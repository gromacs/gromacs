/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "atomprop.h"
#include "grompp.h"
#include "x2top_nm2type.h"
#include "x2top_eemprops.h"
#include "gmx_xml.h"
#include "gmx_elements.h"
#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>

extern int xmlDoValidityCheckingDefaultValue;
	
typedef struct {
  char *type,*ref;
  real J0,w,chi0;
} t_eemrec; 

typedef struct {
  int nbonds;
  char **elem;
  real *blen;
} t_bondrec;

typedef struct {
  char *type,*ref;
  real alpha;
  int  nh;
} t_polrec;

typedef struct {
  char *type,*ref;
  int  nparam;
  real *param;
} t_ffrec;

typedef struct {
  char      *sname,*lname;
  int       atomnumber,maxbonds;
  real      mass,bondlength;
  int       neemrec;
  t_eemrec  *eemrec;
  int       nbondrec;
  t_bondrec *bondrec;
  int       npolrec;
  t_polrec  *polrec;
  int       nffrec;
  t_ffrec   *ffrec;
} t_elemrec;

typedef struct {
  int       nelem;
  t_elemrec *er;
} gmx_elements_t;

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
  exmlELEMENT, exmlELEM_SNAME, exmlELEM_LNAME,
  exmlATOMNUMBER, exmlMASS, exmlBONDLENGTH,
  exmlEEMREC, exmlEEMTYPE, exmlEEMREF, 
  exmlEEMJ0, exmlEEMWEIGHT, exmlEEMCHI0,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "element", "shortname", "longname",
  "atomicnumber", "mass", "bondlength",
  "electronegativity", "eem_type", "eem_ref", 
  "J0", "weight", "chi0"
};

static void fill_elemrec(void *atomprop,
			 void *eem,
			 x2top_nm2t nm2t,
			 t_elemrec *elem,
			 int atomnumber,
			 char *eshort,
			 char *elong)
{
  int  j,nn;
  real wj,val;
  
  elem->sname = strdup(eshort);
  elem->lname = strdup(elong);
  elem->atomnumber = atomnumber;
  if (query_atomprop(atomprop,epropMass,"???",eshort,&val))
    elem->mass = val;
  for(j=0; (j<eqgNR); j++) 
    if ((nn = eem_get_index(eem,atomnumber,j)) != -1) {
      elem->neemrec++;
      srenew(elem->eemrec,elem->neemrec);
      elem->eemrec[elem->neemrec-1].type = strdup(get_eemtype_name(j));
      elem->eemrec[elem->neemrec-1].ref  = strdup(get_eemtype_reference(j));
      elem->eemrec[elem->neemrec-1].J0   = lo_get_j00(eem,nn,&wj,0);
      elem->eemrec[elem->neemrec-1].w    = wj;
      elem->eemrec[elem->neemrec-1].chi0 = eem_get_chi0(eem,nn);
    }
}

gmx_elements gather_element_information(void *atomprop,
					x2top_nm2t nm2t)
{
  gmx_elements_t *elements;
  int            i,j,k,nstr;
  real           val;
  char           **str=NULL;
  char           eshort[32],elong[64];
  char           *db = "elemnew.dat";
  void           *eem;
					
  eem = read_eemprops(NULL,-1,atomprop);
  nstr = get_file(db,&str);
  snew(elements,1);
  elements->nelem = nstr;
  snew(elements->er,elements->nelem);
  for(i=0; (i<nstr); i++) {
    if (sscanf(str[i],"%s%s%d",eshort,elong,&k) == 3) {
      range_check(k,1,nstr+1);
      fill_elemrec(atomprop,eem,nm2t,&(elements->er[k-1]),
		   k,eshort,elong);
    }
    else
      gmx_fatal(FARGS,"Error in file %s line %d",db,i);
    sfree(str[i]);
  }
  sfree(str);
  
  return (gmx_elements) elements;
}

static void add_xml_eemrecs(xmlNodePtr mother,t_elemrec *elem)
{
  int i;
  xmlNodePtr self;
  
  for(i=0; (i<elem->neemrec); i++) {
    self = add_xml_child(mother,(xmlChar *)exml_names[exmlEEMREC]);
    add_xml_char(self,(xmlChar *)exml_names[exmlEEMTYPE],
		 elem->eemrec[i].type);
    add_xml_char(self,(xmlChar *)exml_names[exmlEEMREF],
		 elem->eemrec[i].ref);
    add_xml_double(self,(xmlChar *)exml_names[exmlEEMJ0],
		   elem->eemrec[i].J0);
    add_xml_double(self,(xmlChar *)exml_names[exmlEEMWEIGHT],
		   elem->eemrec[i].w);
    add_xml_double(self,(xmlChar *)exml_names[exmlEEMCHI0],
		   elem->eemrec[i].chi0);
  }
}

static void add_xml_element(xmlNodePtr myroot,t_elemrec *elem)
{
  xmlNodePtr self = add_xml_child(myroot,exml_names[exmlELEMENT]);
  add_xml_char(self,(xmlChar *)exml_names[exmlELEM_SNAME],elem->sname);
  add_xml_char(self,(xmlChar *)exml_names[exmlELEM_LNAME],elem->lname);
  add_xml_int(self,(xmlChar *)exml_names[exmlATOMNUMBER],elem->atomnumber);
  add_xml_double(self,(xmlChar *)exml_names[exmlMASS],elem->mass);
  add_xml_double(self,(xmlChar *)exml_names[exmlBONDLENGTH],elem->bondlength);
  add_xml_eemrecs(self,elem);
}

void write_gmx_elements(char *fn,gmx_elements elem,bool bCompress)
{
  Gmx_elements_t *elements = (gmx_elements_t *) elem;
  xmlDocPtr  doc;
  xmlDtdPtr  dtd;
  xmlNodePtr myroot;
  int        i,nmt;
  xmlChar    *libdtdname,*dtdname,*gmx;
  
  gmx        = (xmlChar *) "chem_elements";
  dtdname    = (xmlChar *) "chem_elements.dtd";
  libdtdname = dtdname;
  
  if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
    gmx_fatal(FARGS,"Creating XML document");
    
  if ((dtd = xmlCreateIntSubset(doc,dtdname,libdtdname,dtdname)) == NULL)
    gmx_fatal(FARGS,"Creating XML DTD");
    
  if ((myroot = xmlNewDocNode(doc,NULL,gmx,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating root element");
  dtd->next = myroot;
  myroot->prev = (xmlNodePtr) dtd;
    
  /* Add molecule definitions */
  for(i=0; (i<elements->nelem); i++)
    add_xml_element(myroot,&(elements->er[i]));

  xmlSetDocCompressMode(doc,bCompress ? 1 : 0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving file %s",fn);
  xmlFreeDoc(doc);
}

#else 

void write_gmx_elements(char *fn,gmx_elements elem,bool bCompress)
{
  ;
}
gmx_elements gather_element_information(void *atomprop,
					x2top_nm2t nm2t)
{
  ;
} 
#endif
