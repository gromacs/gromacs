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
#include "poldata_xml.h"
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
  exmlGENTOP,
  exmlSMATOMS, exmlPOLAR_UNIT, exmlBLENGTH_UNIT,
  exmlSMATOM, exmlELEM, 
  exmlSMNAME, exmlMILLER_EQUIV, exmlNHYDROGEN, exmlCHARGE,
  exmlHYBRIDIZATION, exmlPOLARIZABILITY, exmlBLENGTH, 
  exmlSMBONDS, exmlSMBOND,
  exmlATOM1, exmlATOM2, 
  exmlBSATOMS, exmlBSATOM,
  exmlMILATOMS, exmlTAU_UNIT, exmlAHP_UNIT,
  exmlMILATOM, exmlMILNAME, exmlSPOEL_EQUIV,
  exmlATOMNUMBER, exmlTAU_AHC, exmlALPHA_AHP,
  exmlSYMMETRIC_CHARGES, exmlSYM_CHARGE,
  exmlBASE, exmlATOM3,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "gentop",
  "smatoms", "polarizability_unit", "blength_unit",
  "smatom", "elem", 
  "smname", "miller_equiv", "nhydrogen", "charge",
  "hybridization", "polarizability", "blength", 
  "smbonds", "smbond",
  "atom1", "atom2", 
  "bsatoms", "bsatom",
  "milatoms", "tau_ahc_unit", "alpha_ahp_unit",
  "milatom", "milname", "spoel_equiv",
  "atomnumber", "tau_ahc", "alpha_ahp",
  "symmetric_charges", "sym_charge",
  "base", "atom3"
};

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
			  int indent,gmx_poldata_t pd)
{
  char *attrname,*attrval;
  char buf[100];
  int  i,kkk;
  char *xbuf[exmlNR];
  
  for(i=0; (i<exmlNR); i++)
    xbuf[i] = NULL;    
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
    kkk = find_elem(attrname,exmlNR,exml_names);
    if (attrval != NULL)
      xbuf[kkk] = strdup(attrval);
      
    if (fp)
      fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent,buf,99),
	      attrname,attrval);
    attr = attr->next;
#undef atest
  }
  /* Done processing attributes for this element. Let's see if we still need
   *  to interpret them.
   */
  if (xbuf[exmlPOLAR_UNIT] && xbuf[exmlBLENGTH_UNIT]) 
    gmx_poldata_set_spoel_units(pd,xbuf[exmlPOLAR_UNIT],xbuf[exmlBLENGTH_UNIT]);
  else if (xbuf[exmlPOLAR_UNIT]) 
    gmx_poldata_set_bosque_units(pd,xbuf[exmlPOLAR_UNIT]);
  else if (xbuf[exmlTAU_UNIT] && xbuf[exmlAHP_UNIT]) 
    gmx_poldata_set_miller_units(pd,xbuf[exmlTAU_UNIT],xbuf[exmlAHP_UNIT]);
  else if (xbuf[exmlELEM] && xbuf[exmlMILLER_EQUIV] && xbuf[exmlNHYDROGEN] &&
	   xbuf[exmlCHARGE] && xbuf[exmlHYBRIDIZATION] && 
	   xbuf[exmlPOLARIZABILITY] && xbuf[exmlBLENGTH])
    gmx_poldata_add_spoel(pd,xbuf[exmlELEM],xbuf[exmlMILLER_EQUIV],
			  atoi(xbuf[exmlNHYDROGEN]),
			  atoi(xbuf[exmlCHARGE]),
			  atoi(xbuf[exmlHYBRIDIZATION]),
			  atof(xbuf[exmlPOLARIZABILITY]),
			  atof(xbuf[exmlBLENGTH]));
  else if (xbuf[exmlMILNAME] && xbuf[exmlATOMNUMBER] && xbuf[exmlTAU_AHC] &&
	   xbuf[exmlALPHA_AHP]) 
    gmx_poldata_add_miller(pd,xbuf[exmlMILNAME],atoi(xbuf[exmlATOMNUMBER]),
			   atof(xbuf[exmlTAU_AHC]),atof(xbuf[exmlALPHA_AHP]),xbuf[exmlSPOEL_EQUIV]);
  else if (xbuf[exmlELEM] && xbuf[exmlPOLARIZABILITY]) 
    gmx_poldata_add_bosque(pd,xbuf[exmlELEM],atof(xbuf[exmlPOLARIZABILITY]));
  else if (debug) {
    fprintf(debug,"Unknown combination of attributes:\n");
    for(i=0; (i<exmlNR); i++)
      if (xbuf[i] != NULL)
	fprintf(debug,"%s = %s\n",exml_names[i],xbuf[i]);
  }
    
  /* Clean up */
  for(i=0; (i<exmlNR); i++)
    if (xbuf[i] != NULL)
      sfree(xbuf[i]);

}

static void process_element(FILE *fp,xmlNodePtr tree,int indent,gmx_poldata_t pd)
{
  int elem;
  char buf[100];
  
  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  switch (elem) {

  default:
    break;
  }
  if (elem != exmlGENTOP)
    process_attr(fp,tree->properties,elem,indent+2,pd);
}
 
static void process_tree(FILE *fp,xmlNodePtr tree,int indent,gmx_poldata_t pd)
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
      process_element(fp,tree,indent+2,pd);
      
      if (tree->children)
	process_tree(fp,tree->children,indent+2,pd);
      break;
    default:
      break;
    }
    tree = tree->next;
  }
}

gmx_poldata_t gmx_poldata_read(char *fn)
{
  xmlDocPtr     doc;
  int           i,npd;
  gmx_poldata_t pd;
  
  xmlDoValidityCheckingDefaultValue = 0;
  if ((doc = xmlParseFile(fn)) == NULL) {
    fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.",
	    fn);
    exit(1);
  }

  pd = gmx_poldata_init();
  process_tree(debug,doc->children,0,pd);

  xmlFreeDoc(doc);
  
  return pd;
}

static void add_xml_poldata(xmlNodePtr parent,gmx_poldata_t pd)
{
  xmlNodePtr child,grandchild,comp;
  int    i,nhydrogen,charge,hybridization,atomnumber;
  char *elem,*miller_equiv,*name,*spoel_equiv;
  double polarizability,blength,tau_ahc,alpha_ahp;
  
  child = add_xml_child(parent,exml_names[exmlSMATOMS]);
  add_xml_char(child,exml_names[exmlPOLAR_UNIT],"A3");
  add_xml_char(child,exml_names[exmlBLENGTH_UNIT],"nm");
  
  while ((name = gmx_poldata_get_spoel(pd,NULL,&elem,&miller_equiv,
				       &nhydrogen,&charge,&hybridization,
				       &polarizability,&blength)) != NULL) {
    grandchild = add_xml_child(child,exml_names[exmlSMATOM]);
    add_xml_char(grandchild,exml_names[exmlELEM],elem);
    add_xml_char(grandchild,exml_names[exmlSMNAME],name);
    add_xml_char(grandchild,exml_names[exmlMILLER_EQUIV],miller_equiv);
    add_xml_int(grandchild,exml_names[exmlNHYDROGEN],nhydrogen);
    add_xml_int(grandchild,exml_names[exmlCHARGE],charge);
    add_xml_int(grandchild,exml_names[exmlHYBRIDIZATION],hybridization);
    add_xml_double(grandchild,exml_names[exmlPOLARIZABILITY],polarizability);
    add_xml_double(grandchild,exml_names[exmlBLENGTH],blength);
  }

  child = add_xml_child(parent,exml_names[exmlBSATOMS]);
  add_xml_char(child,exml_names[exmlPOLAR_UNIT],"A3");
  
  while ((name = gmx_poldata_get_bosque(pd,NULL,&polarizability)) != NULL) {
    grandchild = add_xml_child(child,exml_names[exmlBSATOM]);
    add_xml_char(grandchild,exml_names[exmlELEM],name);
    add_xml_double(grandchild,exml_names[exmlPOLARIZABILITY],polarizability);
  }
  child = add_xml_child(parent,exml_names[exmlMILATOMS]);
  add_xml_char(child,exml_names[exmlTAU_UNIT],"A3/2");
  add_xml_char(child,exml_names[exmlAHP_UNIT],"A3");
  
  while ((name = gmx_poldata_get_miller(pd,NULL,&atomnumber,&tau_ahc,&alpha_ahp,&spoel_equiv)) != NULL) {
    grandchild = add_xml_child(child,exml_names[exmlMILATOM]);
    add_xml_char(grandchild,exml_names[exmlMILNAME],name);
    add_xml_int(grandchild,exml_names[exmlATOMNUMBER],atomnumber);
    add_xml_double(grandchild,exml_names[exmlTAU_AHC],tau_ahc);
    add_xml_double(grandchild,exml_names[exmlALPHA_AHP],alpha_ahp);
    if (spoel_equiv)
      add_xml_char(grandchild,exml_names[exmlSPOEL_EQUIV],spoel_equiv);
  }
  
}

void gmx_poldata_write(char *fn,gmx_poldata_t pd)
{
  xmlDocPtr  doc;
  xmlDtdPtr  dtd;
  xmlNodePtr myroot;
  int        i,nmt;
  xmlChar    *libdtdname,*dtdname,*gmx;
  
  gmx        = (xmlChar *) "gentop";
  dtdname    = (xmlChar *) "gentop.dtd";
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
  add_xml_poldata(myroot,pd);

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving file",fn);
  xmlFreeDoc(doc);
}



