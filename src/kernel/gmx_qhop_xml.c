/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata.c,v 1.20 2009/05/17 13:56:55 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
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
#include <string.h>
#include "gmx_fatal.h"
#include "smalloc.h"
#include "macros.h"
#include "futil.h"
#include "string2.h"
#include "types/gmx_qhop_types.h"
#include "gmx_qhop_parm.h"
#include "gmx_qhop_db.h"
#include "gmx_xml.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
	
typedef struct xmlrec *xmlrec_t;

typedef struct xmlrec {
  int        nqh;
  qhop       **gqh;         /* Hopping parameters */
  qhop_resblocks_t    rb;   /* Bunches of related residues
			     * and their interaction parameters */
  t_symtab   tab;
} xmlrec;


extern int xmlDoValidityCheckingDefaultValue;
	
#define NN(x) (NULL != (x))
#define CPYLEN 128 /* This is the default length used for strncpy.
		    * I avoid using strcpy becase of potential buffer-overflows
		    * Note that the careless use of CPYLEN isn't safe either, but
		    * as the code maturates, CPYLEN may be replaced by explicit numbers.
		    * - Erik Marklund */
#define RNMLEN 6 /* Ditto for resnames. This can probabaly be a bit shorter,
		  * but let's keep it for now. */
#define ANMLEN 6 /* Ditto for atomnames. */

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
  exmlTITRATIONS,
  exmlRESBLOCKS,
  exmlRESIDUE_TYPE, exmlFILE, exmlRES, exmlPROTON, exmlPRODUCT,
  exmlBWATER, exmlTITRATION, exmlDONOR, exmlACCEPTOR,
  exmlPARAM, exmlNAME, exmlPROTONATED, exmlVALUE, exmlDESCRIPTION,
  exmlUNIT, exmlDON_ATOM, exmlACC_ATOM,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "titrations",
  "resblocks",
  "residue_type", "file", "res", "proton", "product",
  "water", "titration", "donor", "acceptor", "parameter", 
  "name", "protonated", "value", "description", "unit", "don_atom", "acc_atom"
};

static char *strip_spaces(char *s)
{
  int i, j;

  i = 0;
  while (s[i] != '\0')
    {
      if (s[i] != ' ')
	{
	  break;
	}
      i++;
    }

  if (s[i] == '\0')
    {
      return s;
    }

  j = i;
  while (s[j] != '\0')
    {
      if (s[j] == ' ')
	{
	  s[j] = '\0';
	  return &(s[i]);
	}
      j++;
    }

  gmx_fatal(FARGS, "String could not be stripped of spaces.");
  
  return (char *)NULL;
}

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

/* Finds the first matching string in an array of strings. Returns -1 if not found. */
static int qhop_find_name(char **ss, const char *name, int n)
{
  int i;
  if (ss != NULL)
    {
      for (i=0; i<n; i++)
	{
	  if (strcmp(ss[i], name) == 0)
	    {
	      return i;
	    }
	}
    }

  return -1;/* not found */
}

/* Adds a file with name *f to a qhop_resblocks*/
static void rb_add_file(qhop_resblocks *rb, char *f, t_symtab *tab)
{
  if NN(f)
    {
      if (qhop_find_name(rb->files, f, rb->nf) < 0) /* New file? */
	{
	  printf("Adding file %s\n", f);
	  /* Copy the filename to the array in rb.*/
	  srenew(rb->files, rb->nf+1);
	  rb->files[rb->nf++] = *(put_symtab(tab, f));
	}
    }
  else
    gmx_fatal(FARGS, "No filename found.");
}

static void rb_add_restype(qhop_resblocks_t rb, currentRes *ri, 
			   char *description,char *canonical, 
			   char *protonated,char *water, t_symtab *tab)
{
  char *s;
  int  i;
  
  if NN(canonical)
    {
      if (NULL != debug)
	fprintf(debug,"Adding restype %s (%s)\n", canonical,description);
      /* Is this a new resblock? If not, make a new one. */
      for(i=0; (i<rb->nrestypes); i++) {
	if (strcmp(rb->qrt[i].canonical,canonical) == 0)
	  break;
      }
      if (i == rb->nrestypes) 
	{
	  srenew(rb->qrt, rb->nrestypes+1);       /* Make new resblock */
	  memset(&rb->qrt[rb->nrestypes],0,sizeof(rb->qrt[rb->nrestypes]));
	  
	  ri->rt = rb->nrestypes;            /* Make this resblock the current one */
	  ri->r  = -1;                       /* just to be sure... */
	  ri->da = -1;
	  rb->qrt[rb->nrestypes].canonical = strdup(canonical); 
	  rb->qrt[rb->nrestypes].protonated = strdup(protonated); 
	  rb->qrt[rb->nrestypes].description = strdup(description); 
	  rb->qrt[rb->nrestypes].bWater = (strcasecmp(water, "TRUE") == 0);
	  rb->nrestypes++;
	}
    }
  else
    gmx_fatal(FARGS, "No name found for this resblock.");
}

static void rb_add_res(qhop_resblocks_t rb, currentRes *ri, char *name, t_symtab *tab)
{
  gmx_bool found;
  int i;
  if NN(name)
    {
      if (NULL != debug)
	fprintf(debug,"Adding res %s\n", name);
      /* Have we seen this res before for this resblock? */
      found = FALSE;
      for (i=0; i < rb->qrt[ri->rt].nsubres && !found; i++)
	if (strcmp(rb->qrt[ri->rt].subres[i].name, name) == 0)
	  found = TRUE;

      if (!found)
	{
	  ri->r = rb->qrt[ri->rt].nsubres;
	  if (ri->r == 0)
	    snew(rb->qrt[ri->rt].subres, 1);
	  else
	    {	  
	      srenew(rb->qrt[ri->rt].subres, rb->qrt[ri->rt].nsubres+1);
	      memset(&(rb->qrt[ri->rt].subres[ri->r]),0,sizeof(rb->qrt[ri->rt].subres[ri->r]));
	    }
	  rb->qrt[ri->rt].subres[rb->qrt[ri->rt].nsubres++].name = strdup(name);
	  rb->qrt[ri->rt].subres[ri->r].irtp = NOTSET;
	  /*add_to_record(name, &(rb->res[ri->rt][ri->r].name),
	    &(rb->nres[ri->rt]), RNMLEN, eUPDATE_INDEX);*/
	}
      else
	{
	  ri->r = i;
	}
    }
  else
    gmx_fatal(FARGS, "No name found for this residue.");

}

static void rb_add_proton(qhop_subres **res, currentRes *ri, char *name,
			  t_symtab *tab)
{
  qhop_reactant *r;
  int nr;
 
  r = (ri->Acc == exmlACCEPTOR) ?
    &(res[ri->rt][ri->r].acc[ri->da]) :
    &(res[ri->rt][ri->r].don[ri->da]);
      
  if NN(name)
    {
      printf("Adding proton %s\n", name);

      if (qhop_find_name(r->H, name, r->nH) < 0) /* New proton? */
	{
	  if (r->nH != 0)
	    srenew(r->H, r->nH + 1);
	  else
	    snew(r->H, 1);

	  r->H[r->nH++] = *(put_symtab(tab, name)); /*trim_strndup(name,6);*/
	  /*add_to_record(name,
	    rb->res[ri->rt][ri->r].acc[ri->da].H,
	    &(rb->res[ri->rt][ri->r].acc[ri->da].nH),
	    ANMLEN, eUPDATE_INDEX);*/
	}
    }
}

/* adds a qhop_reactant with name and product to a qhop_res 'res'.
 * If acc_or_don == exmlACCEPTOR then this is an acceptor, otherwise a donor. */
static void rb_add_name_product(qhop_subres_t res, const char *name,
				const char *product, currentRes *ri, 
				t_symtab *tab)
{
  /* Now, the argument 'name' may be a list of atoms,
   * because of proton tautomerism. */
  int i, j, k, nda, nnames;
  gmx_bool acc, w;
  char **s, *buf;
  qhop_reactant *da;
  
  if (NN(name) && NN(product))
    {
      if (NULL != debug)
	fprintf(debug,"Adding name(s) %s and product %s \n", name, product);

      /* Split the name line. */
      buf=strdup(name);
      j = strlen(buf);
      nnames = 0;
      snew(s,1);
      w = FALSE; /* Flags whether we're in a word (TRUE) or a whitespace region (FALSE) */
      for (i=0; i<j; i++)
	{
	  if (buf[i] == ' '  ||
	      buf[i] == '\t' ||
	      buf[i] == '\n') /* Allowed separators */
	    {
	      w = FALSE; /* We're not within a word. */
	      buf[i] = '\0'; /* So we don't have to save the lengths or end positions of the words. */
	    }
	  else
	    {
	      if (w == FALSE) /* We've hit the start of another word! */
		{		  
		  srenew(s,nnames+1);
		  s[nnames++] = &(buf[i]);
		  w = TRUE;
		}
	    }
	}

      if (ri->Acc == exmlACCEPTOR)
	{
	  nda = res->na;
	  acc = TRUE;
	  da = res->acc;
	}
      else
	{
	  nda = res->nd;
	  acc = FALSE;
	  da = res->don;
	}

      /* Have we seen it before?*/
      for(i=0; i < nda; i++)
	for (j=0; j < da->nname; j++) /* Loop over da[i].name[j] */
	  for (k=0; k < nnames; k++)  /* Loop over name[k] */
	    if (strcmp(da[i].name[j], name) == 0)
	      gmx_fatal(FARGS, "%s %s has already been encountered.",
			acc ? "Acceptor":"Donor", s[k]);

      if (nda != 0)
	{
	  srenew(da, nda+1);
	  memset(&(da[nda]), 0, sizeof(qhop_reactant));
	}
      else
	{
	  snew(da, 1);
	}

      srenew(da[nda].name, nnames);
      da[nda].nname = nnames;
      for (i=0; i<nnames; i++)
	da[nda].name[i] = *(put_symtab(tab, s[i]));
      da[nda].product = *(put_symtab(tab,product));
      da[nda].nH = 0;
      ri->da = (nda)++;
      if(acc)
	{
	  res->na = nda;
	  res->acc = da;
	}
      else
	{
	  res->nd = nda;
	  res->don = da;
	}
      sfree(s);
      sfree(buf);
    }
  else
    gmx_fatal(FARGS, "%s without %s%s%s. Can't have that!",
	      ri->Acc==exmlACCEPTOR   ? "Acceptor" : "Donor",
	      NN(name)                ? "": "name",
	      NN(name) || NN(product) ? "":" and ",
	      NN(product)             ? "": "product");
}

/* rt: residue type/resblock
 * r:  res in resblock. */
static void qhop_process_attr(FILE *fp,xmlAttrPtr attr,int parent,
			      int elem,int indent, xmlrec_t xml, currentRes *ri)
{
  char *attrname,*attrval;
  char buf[100];
  int  i,kkk,eprop;
  gmx_bool found;
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
  
  switch (elem) {
  case exmlTITRATION:
    if (NN(xbuf[exmlDONOR]) && NN(xbuf[exmlACCEPTOR])) {
      if (NULL != debug)
	fprintf(debug,"Adding qhop: ");
      qhop_set_donor(xml->gqh[xml->nqh-1],xbuf[exmlDONOR]);
      qhop_set_acceptor(xml->gqh[xml->nqh-1],xbuf[exmlACCEPTOR]);
    }
    break;
  case exmlPARAM:
    if (NN(xbuf[exmlNAME]) && NN(xbuf[exmlUNIT]) && NN(xbuf[exmlVALUE])) {
      qhop_add_param(xml->gqh[xml->nqh-1],xbuf[exmlNAME],xbuf[exmlVALUE],
		     xbuf[exmlUNIT]);
    }
    break;
    /* Here's where some resblocks stuff needs to be implemented. */
  case exmlRESIDUE_TYPE: /* 'name' */
    rb_add_restype(xml->rb, ri, xbuf[exmlDESCRIPTION], xbuf[exmlNAME], 
		   xbuf[exmlPROTONATED],xbuf[exmlBWATER], &(xml->tab));
    break;
  case exmlRES: /* 'name' */
    rb_add_res(xml->rb, ri, xbuf[exmlNAME],&(xml->tab));
    break;

  case exmlACCEPTOR: /* 'name', 'product' */
    rb_add_name_product(&(xml->rb->qrt[ri->rt].subres[ri->r]),
			xbuf[exmlNAME], xbuf[exmlPRODUCT], ri, &(xml->tab));
    break;
  case exmlDONOR: /* 'name', 'product' */
    rb_add_name_product(&(xml->rb->qrt[ri->rt].subres[ri->r]),
			xbuf[exmlNAME], xbuf[exmlPRODUCT], ri, &(xml->tab));
    break;

  default:
    break;
  }
  

  for(i=0; (i<exmlNR); i++)
    if (NN(xbuf[i]))
      sfree(xbuf[i]);
}

static void qhop_process_element(FILE *fp, xmlNodePtr tree, int parent,
				 int indent, xmlrec_t xml, currentRes *ri)
{
  int elem, i;
  char buf[100], tmp[100],*ptr;

  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  switch(elem)
    {
    case exmlTITRATION:
      xml->nqh++;
      srenew(xml->gqh,xml->nqh);
      xml->gqh[xml->nqh-1] = qhop_init();
      break;

    case exmlRESBLOCKS:
      /* Nothing. */
      break;
    case exmlRESIDUE_TYPE:
      /* A new resblock will be created if this is indeed a new one.
       * That all happens in qhop_process_attr below. *rt will also be set accordingly.
       * Thus, nothing need to be done here.*/
      break;

    case exmlFILE:
      rb_add_file(xml->rb, (char*) (tree->children->content), &(xml->tab));
      break;

    case exmlRES:
      /* The identity of the resblock was determined at
       * the last encounter of an exmlRESIDUE_TYPE and the
       * subsequent call to qhop_process_attr.
       * r will be set in the next call to qhop_process_attr.
       * Thus, nothing to be done here */
      break;

      /* ri->da will be set when any of the following two are encountered. */
    case exmlACCEPTOR:
      /* A new acceptor will be added in the next call
       * to qhop_process_attr. */
      ri->Acc = exmlACCEPTOR;
      break;

    case exmlDONOR:
      /* A new donor will be added in the next call
       * to qhop_process_attr. */
      ri->Acc = exmlDONOR;
      break;

    case exmlPROTON:
      /* Now add this proton to the current acceptor/donor. */
      /* Explicit protons are not in the xml-file anymore.
      rb_add_proton(xml->rb->res, ri, (char*)(tree->children->content));
      */
      break;

    case exmlDON_ATOM:
      /* This is a child of qhop. */
      sprintf(tmp, "%s", (char*)(tree->children->content));
      while (NULL != (ptr = strchr(tmp,'\n'))) {
	*ptr = ' ';
      }
      trim(tmp);
      qhop_set_don_atom(xml->gqh[xml->nqh-1],tmp);
      break;

    case exmlACC_ATOM:
      /* This is a child of qhop. */
      sprintf(tmp, "%s", (char*)(tree->children->content));
      while (NULL != (ptr = strchr(tmp,'\n'))) {
	*ptr = ' ';
      }
      trim(tmp);
      qhop_set_acc_atom(xml->gqh[xml->nqh-1],tmp);
      break;

    default:
      break;
    }
  /*  if (elem != exmlTITRATIONS)*/
  if NN(tree->properties)
    qhop_process_attr(fp,tree->properties,parent,
		      elem,indent+2,/*xml->gqh[xml->nqh-1],*/
		      xml, ri);
}
 
static void qhop_process_tree(FILE *fp, xmlNodePtr tree, int parent,
			      int indent, xmlrec_t xml, currentRes *ri)
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
      qhop_process_element(fp, tree, parent, indent+2, xml, ri);
      
      if (tree->children) {
	elem = find_elem((char *)tree->name,exmlNR,exml_names);
	qhop_process_tree(fp,tree->children,elem,indent+2,xml,ri);
      }
      break;
    default:
      break;
    }
    tree = tree->next;
  }
}

qhop_db_t qhops_read(char *fn)
{
  xmlDocPtr     doc;
  int           i,npd;
  xmlrec        *xml;
  const char    *db="qhops.dat";
  gmx_bool      fna=FALSE;
  currentRes    ri;
  qhop_db_t     qdb;
  
  ri.Acc = FALSE;
  ri.rt  = -1;
  ri.r   = -1;
  ri.da  = -1;

  xmlDoValidityCheckingDefaultValue = 0;
  if (NULL == fn) 
    {
      fn = (char *)gmxlibfn(db);
      fna=TRUE;
    }
  if ((doc = xmlParseFile(fn)) == NULL)
    {
      fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.", fn);
      exit(1);
    }

  snew(qdb,1);
  open_symtab(&(qdb->tab));
  snew(xml,1);
  xml->rb = &(qdb->rb);
  xml->rb->nf = 0;
  xml->tab = qdb->tab;
  qhop_process_tree(NULL,doc->children,0,0,xml,&ri);
  
  xmlFreeDoc(doc);

  qdb->ngqh = xml->nqh;
  qdb->gqh = xml->gqh;
  
  return qdb;
}

static void add_xml_qhop(xmlNodePtr parent,qhop_t qht)
{
  xmlNodePtr ptr,child,grandchild,comp;
  char   *name,*type,*value,*unit;
  
  ptr = add_xml_child(parent,(char *)exml_names[exmlTITRATION]);
  add_xml_char(ptr,exml_names[exmlDONOR],qhop_get_donor(qht));
  add_xml_char(ptr,exml_names[exmlACCEPTOR],qhop_get_acceptor(qht));
  child = add_xml_child_val(ptr,exml_names[exmlDON_ATOM],
			    qhop_get_don_atom(qht));
  child = add_xml_child_val(ptr,exml_names[exmlACC_ATOM],
			    qhop_get_acc_atom(qht));
  while (qhop_get_param(qht,&name,&value,&unit) == 1) {
    child = add_xml_child(ptr,exml_names[exmlPARAM]);
    add_xml_char(child,exml_names[exmlNAME],name);
    add_xml_char(child,exml_names[exmlVALUE],value);
    add_xml_char(child,exml_names[exmlUNIT],unit);
    sfree(name);
    sfree(value);
    sfree(unit);
  }
}

static void add_xml_resblocks(xmlNodePtr parent,qhop_resblocks_t qdb)
{
  xmlNodePtr rb,res,subres,ad;
  int        i,j,k,l;
  
  rb = add_xml_child(parent,exml_names[exmlRESBLOCKS]);
  for(i=0; (i<qdb->nrestypes); i++) 
    {
      res = add_xml_child(rb,exml_names[exmlRESIDUE_TYPE]);
      add_xml_char(res,exml_names[exmlDESCRIPTION],qdb->qrt[i].description);
      add_xml_char(res,exml_names[exmlNAME],qdb->qrt[i].canonical);
      add_xml_char(res,exml_names[exmlPROTONATED],qdb->qrt[i].protonated);
      if (qdb->qrt[i].bWater)
	add_xml_char(res,exml_names[exmlBWATER],"TRUE");
      else
	add_xml_char(res,exml_names[exmlBWATER],"FALSE");
      for(j=0; (j<qdb->qrt[i].nsubres); j++) 
	{
	  subres = add_xml_child(res,exml_names[exmlRES]);
	  add_xml_char(subres,exml_names[exmlNAME],qdb->qrt[i].subres[j].name);
	  for(k=0; (k<qdb->qrt[i].subres[j].na); k++) 
	    {
	      for(l=0; (l<qdb->qrt[i].subres[j].acc[k].nname); l++) 
		{
		  ad = add_xml_child(subres,exml_names[exmlACCEPTOR]);
		  add_xml_char(ad,exml_names[exmlNAME],qdb->qrt[i].subres[j].acc[k].name[l]);
		  add_xml_char(ad,exml_names[exmlPRODUCT],qdb->qrt[i].subres[j].acc[k].product);
		}
	    }
	  for(k=0; (k<qdb->qrt[i].subres[j].nd); k++) 
	    {
	      for(l=0; (l<qdb->qrt[i].subres[j].don[k].nname); l++) 
		{
		  ad = add_xml_child(subres,exml_names[exmlDONOR]);
		  add_xml_char(ad,exml_names[exmlNAME],qdb->qrt[i].subres[j].don[k].name[l]);
		  add_xml_char(ad,exml_names[exmlPRODUCT],qdb->qrt[i].subres[j].don[k].product);
		}
	    }
	}
    }
}

void qhops_write(char *fn,qhop_db_t qdb)
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
  add_xml_resblocks(myroot,&(qdb->rb));
  for(i=0; (i<qdb->ngqh); i++)
    add_xml_qhop(myroot,qdb->gqh[i]);

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving qhop xml file",fn);
  xmlFreeDoc(doc);
}



#else

qhop_t qhops_read(char *fn,int *nqhop)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function qhops_read to work");
  return NULL;
}

void qhops_write(char *fn,int nqhop,qhop_t qht)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function qhops_write to work");
}

#endif
