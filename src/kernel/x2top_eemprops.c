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
#include "x2top_qgen.h"
#include "x2top_eemprops.h"

typedef struct {
  char *name,*opts;
  int  eemtype,elem,row;
  real J0,w,chi0; 
  /* J0 in Yang & Sharp corresponds to n (eta) in Bultinck */
} t_eemprops;

typedef struct {
  int        nep;
  t_eemprops *eep;
} t_eemrecord;

typedef struct {
  char *name,*ref;
  bool bWeight;
} t_eemtype_props;

static t_eemtype_props eemtype_props[eqgNR] = { 
  { "None",     "None",	         FALSE },   
  { "Yang", 	"Yang2006b",     TRUE },    
  { "Bultinck", "Bultinck2002a", FALSE },
  { "Rappe",    "Rappe1991a",	 TRUE },  
  { "SMp", 	"Spoel2008b",	 FALSE },    
  { "SMpp", 	"Spoel2008b",    FALSE },
  { "SMs", 	"Spoel2008b",    TRUE },
  { "SMps", 	"Spoel2008b",    TRUE },
  { "SMg", 	"Spoel2008b",    TRUE },
  { "SMpg", 	"Spoel2008b",    TRUE }
};

int name2eemtype(char *name)
{
  int i;
  
  for(i=0; (i<eqgNR); i++) {
    if (strcasecmp(name,eemtype_props[i].name) == 0)
      return i;
  }
  return -1;
}

char *get_eemtype_name(int eem)
{
  int i;
  
  if ((eem >= 0) && (eem < eqgNR))
    return eemtype_props[eem].name;
    
  return NULL;
}

char *get_eemtype_reference(int eem)
{
  int i;
  
  if ((eem >= 0) && (eem < eqgNR))
    return eemtype_props[eem].ref;
    
  return NULL;
}

void *read_eemprops(char *fn,int eemtype,void *atomprop)
{
  t_eemrecord *eem=NULL;
  char   buf[STRLEN],**strings,*ptr;
  int    i,n,narg,nn=0;
  char   nmbuf[32],algbuf[32],optbuf[32];
  int    elem,row;
  real   value;
  double J0,w,chi0;
  
  if (fn == NULL) 
    sprintf(buf,"eemprops.dat");
  else
    strcpy(buf,fn);
  n  = get_file(buf,&strings);
  if (n > 0) {
    snew(eem,1);
    snew(eem->eep,n);
    for(i=0; (i<n); i++) {
      ptr = strings[i];
      while (*ptr && isspace(*ptr))
	ptr++;
      optbuf[0] = '\0';
      if (((ptr) && (*ptr != ';')) &&
	  ((narg = sscanf(strings[i],"%s%s%lf%lf%lf%s",nmbuf,algbuf,
			  &J0,&chi0,&w,optbuf)) >= 4))  {
	if ((eem->eep[nn].eemtype = name2eemtype(algbuf)) == -1)
	  fprintf(stderr,"Warning in %s on line %d, unknown algorithm '%s'\n",
		  buf,i+1,algbuf);
	else if ((eemtype == -1) || (eem->eep[nn].eemtype == eemtype)) {
	  eem->eep[nn].name    = strdup(nmbuf);
	  if (!query_atomprop(atomprop,epropElement,"???",nmbuf,&value))
	    gmx_fatal(FARGS,"Can not find element type for atom %s",nmbuf);
	  elem = gmx_nint(value);
	  eem->eep[nn].elem  = elem;
	  /* Compute which row in the periodic table is this element */
	  if (elem <= 2)
	    row = 1;
	  else if (elem <= 10)
	    row = 2;
	  else if (elem <= 18)
	    row = 3;
	  else if (elem <= 36)
	    row = 4;
	  else if (elem <= 54)
	    row = 5;
	  else if (elem <= 86)
	    row = 6;
	  else
	    row = 7;

	  eem->eep[nn].row   = row;
	  eem->eep[nn].J0    = J0;
	  eem->eep[nn].chi0  = chi0;
	  
	  if (narg > 4)
	    eem->eep[nn].w     = w;
	  else
	    eem->eep[nn].w     = 0;
	  
	  if (strlen(optbuf) > 0)
	    eem->eep[nn].opts = strdup(optbuf);
	  else
	    eem->eep[nn].opts = NULL;
	  nn++;
	}
      }
    }
    eem->nep = nn;
  }
    
  return eem;
}

static void write_eemprops_header(FILE *fp,int eemtype)
{
  switch (eemtype) {
  case eqgYang:
    fprintf(fp,";\n");
    fprintf(fp,"; Parameters for Yang & Sharp, JCTC 2 (2006) 1152\n");
    fprintf(fp,"; Atom      Model    J_aa (eV)  chi_a (eV)    R_a (A)\n");
    break;
  case eqgRappe:
    fprintf(fp,";\n");
    fprintf(fp,"; Parameters for Rappe & Goddard, JPC 95 (1991) 3358\n");
    fprintf(fp,"; Atom      Model    J_aa (eV)  chi_a (eV) z_a (1/au)\n");
    break;
  case eqgBultinck:
    fprintf(fp,";\n");
    fprintf(fp,"; Parameters for Bultinck et al., JPCA 106 (2002) 7887\n");
    fprintf(fp,"; J_aa equals Hardness, this is multiplied by 2 in x2top.\n");
    fprintf(fp,"; Atom      Model    J_aa (eV) chi_a (eV)\n");
    break;
  
  case eqgSMp:
  case eqgSMpp:
    fprintf(fp,";\n");
    fprintf(fp,"; Parameters for van Maaren & van der Spoel\n");
    fprintf(fp,"; Atom      Model   J_aa (eV)  chi_a (eV)\n");
    break;
  case eqgSMg:
  case eqgSMpg:
  case eqgSMs:
  case eqgSMps:
    fprintf(fp,";\n");
    fprintf(fp,"; Parameters for van Maaren & van der Spoel\n");
    fprintf(fp,"; Atom      Model   J_aa (eV)  chi_a (eV)  z_a (1/nm)\n");
    break;
  default:
    break;
  }
}

void write_eemprops(FILE *fp,void *eem)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  double w;
  int i;

  fprintf(fp,"; Parameters for electronegativity algorithms. This file is used by x2top.\n");
  fprintf(fp,"; Note that parameters may have different meaning and different units.\n");
  for(i=0; (i<er->nep); i++) {
    if ((i == 0) || 
	((i > 0) && (er->eep[i].eemtype != er->eep[i-1].eemtype)))
      write_eemprops_header(fp,er->eep[i].eemtype);
    w = er->eep[i].w;
    fprintf(fp,"%-5s  %10s  %10.4f  %10.5f",
	    er->eep[i].name,eemtype_props[er->eep[i].eemtype].name,
	    er->eep[i].J0,er->eep[i].chi0);
    if (eemtype_props[er->eep[i].eemtype].bWeight)
      fprintf(fp,"  %10.5f\n",w);
    else
      fprintf(fp,"\n");
  }
}

int eem_get_numprops(void *eem,int eemtype)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i,n=0;
  
  for(i=0; (i<er->nep); i++) {
    if (er->eep[i].eemtype == eemtype)
      n++;
  }
  return n;
}

int eem_get_index(void *eem,int atomicnumber,int eemtype)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i;
  
  for(i=0; (i<er->nep); i++) 
    if ((er->eep[i].elem == atomicnumber) && 
	(er->eep[i].eemtype == eemtype))
      return i;
  return -1;
}

int eem_get_elem_index(void *eem,int atomicnumber,int eemtype)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i;
  
  for(i=0; (i<er->nep); i++) 
    if ((er->eep[i].elem == atomicnumber) && 
	(er->eep[i].eemtype == eemtype))
      return i;
  return -1;
}

int eem_get_row(void *eem,int index)
{
  t_eemrecord *er = (t_eemrecord *) eem;

  range_check(index,0,er->nep);
  
  return er->eep[index].row;
}

real eem_get_j00(void *eem,int index,real *wj,real q)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  double Ra,zeta,zeta0,j0;

  /* Bohr is in nm */
#define BOHR  (0.052917)
  /* Weird fudge factor lamba, see Rappe & Goddard, 
     JPC 95 (1991) p.3358 */
#define lambda 0.4913

  range_check(index,0,er->nep);
  j0 = er->eep[index].J0;
  if (er->eep[index].eemtype == eqgYang) {
    /* Convert Ra from Angstrom to au */
    Ra = er->eep[index].w/(10*BOHR);
    if (er->eep[index].elem == 1) {
      zeta0 = (3/(4*Ra));
      j0    = (1+q/zeta0)*j0;
      zeta  = zeta0 /*+ q*/;
    }
    else {
      zeta = lambda*(2*er->eep[index].row+1)/(2*Ra);
    }
    *wj = zeta/BOHR;
  }
  else if (er->eep[index].eemtype == eqgRappe) {
    zeta = er->eep[index].w;
    if (er->eep[index].elem == 1) {
      j0 = (1+q/zeta)*j0;
    }
    *wj = zeta/BOHR;
  }
  else if ((er->eep[index].eemtype == eqgSMs) || 
	   (er->eep[index].eemtype == eqgSMps) ||
	   (er->eep[index].eemtype == eqgSMg) ||
	   (er->eep[index].eemtype == eqgSMpg)) {
    zeta = er->eep[index].w*BOHR;
    if (er->eep[index].elem == 1) {
      j0 = (1+q/zeta)*j0;
    }
    *wj = zeta/BOHR;
  }
  else
    *wj = 0;
    
  return j0;
}

char *eem_get_opts(void *eem,int index)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);
  
  return er->eep[index].opts;
}

real eem_get_chi0(void *eem,int index)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);
  
  return er->eep[index].chi0;
}

real eem_get_w(void *eem,int index)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);
  
  return er->eep[index].w;
}

int eem_get_elem(void *eem,int index)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);
  
  return er->eep[index].elem;
}

void eem_set_props(void *eem,int index,real J0,real w,real chi0)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);
  
  er->eep[index].J0   = J0;
  er->eep[index].w    = w;
  er->eep[index].chi0 = chi0;
}

void *copy_eem(void *eem_dst,void *eem_src)
{
  t_eemrecord *dst = (t_eemrecord *) eem_dst;
  t_eemrecord *src = (t_eemrecord *) eem_src;
  int i;
  
  if (dst == NULL) {
    snew(dst,1);
    dst->nep = src->nep;
    snew(dst->eep,dst->nep);
    for(i=0; (i<dst->nep); i++) {
      dst->eep[i].name = strdup(src->eep[i].name);
      if (src->eep[i].opts)
	dst->eep[i].opts = strdup(src->eep[i].opts);
      else
	dst->eep[i].opts = NULL;
    }
  }
  else if (dst->nep != src->nep) 
    gmx_fatal(FARGS,"src->nep = %d is different from dst->nep = %d",
	      src->nep,dst->nep);
  for(i=0; (i<dst->nep); i++) {
    dst->eep[i].eemtype = src->eep[i].eemtype;
    dst->eep[i].elem    = src->eep[i].elem;
    dst->eep[i].J0      = src->eep[i].J0;
    dst->eep[i].w       = src->eep[i].w;
    dst->eep[i].chi0    = src->eep[i].chi0;
  }
  return dst;
}
