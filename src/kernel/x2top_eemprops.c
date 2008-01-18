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
#include "vec.h"
#include "atomprop.h"
#include "grompp.h"
#include "x2top_qgen.h"
#include "x2top_eemprops.h"

typedef struct {
  char *name;
  int  eemtype,elem,row;
  real J0,w,chi0; 
  /* J0 in Yang & Sharp corresponds to n (eta) in Bultinck */
} t_eemprops;

typedef struct {
  int        nep;
  t_eemprops *eep;
} t_eemrecord;

static char *eemtype_name[eqgNR] = { 
  "None", "Linear", "Yang", "Bultinck", "SM1", "SM2", "SM3", "SM4" 
};

int name2eemtype(char *name)
{
  int i;
  
  for(i=0; (i<eqgNR); i++) {
    if (strcasecmp(name,eemtype_name[i]) == 0)
      return i;
  }
  return -1;
}

void *read_eemprops(char *fn,int eemtype)
{
  t_eemrecord *eem=NULL;
  char   buf[STRLEN],**strings,*ptr;
  int    i,n,nn=0;
  char   nmbuf[32],algbuf[32];
  int    elem,row;
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
      if (((ptr) && (*ptr != ';')) &&
	  (sscanf(strings[i],"%s%s%d%d%lf%lf%lf",nmbuf,algbuf,&elem,&row,
		  &J0,&w,&chi0) == 7))  {
	if ((eem->eep[nn].eemtype = name2eemtype(algbuf)) == -1)
	  fprintf(stderr,"Warning in %s on line %d, unknown algorithm '%s'\n",
		  buf,i+1,algbuf);
	else if ((eemtype == -1) || (eem->eep[nn].eemtype == eemtype)) {
	  eem->eep[nn].name    = strdup(nmbuf);
	  eem->eep[nn].elem    = elem;
	  eem->eep[nn].row     = row;
	  eem->eep[nn].J0      = J0;
	  eem->eep[nn].w  = w;
	  eem->eep[nn].chi0    = chi0;
	  nn++;
	}
      }
    }
    eem->nep = nn;
  }
    
  return eem;
}

void write_eemprops(FILE *fp,void *eem)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i;

  fprintf(fp,"; Electronegativity parameters. J_aa and Chi_a are in eV, w_a in nm\n");
  fprintf(fp,"; Atom      Model   Nr  Row        J_aa         w_a       Chi_a\n");
  for(i=0; (i<er->nep); i++)
    fprintf(fp,"%-5s  %10s  %3d  %3d  %10.4f  %10.4f  %10.4f\n",
	    er->eep[i].name,eemtype_name[er->eep[i].eemtype],
	    er->eep[i].elem,er->eep[i].row,er->eep[i].J0,
	    er->eep[i].w,er->eep[i].chi0);
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

int eem_get_index(void *eem,char *resname,char *aname,int eemtype)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i;
  
  for(i=0; (i<er->nep); i++) 
    if ((strstr(aname,er->eep[i].name) == aname) && 
	(er->eep[i].eemtype == eemtype))
      return i;
  return -1;
}

int eem_get_elem_index(void *eem,int elem,int eemtype)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  int i;
  
  for(i=0; (i<er->nep); i++) 
    if ((er->eep[i].elem == elem) && (er->eep[i].eemtype == eemtype))
      return i;
  return -1;
}

real lo_get_j00(void *eem,int index,real *wj,real qH)
{
  t_eemrecord *er = (t_eemrecord *) eem;
  
  range_check(index,0,er->nep);

  if (er->eep[index].eemtype == eqgYang) {
    if (er->eep[index].elem == 1) 
      *wj = (3/(4*er->eep[index].w)+10*qH);
    else 
      *wj = (3/(4*er->eep[index].w));
  }
  else if ((er->eep[index].eemtype == eqgSM2) || 
	   (er->eep[index].eemtype == eqgSM3) ||
	   (er->eep[index].eemtype == eqgSM4))
    *wj = 1.0/er->eep[index].w;
  else
    *wj = 0;
    
  return er->eep[index].J0;
}

real eem_get_j00(void *eem,char *resname,char *aname,real *wj,real qH,int eemtype)
{
  int k = eem_get_index(eem,resname,aname,eemtype);

  return lo_get_j00(eem,k,wj,qH);
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
