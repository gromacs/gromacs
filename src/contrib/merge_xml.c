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
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "smalloc.h"
#include "tune_pol.h"

static int comp_mp(const void *a,const void *b)
{
  t_molprop *ma = (t_molprop *)a;
  t_molprop *mb = (t_molprop *)b;
  
  return strcasecmp(ma->molname,mb->molname);
}

static int comp_mp_formula(const void *a,const void *b)
{
  int r;
  
  t_molprop *ma = (t_molprop *)a;
  t_molprop *mb = (t_molprop *)b;
  
  r = strcasecmp(ma->formula,mb->formula);
  if (r == 0) 
    return strcasecmp(ma->molname,mb->molname);
  else 
    return r;
}

static void copy_molprop(t_molprop *dst,t_molprop *src)
{
  int  i;
  char *ptr;
  
  memcpy(dst,src,sizeof(*src));
  dst->molname = strdup(src->molname);
  while ((ptr = strchr(dst->molname,' ')) != NULL)
    *ptr = '-';
  dst->formula = strdup(src->formula);
  snew(dst->experiment,dst->nexperiment);
  snew(dst->reference,dst->nexperiment);
  snew(dst->pname,dst->nexperiment);
  for(i=0; (i<dst->nexperiment); i++) {
    dst->experiment[i] = src->experiment[i];
    dst->reference[i] = strdup(src->reference[i]);
    dst->pname[i] = strdup(src->pname[i]);
  }
}

static void merge_molprop(t_molprop *dst,t_molprop *src)
{
  int i,nd,ns;
  
  srenew(dst->experiment,dst->nexperiment+src->nexperiment);
  srenew(dst->reference,dst->nexperiment+src->nexperiment);
  srenew(dst->pname,dst->nexperiment+src->nexperiment);
  for(i=dst->nexperiment; (i<dst->nexperiment+src->nexperiment); i++) {
    dst->experiment[i] = src->experiment[i-dst->nexperiment];
    dst->reference[i] = strdup(src->reference[i-dst->nexperiment]);
    dst->pname[i] = strdup(src->pname[i-dst->nexperiment]);
  }
  dst->nexperiment+=src->nexperiment;
  /* Compare src and dst for composition etc. */
  nd = ns = 0;
  for(i=0; (i<eatNR+eatExtra); i++) {
    if (src->frag_comp[i] > 0)
      ns++;
    if (dst->frag_comp[i] > 0)
      nd++;
  }
  for(i=0; (i<eelemNR); i++) {
    if (src->elem_comp[i] > 0)
      ns++;
    if (dst->elem_comp[i] > 0)
      nd++;
  }
  for(i=0; (i<emlNR); i++) {
    if (src->emil_comp[i] > 0)
      ns++;
    if (dst->emil_comp[i] > 0)
      nd++;
  }
  if (ns > 0) {
    if (nd == 0) {
      for(i=0; (i<eatNR+eatExtra); i++) 
	dst->frag_comp[i] = src->frag_comp[i];
      for(i=0; (i<eelemNR); i++) 
	dst->elem_comp[i] = src->elem_comp[i];
      for(i=0; (i<emlNR); i++) 
	dst->emil_comp[i] = src->emil_comp[i];
    }
    else
      printf("Both src and dst for %s contain composition entries. Not changing anything\n",dst->molname);
  }
}

static void clear_molprop(t_molprop *dst)
{
  int i;
  
  for(i=0; (i<dst->nexperiment); i++) {
    sfree(dst->reference[i]);
    sfree(dst->pname[i]);
  }
  if (dst->nexperiment > 0) {
    sfree(dst->experiment);
    sfree(dst->reference);
    sfree(dst->pname);
  }
  dst->nexperiment = 0;
}

static void merge_doubles(int *np,t_molprop mp[],char *doubles)
{
  int i,j,ndouble=0;
  FILE *fp;
  
  fp = fopen(doubles,"w");
  for(i=1; (i<*np); i++) {
    if (strcasecmp(mp[i].molname,mp[i-1].molname) == 0) {
      if (strcasecmp(mp[i].formula,mp[i-1].formula) == 0) {
	fprintf(fp,"%5d  %s\n",ndouble+1,mp[i-1].molname);
	merge_molprop(&(mp[i-1]),&(mp[i]));
	for(j=i+1; (j<*np); j++) {
	  clear_molprop(&(mp[j-1]));
	  copy_molprop(&(mp[j-1]),&(mp[j]));
	}
	ndouble++;
	(*np)--;
      }
      else {
	printf("Molecules %s, %s have formulae %s resp. %s\n",
	       mp[i].molname,mp[i-1].molname,
	       mp[i].formula,mp[i-1].formula);
      }
      
    }
  }
  fclose(fp);
  printf("There were %d double entries\n",ndouble);
}

static void dump_mp(int np,t_molprop mp[])
{
  FILE *fp;
  int  i,j,k;
  
  fp = fopen("dump_mp.dat","w");
  
  for(i=0; (i<np); ) {
    for(j=i; (j<np-1) && (strcasecmp(mp[i].formula,mp[j+1].formula) == 0); j++)
      ;
    if (j > i) {
      for(k=i; (k<=j); k++)
	fprintf(fp,"%-20s  %s\n",mp[k].formula,mp[k].molname);
      fprintf(fp,"\n");
    }
    i=j+1;
  }
  
  fclose(fp);
}

t_molprop *merge_xml(int argc,char *argv[],char *outf,
		     char *sorted,char *doubles,int *nmolprop)
{
  t_molprop *mp=NULL,*mpout=NULL;
  int       i,j,np,npout=0;
  char      buf[100];
  
  for(i=1; (i<argc); i++) {
    np = read_molprops(argv[i],&mp,1);
    srenew(mpout,npout+np);
    for(j=0; (j<np); j++)
      copy_molprop(&(mpout[npout+j]),&(mp[j]));
    npout += np;
  }
  
  qsort(mpout,npout,sizeof(mp[0]),comp_mp);
  merge_doubles(&npout,mpout,doubles);
      
  if (outf) {
    printf("There are %d entries to store in output file %s\n",npout,outf);
    write_molprops(outf,npout,mpout);
  }
  if (sorted) {
    qsort(mpout,npout,sizeof(mp[0]),comp_mp_formula);
    write_molprops(sorted,npout,mpout);
    dump_mp(npout,mpout);
  }
  
  *nmolprop = npout; 
  
  return mpout;
}
