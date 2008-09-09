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
#include "molprop.h"
#include "molprop_xml.h"

static int comp_mp(const void *a,const void *b)
{
  gmx_molprop_t ma = (gmx_molprop_t)a;
  gmx_molprop_t mb = (gmx_molprop_t)b;
  
  return strcasecmp(gmx_molprop_get_molname(ma),
		    gmx_molprop_get_molname(mb));
}

static int comp_mp_formula(const void *a,const void *b)
{
  int r;
  gmx_molprop_t ma = (gmx_molprop_t)a;
  gmx_molprop_t mb = (gmx_molprop_t)b;
  
  r = strcasecmp(gmx_molprop_get_formula(ma),
		 gmx_molprop_get_formula(mb));
  
  if (r == 0) 
    return comp_mp(a,b);
  else 
    return r;
}

static void merge_doubles(int *np,gmx_molprop_t mp[],char *doubles)
{
  int i,j,ndouble=0;
  FILE *fp;
  
  fp = fopen(doubles,"w");
  for(i=1; (i<*np); i++) {
    if (strcasecmp(gmx_molprop_get_molname(mp[i]),
		   gmx_molprop_get_molname(mp[i-1])) == 0) {
      if (strcasecmp(gmx_molprop_get_formula(mp[i]),
		     gmx_molprop_get_formula(mp[i-1])) == 0) {
	fprintf(fp,"%5d  %s\n",ndouble+1,
		gmx_molprop_get_molname(mp[i-1]));
	gmx_molprop_merge(mp[i-1],mp[i]);
	for(j=i+1; (j<*np); j++) {
	  gmx_molprop_delete(mp[j-1]);
	  mp[j-1] = gmx_molprop_copy(mp[j]);
	}
	ndouble++;
	(*np)--;
      }
      else {
	printf("Molecules %s, %s have formulae %s resp. %s\n",
	       gmx_molprop_get_molname(mp[i]),
	       gmx_molprop_get_formula(mp[i]));
      }
      
    }
  }
  fclose(fp);
  printf("There were %d double entries\n",ndouble);
}

static void dump_mp(int np,gmx_molprop_t mp[])
{
  FILE *fp;
  int  i,j,k;
  
  fp = fopen("dump_mp.dat","w");
  
  for(i=0; (i<np); ) {
    for(j=i; (j<np-1) && (strcasecmp(gmx_molprop_get_formula(mp[i]),
				     gmx_molprop_get_formula(mp[j+1])) == 0); j++)
      ;
    if (j > i) {
      for(k=i; (k<=j); k++)
	fprintf(fp,"%-20s  %s\n",
		gmx_molprop_get_formula(mp[k]),
		gmx_molprop_get_molname(mp[k]));
      fprintf(fp,"\n");
    }
    i=j+1;
  }
  
  fclose(fp);
}

gmx_molprop_t *merge_xml(int argc,char *argv[],char *outf,
			 char *sorted,char *doubles,int *nmolprop)
{
  gmx_molprop_t *mp=NULL,*mpout=NULL;
  int       i,j,np,npout=0;
  char      buf[100];
  
  for(i=1; (i<argc); i++) {
    mp = read_molprops(argv[i],&np);
    srenew(mpout,npout+np);
    for(j=0; (j<np); j++)
      mpout[npout+j] = gmx_molprop_copy(mp[j]);
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
