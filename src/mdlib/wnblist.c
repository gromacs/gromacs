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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include "string2.h"
#include "force.h"
#include "smalloc.h"
#include "ns.h"
#include "nrnb.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "futil.h"
#include "names.h"
#include "domdec.h"

#define header "Neighborlist:"

static void write_nblist(FILE *out,gmx_domdec_t *dd,t_nblist *nblist,int nDNL)
{
  int i,j,ci,cj0,cj1,aj,cj,nj;
  int ca1[DD_MAXCELL],np[DD_MAXCELL];

  if (nblist->nri > 0) {  
    fprintf(out,"il_name: %s  Solvent opt: %s\n",
            nrnb_str(nblist->il_code),
            enlist_names[nblist->solvent_opt]);
    fprintf(out,"nri: %d  npair: %d\n",nblist->nri,nblist->nrj);
    if (dd) {
      for(ci=0; ci<dd->ncell; ci++)
	ca1[ci] = dd->cgindex[dd->ncg_cell[ci+1]];
      i = 0;
      for(ci=0; ci<dd->nicell; ci++) {
	cj0 = dd->icell[ci].j0;
	cj1 = dd->icell[ci].j1;
	for(cj=cj0; cj<cj1; cj++)
	  np[cj] = 0;
	while(i < nblist->nri && nblist->iinr[i] < ca1[ci]) {
	  for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++) {
	    aj = nblist->jjnr[j];
	    cj = cj0;
	    while (aj >= ca1[cj])
	      cj++;
	    np[cj]++;
	  }
	  i++;
	}
	fprintf(out,"DD cell %d:",ci);
	for(cj=cj0; cj<cj1; cj++)
	  fprintf(out," %d %d",cj,np[cj]);
	fprintf(out,"\n");
      }
    }
    if (nDNL == 2) {
      for(i=0; i<nblist->nri; i++) {
	nj = nblist->jindex[i+1] - nblist->jindex[i];
	fprintf(out,"i: %d shift: %d gid: %d nj: %d\n",
		glatnr(dd,nblist->iinr[i]),nblist->shift[i],nblist->gid[i],nj);
	for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++)
	  fprintf(out,"  j: %d\n",glatnr(dd,nblist->jjnr[j]));
      }
    }
    fflush(out);
  }
}

static void set_mat(FILE *fp,int **mat,int i0,int ni,int j0,int nj,
		    bool bSymm,int shift)
{
  int i,j;
  
  for(i=i0; (i<i0+ni); i++) {
    for(j=j0; (j<j0+nj); j++) {
      if (mat[i][j] != 0)
	fprintf(fp,"mat[%d][%d] changing from %d to %d\n",
		i,j,mat[i][j],shift+1);
      mat[i][j] = shift+1;
      if (bSymm)
	mat[j][i] = 27-shift;
    }
  }
}

int read_nblist(FILE *in,FILE *fp,int **mat,int natoms,bool bSymm)
{
    bool bNL;
    char buf[256],b1[32],b2[32],solv[256],il_code[256];
    int  i,ii,j,nnbl,full,icmp,nri,isolv;
    int  iatom,nrj,nj,shift,gid,nargs,njtot=0;
    
    do {
        if (fgets2(buf,255,in) == NULL)
            gmx_fatal(FARGS,"EOF when looking for '%s' in logfile",header);
    } while (strstr(buf,header) == NULL);
    
    do {
      do {
        if (fgets2(buf,255,in) == NULL)
	  return njtot;
      } while (strstr(buf,"nri:") == NULL);
      
      if (0) {
	if ((nargs = sscanf(buf,"%*s%s%*s%s",il_code,solv)) != 2) {
	  fprintf(stderr,"Can not find the right il_code\n");
	  return njtot;
	}
        for(isolv=0; (isolv<esolNR); isolv++)
	  if (strstr(esol_names[isolv],solv) != NULL)
	    break;
	
        if (isolv == esolNR) {
	  fprintf(stderr,"Can not read il_code or solv (nargs=%d)\n",nargs);
	  return njtot;
        }
      }
      else
	isolv = enlistATOM;
      
      /* gmx_fatal(FARGS,"Can not read il_code or solv (nargs=%d)",nargs);*/
      if ((nargs = sscanf(buf,"%*s%d%*s%d",&nri,&nrj)) != 2)
	gmx_fatal(FARGS,"Can not read nri or nrj (nargs=%d)",nargs);
      for(ii=0; (ii<nri); ii++) {
	if ((nargs = fscanf(in,"%*s%d%*s%d%*s%d%*s%d",
			    &iatom,&shift,&gid,&nj)) != 4)
	  gmx_fatal(FARGS,"Can not read iatom, shift gid or nj (nargs=%d)",nargs);
	/* Number shifts from 1 to 27 iso 0 to 26 to distinguish uninitialized 
	 * matrix elements.
	 */
	range_check(iatom,0,natoms);
	for(i=0; (i<nj); i++) {
	  if ((nargs = fscanf(in,"%*s%d",&j)) != 1)
	    gmx_fatal(FARGS,"Can not read j");
	  range_check(j,0,natoms);
	  switch (isolv) {
	  case enlistATOM:
	    set_mat(fp,mat,iatom,1,j,1,bSymm,shift);
	    njtot++;
	    break;
	  case enlistWATER:
	    set_mat(fp,mat,iatom,3,j,1,bSymm,shift);
	    njtot+=3;
	    break;
	  case enlistWATERWATER:
	    set_mat(fp,mat,iatom,3,j,3,bSymm,shift);
	    njtot+=9;
	    break;
	  default:
	    gmx_incons("non-existing solvent type");
	  }
	}
      }
      fprintf(fp,"nri = %d  nrj = %d\n",nri,nrj);
    } while (TRUE);
    return -1;
}

void dump_nblist(FILE *out,t_commrec *cr,t_forcerec *fr,int nDNL)
{
  int  n,i;
  
  fprintf(out,"%s\n",header);

  for(n=0; (n<fr->nnblists); n++)
    for(i=0; (i<eNL_NR); i++) 
      write_nblist(out,cr->dd,&fr->nblists[n].nlist_sr[i],nDNL);
}

