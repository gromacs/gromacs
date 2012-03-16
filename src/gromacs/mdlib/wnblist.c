/*
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
#include "gmxfio.h"

#define header "Neighborlist:"

static void write_nblist(FILE *out,gmx_domdec_t *dd,t_nblist *nblist,int nDNL)
{
  int i,nii,ii,j,zi,zj0,zj1,aj,zj,nj;
  int ca1[DD_MAXZONE],np[DD_MAXZONE];
  gmx_domdec_zones_t *dd_zones;

  if (nblist->nri > 0) {  
    fprintf(out,"il_name: %s  Solvent opt: %s\n",
            nrnb_str(nblist->il_code),
            enlist_names[nblist->enlist]);
    fprintf(out,"nri: %d  npair: %d\n",nblist->nri,nblist->nrj);
    if (dd) {
      dd_zones = domdec_zones(dd);

      for(zi=0; zi<dd_zones->n; zi++)
	ca1[zi] = dd->cgindex[dd_zones->cg_range[zi+1]];
      i = 0;
      for(zi=0; zi<dd_zones->nizone; zi++) {
	zj0 = dd_zones->izone[zi].j0;
	zj1 = dd_zones->izone[zi].j1;
	for(zj=zj0; zj<zj1; zj++)
	  np[zj] = 0;
	while(i < nblist->nri && nblist->iinr[i] < ca1[zi]) {
	  for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++) {
	    aj = nblist->jjnr[j];
	    zj = zj0;
	    while (aj >= ca1[zj])
	      zj++;
	    np[zj]++;
	  }
	  i++;
	}
	fprintf(out,"DD zone %d:",zi);
	for(zj=zj0; zj<zj1; zj++)
	  fprintf(out," %d %d",zj,np[zj]);
	fprintf(out,"\n");
      }
    }
    if (nDNL >= 2) {
      for(i=0; i<nblist->nri; i++) {
	nii = 1;
	if (nDNL >= 3 && nblist->enlist != enlistATOM_ATOM)
	  nii = 3;
	nj = nblist->jindex[i+1] - nblist->jindex[i];
	fprintf(out,"i: %d shift: %d gid: %d nj: %d\n",
		ddglatnr(dd,nblist->iinr[i]),
		nblist->shift[i],nblist->gid[i],nj);
	for(ii=0; ii<nii; ii++) {
	  for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++) {
	    fprintf(out,"  i: %5d  j: %5d\n",
		    ddglatnr(dd,nblist->iinr[i]+ii),
		    ddglatnr(dd,nblist->jjnr[j]));
	  }
	}
      }
    }
    fflush(out);
  }
}

static void set_mat(FILE *fp,int **mat,int i0,int ni,int j0,int nj,
		    gmx_bool bSymm,int shift)
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

int read_nblist(FILE *in,FILE *fp,int **mat,int natoms,gmx_bool bSymm)
{
    gmx_bool bNL;
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
	isolv = enlistATOM_ATOM;
      
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
	  case enlistATOM_ATOM:
	    set_mat(fp,mat,iatom,1,j,1,bSymm,shift);
	    njtot++;
	    break;
	  case enlistSPC_ATOM:
	    set_mat(fp,mat,iatom,3,j,1,bSymm,shift);
	    njtot+=3;
	    break;
	  case enlistSPC_SPC:
	    set_mat(fp,mat,iatom,3,j,3,bSymm,shift);
	    njtot+=9;
	    break;
	  case enlistTIP4P_ATOM:
	    set_mat(fp,mat,iatom,4,j,1,bSymm,shift);
	    njtot+=4;
	    break;
	  case enlistTIP4P_TIP4P:
	    set_mat(fp,mat,iatom,4,j,4,bSymm,shift);
	    njtot+=16;
	    break;
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
#if 0
  static FILE *fp=NULL;
  char buf[STRLEN];
  int  n,i;

  if (fp == NULL) {
    if (PAR(cr)) {
      sprintf(buf,"nlist_n%d.txt",cr->nodeid);
    } else {
      sprintf(buf,"nlist.txt");
    }
    fp = gmx_fio_fopen(buf,"w");
  }
  fprintf(fp,"%s\n",header);

  for(n=0; (n<fr->nnblists); n++)
    for(i=0; (i<eNL_NR); i++) 
      write_nblist(fp,cr->dd,&fr->nblists[n].nlist_sr[i],nDNL);
#endif
  char buf[STRLEN];
  int  n,i;

  fprintf(out,"%s\n",header);

  for(n=0; (n<fr->nnblists); n++)
    for(i=0; (i<eNL_NR); i++) 
      write_nblist(out,cr->dd,&fr->nblists[n].nlist_sr[i],nDNL);

}

