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
#include "wnblist.h"
#include "fatal.h"
#include "macros.h"
#include "futil.h"

#define header "Neighborlist:"

static void write_nblist(FILE *out,t_nblist *nblist)
{
  int i,j,j0,k,i_atom,jid,nj;
  fprintf(out,"il_code: %d solvent: %d\n",nblist->il_code,nblist->solvent);
  
  fprintf(out,"nri: %d  nrj: %d\n",nblist->nri,nblist->nrj);
  for(i=0; i<nblist->nri; i++) {
    nj = nblist->jindex[i+1] - nblist->jindex[i];
    fprintf(out,"i: %d shift: %d gid: %d nj: %d\n",
	    nblist->iinr[i],nblist->shift[i],nblist->gid[i],nj);
    for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++)
      fprintf(out,"  j: %d\n",nblist->jjnr[j]);
  }
  fflush(out);
}

int read_nblist(FILE *in,FILE *log,int **mat,int natoms)
{
  bool bNL;
  char buf[256],b1[32],b2[32];
  int  i,ii,j,nnbl,full,icmp,nri,il_code,solv;
  int  iatom,nrj,nj,shift,gid,nargs,njtot=0;
  
  do {
    if (fgets2(buf,255,in) == NULL)
      gmx_fatal(FARGS,"EOF when looking for '%s' in logfile",header);
  } while (strstr(buf,header) == NULL);
  
  do {
    if ((nargs = fscanf(in,"%*s%d%*s%d",&il_code,&solv)) != 2)
      return njtot;
    /* gmx_fatal(FARGS,"Can not read il_code or sol (nargs=%d)",nargs);*/
    if ((nargs = fscanf(in,"%*s%d%*s%d",&nri,&nrj)) != 2)
      gmx_fatal(FARGS,"Can not read nri or nrj (nargs=%d)",nargs);
    for(ii=0; (ii<nri); ii++) {
      if ((nargs = fscanf(in,"%*s%d%*s%d%*s%d%*s%d",
			  &iatom,&shift,&gid,&nj)) != 4)
	gmx_fatal(FARGS,"Can not read iatom, shift gid or nj (nargs=%d)",nargs);
      /* Number shifts from 1 to 27 iso 0 to 26, to distinguish uninitialized 
       * matrix elements.
       */
      shift+=1; 
      if ((iatom < 0) || (iatom >= natoms))
	gmx_fatal(FARGS,"iatom = %d (max %d)\n",iatom,natoms);
      nrj+=nj;
      for(i=0; (i<nj); i++) {
	if ((nargs = fscanf(in,"%*s%d",&j)) != 1)
	  gmx_fatal(FARGS,"Can not read j");
	range_check(j,0,natoms);
	if (mat[iatom][j] != 0)
	  fprintf(log,"mat[%d][%d] changing from %d to %d\n",
		  i,j,mat[iatom][j],shift);
	mat[iatom][j] = shift;
	njtot++;
      }
    }
    fprintf(log,"nri = %d  nrj = %d\n",nri,nrj);
  } while (TRUE);
  return -1;
}

void dump_nblist(FILE *out,t_forcerec *fr,int nDNL)
{
  int  i;
  
  fprintf(out,"%s\n",header);

  for(i=0; (i<eNL_NR); i++) 
    write_nblist(out,&fr->nlist_sr[i]);
}

