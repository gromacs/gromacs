/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_wnblist_c = "$Id$";

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
  int i,j,j0,k,i_atom,jid;
  
  fprintf(out,"nri=%d  nrj=%d\n",nblist->nri,nblist->nrj);
  for(i=0; i<nblist->nri; i++) {
    fprintf(out,"i: %d, shift: %d, gid: %d\n",
	    nblist->iinr[i],nblist->shift[i],nblist->gid[i]);
    for(j=nblist->jindex[i]; (j<nblist->jindex[i+1]); j++)
      fprintf(out,"  j: %d\n",nblist->jjnr[j]);
  }
  fflush(out);
}

void read_nblist(FILE *in,FILE *log,int **mat,int natoms)
{
  char buf[256];
  int  i,ii,j,nnbl,full,icmp,nri;
  int  iatom,nrj,nj,shift;
  
  do {
    if (fgets2(buf,255,in) == NULL)
      fatal_error(0,"EOF when looking for '%s' in logfile",header);
  } while (strstr(buf,header) == NULL);
  if (fscanf(in,"%d",&nnbl) != 1)
    fatal_error(0,"Not enough arguments read line %d",__LINE__);
  
  for(i=0; (i<nnbl); i++) {
    if (fscanf(in,"%*s%d",&nri) != 1)
      fatal_error(0,"Not enough arguments read line %d",__LINE__);
    nrj = 0;
    for(ii=0; (ii<nri); ii++) {
      if (fscanf(in,"%*s%d%*s%d%*s%d",&iatom,&nj,&shift) != 3)
	fatal_error(0,"Not enough arguments read line %d",__LINE__);
      if ((iatom < 0) || (iatom >= natoms))
	fatal_error(0,"iatom = %d (max %d)\n",iatom,natoms);
      nrj+=nj;
      fscanf(in,"%*s%*s");
      for(i=0; (i<nj); i++) {
	if (fscanf(in,"%d",&j) != 1)
	  fatal_error(0,"Not enough arguments read line %d",__LINE__);
	if ((j < 0) || (j >= natoms))
	  fatal_error(0,"iatom = %d (max %d)\n",iatom,natoms);
	if (mat[iatom][j] != 0)
	  fprintf(log,"mat[%d][%d] changing from %d to %d\n",
		  i,j,mat[iatom][j],shift+1);
	mat[iatom][j] = shift+1;
      }
    }
    fprintf(log,"nri = %d  nrj = %d\n",nri,nrj);
  }
}

void dump_nblist(FILE *out,t_forcerec *fr,int nDNL)
{
  int  i;
  
  fprintf(out,"%s\n",header);

  for(i=0; (i<eNL_NR); i++) 
    write_nblist(out,&fr->nlist_sr[i]);
}

