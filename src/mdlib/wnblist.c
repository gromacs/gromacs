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

static void write_nblist(FILE *out,t_nblist *nblist,rvec sv[SHIFTS],
			 bool bFull)
{
  int j,j0,k,i_atom,jid;
  
  fatal_error(0,"write_nblist temporarily out of order");
  fprintf(out,"nri=%8u\n",nblist->nri);
  for(j=0; j<(int)nblist->nri; j++) {
    /*  i_atom=nblist->nl_i[j].i_atom;
    fprintf(out,"i_atom=%5d  nj=%8u  shift=%4d\n",
	    i_atom,nblist->nl_i[j].nj,nblist->nl_i[j].shift);
    fprintf(out,"%5s  %4s\n","jid","grp");
    j0=nblist->nl_i[j].j_index;
    for (k=0; k<(int)nblist->nl_i[j].nj; k++) {
      jid=nblist->nl_j[j0+k];
      fprintf(out,"%5d\n",jid);
    }
    */
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
  rvec *sv;
  
  fprintf(out,"%s\n",header);
  /*  fprintf(out,"%d\n",fr->nn*2);*/

  sv=fr->shift_vec;
  /*for(i=0; (i<fr->nn); i++) {
    write_nblist(out,&fr->coul[i],sv,(nDNL > 1));
    write_nblist(out,&fr->vdw[i],sv,(nDNL > 1));
    }*/
}

