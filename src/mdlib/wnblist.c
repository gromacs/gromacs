/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
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

void read_nblist(FILE *in,FILE *log,int **mat,int natoms)
{
  bool bNL;
  char buf[256],b1[32],b2[32];
  int  i,ii,j,nnbl,full,icmp,nri,il_code,solv;
  int  iatom,nrj,nj,shift,gid;
  
  do {
    if (fgets2(buf,255,in) == NULL)
      fatal_error(0,"EOF when looking for '%s' in logfile",header);
  } while (strstr(buf,header) == NULL);
  
  do {
    if (fscanf(in,"%*s%d%*s%d",&il_code,&solv) != 2)
      break;
    if (fscanf(in,"%*s%d%*s%d",&nri,&nrj) != 2)
      fatal_error(0,"Not enough arguments read line %d",__LINE__);
    for(ii=0; (ii<nri); ii++) {
      if (fscanf(in,"%*s%d%*s%d%*s%d%*s%d",&iatom,&gid,&shift,&nj) != 4)
	fatal_error(0,"Not enough arguments read line %d",__LINE__);
      /* Number shifts from 1 to 27 iso 0 to 26, to distinguish uninitialized 
       * matrix elements.
       */
      shift+=1; 
      if ((iatom < 0) || (iatom >= natoms))
	fatal_error(0,"iatom = %d (max %d)\n",iatom,natoms);
      nrj+=nj;
      for(i=0; (i<nj); i++) {
	if (fscanf(in,"%*s%d",&j) != 1)
	  fatal_error(0,"Not enough arguments read line %d",__LINE__);
	if ((j < 0) || (j >= natoms))
	  fatal_error(0,"iatom = %d (max %d)\n",j,natoms);
	if (mat[iatom][j] != 0)
	  fprintf(log,"mat[%d][%d] changing from %d to %d\n",
		  i,j,mat[iatom][j],shift);
	mat[iatom][j] = shift;
      }
    }
    fprintf(log,"nri = %d  nrj = %d\n",nri,nrj);
  } while (TRUE);
}

void dump_nblist(FILE *out,t_forcerec *fr,int nDNL)
{
  int  i;
  
  fprintf(out,"%s\n",header);

  for(i=0; (i<eNL_NR); i++) 
    write_nblist(out,&fr->nlist_sr[i]);
}

