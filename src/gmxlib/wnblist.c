/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
#include <string.h>
#include <stdio.h>
#include "force.h"
#include "smalloc.h"
#include "wnblist.h"
#include "fatal.h"
#include "macros.h"
#include "futil.h"

static void write_nblist(FILE *out,t_nblist *nblist,rvec sv[SHIFTS],
			 bool bFull)
{
  int i,j,j0,k,i_atom,jid,grp;
  

  fprintf(out,"nri=%8d\n",nblist->nri);
  for(j=0; j<(int)nblist->nri; j++) {
    i_atom=nblist->nl_i[j].i_atom;
    fprintf(out,"i_atom=%5d  nj=%8d  shift=%4d\n",i_atom,
	    nblist->nl_i[j].nj,nblist->nl_i[j].shift);
    fprintf(out,"%5s  %4s\n","jid","grp");
    j0=nblist->nl_i[j].j_index;
    for (k=0; k<(int)nblist->nl_i[j].nj; k++) {
      jid=nblist->nl_j[j0+k];
      fprintf(out,"%5d\n",jid);
    }
  }
  fflush(out);
}

/*void low_readnblist(FILE *in,t_nblist *nbl,int full)
{
  char buf[256];
  int  i,j,ni,ia,nj,t,ja,grp,nrj;
  
  grp=0;
  nrj=0;
  if (fscanf(in,"%*s%d",&ni) != 1)
    fatal_error(0,"Not enough arguments read line %d",__LINE__);
  nbl->nri=ni;
  snew(nbl->nl_i,ni);
  for(i=0; (i<ni); i++) {
    fprintf(stderr,"\rnri: %d",i);
    if (fscanf(in,"%*s%d%*s%d%*s%d",&ia,&nj,&t) != 3) 
      fatal_error(0,"Not enough arguments read line %d",__LINE__);
    nbl->nl_i[i].i_atom=ia;
    nbl->nl_i[i].nj=nj;
    nbl->nl_i[i].shift=t;
    snew(nbl->nl_i[i].nlj,2*nj);
    if (fscanf(in,"%*s%s",buf) != 1)
      fatal_error(0,"Not enough arguments read line %d",__LINE__);
    for(j=0; (j<nj); j++) {
      if (full) {
	if (fscanf(in,"%d%d",&ja,&grp) != 2)
	  fatal_error(0,"Not enough arguments read line %d",__LINE__);
      }
      else {
	if (fscanf(in,"%d",&ja) != 1)
	  fatal_error(0,"Not enough arguments read line %d",__LINE__);
      }
      nbl->nl_i[i].nlj[2*j]=ja;
      nbl->nl_i[i].nlj[2*j+1]=grp;    
      nrj++;
    }
  }
  fprintf(stderr,", nrj: %d\n",nrj);
}
*/
void read_nblist(FILE *in,bool **matje)
{
  char     buf[256];
  int      i,ia,ja,nnbl,ii,jj,full;
  t_nblist nbl[3];

  for(i=0; (i<2); i++)  
    do {
      if (eof(in))
	fatal_error(0,"EOF when looking for lj-qq in logfile");
      fscanf(in,"%s",buf);
    } while ((int)strcmp(buf,"lj-qq") != 0);
    
  if (fscanf(in,"%d%d",&nnbl,&full) != 2)
    fatal_error(0,"Not enough arguments read line %d",__LINE__);
  for(i=0; (i<nnbl); i++) {
    /*low_readnblist(in,&nbl[i],full);
      for(ia=0; (ia<nbl[i].nri); ia++) {
      ii=nbl[i].nl_i[ia].i_atom;
      for(ja=0; (ja<nbl[i].nl_i[ia].nj); ja++) {
      jj=nbl[i].nl_i[ia].nlj[2*ja];
      if (ii < jj)
      matje[ii][jj]=TRUE;
      else
      matje[jj][ii]=TRUE;
      }
      }*/
  }
}

void read_nblistshift(FILE *in,int **matje,int maxatom)
{
  static   int count=0;
  FILE     *out;
  char     buf[256];
  int      t,li,ia,ja,nnbl,tempi,ii,jj,swap,gid,full,nrI,nrJ;
  t_nblist nbl[3],*list;
  t_nl_i   *nli;
  int  ntw=0;
  
  do {
    fscanf(in,"%s",buf);
  } while ((int)strcmp(buf,"lj-qq") != 0);
  if (fscanf(in,"%d%d",&nnbl,&full) != 2)
    fatal_error(0,"Not enough arguments read line %d",__LINE__);
  if (full)
    fprintf(stderr,"Also reading grp info...\n");
  for(li=0; (li<nnbl); li++) {
    list=&(nbl[li]);
    /*low_readnblist(in,list,full);
      nrI=list->nri;
      for(ia=0; (ia<nrI); ia++) {
      nli=&(list->nl_i[ia]);
      ii=nli->i_atom;
      t=nli->shift;
      if (ii < maxatom) {
      nrJ=nli->nj;
      for(ja=0; (ja<nrJ); ja++) {
      jj  = nli->nlj[2*ja];
      gid = nli->nlj[2*ja+1]+1;
      if (jj<maxatom) {
      swap=0;
      tempi=ii;
      if (tempi > jj) {
      swap=jj;
      jj=tempi;
      tempi=swap;
      swap=-1;
      }
      if (matje[tempi][jj] != 0) {
      ntw++;
      printf("entry %d-%d occurs twice, grps %d - %d\n",
      tempi,jj,matje[tempi][jj],gid);
      }
      matje[tempi][jj]=t+1;
      }
      }
      }
      }*/
  }
  if (ntw > 0)
    fprintf(stderr,"# twice=%d\n",ntw);
}

void dump_nblist(FILE *out,t_forcerec *fr,int nDNL)
{
  int  i;
  rvec *sv;
  
  fprintf(out,"Neighborlist:\n");
  fprintf(out,"%d\n",fr->nn*2);

  sv=fr->shift_vec;
  for(i=0; (i<fr->nn); i++) {
    write_nblist(out,&fr->coul[i],sv,(nDNL > 1));
    write_nblist(out,&fr->vdw[i],sv,(nDNL > 1));
  }
}

