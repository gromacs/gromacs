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
 * Good ROcking Metal Altar for Chronical Sinners
 */
#include "sysstuff.h"
#include "string2.h"
#include "smalloc.h"
#include "futil.h"
#include "protein.h"

static int comprb(const void *a,const void *b)
{
  t_resbond *r1,*r2;

  r1=(t_resbond *)a;
  r2=(t_resbond *)b;

  return strcasecmp(r1->resname,r2->resname);
}

t_resbond *search_rb(char *key,int nrb,t_resbond rb[])
{
  t_resbond rbkey;
  
  rbkey.resname=key;
  return (t_resbond *)bsearch((void *)&rbkey,rb,nrb,sizeof(rbkey),comprb);
}

int read_b_db(char *fn,t_resbond **rb)
{
  FILE *in;
  char buf[256],ai[12],aj[12];
  int  i,nb;
  int  nrb=0;
  t_resbond *rrb=NULL;

  in=ffopen(fn,"r");

  while(TRUE) {
    if (fscanf(in,"%s",buf) != 1)
      break;
    srenew(rrb,++nrb);
    rrb[nrb-1].resname=strdup(buf);
    fscanf(in,"%d",&nb);
    rrb[nrb-1].nb=nb;
    snew(rrb[nrb-1].rbond,nb);
    for(i=0; (i<nb); i++) {
      fscanf(in,"%s%s",ai,aj);
      rrb[nrb-1].rbond[i].ai=strdup(ai);
      rrb[nrb-1].rbond[i].aj=strdup(aj);
    }
  }
  fclose(in);

  qsort(rrb,nrb,sizeof(rrb[0]),comprb);

  *rb=rrb;
  return nrb;
}

void print_b_db(FILE *out,int nb,t_resbond rb[])
{
  int i,j;

  for(i=0; (i<nb); i++) {
    fprintf(out,"%s\t%d\n",rb[i].resname,rb[i].nb);
    for(j=0; (j<rb[i].nb); j++)
      fprintf(out,"\t%s\t%s\n",rb[i].rbond[j].ai,rb[i].rbond[j].aj);
  }
}

