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
 * Gyas ROwers Mature At Cryogenic Speed
 */
#include "sysstuff.h"
#include "smalloc.h"
#include "futil.h"
#include "protein.h"

int icomp(const void *a,const void *b)
{
  t_idihres *ra,*rb;

  ra=(t_idihres *)a;
  rb=(t_idihres *)b;
  
  return strcasecmp(ra->resname,rb->resname);
}

t_idihres *search_idih(char *key,int nrdh,t_idihres ires[])
{
  t_idihres rkey;

  rkey.resname=key;

  return (t_idihres *)bsearch((void *)&key,ires,nrdh,sizeof(rkey),icomp);
}

int read_idih_db(t_idihres **ires)
{
  FILE      *in;
  t_idihres *idhr=NULL;
  int       ndhr=0;
  int       i,j,nidih;
  char      idb[256],name[STRLEN],ai[12];

  sprintf(idb,"%s/impropers.db",getenv("GMXLIB"));
  fprintf(stderr,"Reading improper dihedral database from %s\n",idb);
  in=ffopen(idb,"r");
  while (TRUE) {
    if (fscanf(in,"%s%d",name,&nidih) != 2)
      break;
    srenew(idhr,++ndhr);
    idhr[ndhr-1].resname=strdup(name);
    idhr[ndhr-1].nidih=nidih;
    snew(idhr[ndhr-1].idih,nidih);
    for(i=0; (i<nidih); i++) {
      for(j=0; (j<4); j++) {
	fscanf(in,"%s",ai);
	idhr[ndhr-1].idih[i].ai[j]=strdup(ai);
      }
    }
  }
  fclose(in);

  qsort(idhr,ndhr,sizeof(idhr[0]),icomp);

  *ires=idhr;
  return ndhr;
}

void print_idih_db(char *fn,int nrdh,t_idihres ires[])
{
  FILE *out;
  int  i,j,k;
  
  fprintf(stderr,"Writing improper dihedral database to %s\n",fn);
  out=ffopen(fn,"w");
  for(i=0; (i<nrdh); i++) {
    fprintf(out,"%s\t%d\n",ires[i].resname,ires[i].nidih);
    for(j=0; (j<ires[i].nidih); j++) {
      for(k=0; (k<4); k++)
	fprintf(out,"\t%s",ires[i].idih[j].ai[k]);
      fprintf(out,"\n");
    }
  }
  fclose(out);
}

