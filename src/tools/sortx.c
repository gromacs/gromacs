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
 * GROtesk MACabre and Sinister
 */
#include <stdio.h>
#include <stdlib.h>
#include <typedefs.h>
#include <smalloc.h>
#include <confio.h>
#include <string2.h>

void dump_ren(char *outfile,char *title,
	      int natoms,int nmol,rvec x[],rvec v[],matrix box,
	      char *atomname[],char *resname[],
	      int renum[])
{
  FILE *out;
  int  i,n,m,r;

  if ((out=fopen(outfile,"w")) == NULL) {
    perror (outfile);
    exit (1);
  }

  /* title */
  fprintf(out,"%s\n",title);
  fprintf(out,"%5d\n",natoms);

  for(i=0; (i<nmol); i++) {
    r=3*renum[i];
    for(m=0; (m<3); m++) {
      fprintf(out,"%5d%5s%5s%5d",i+1,resname[r+m],atomname[r+m],i*3+m+1);
      for(n=0; (n<DIM); n++)
	fprintf(out,"%8.3f",x[r+m][n]);
      for(n=0; (n<DIM); n++)
	fprintf(out,"%8.4f",v[r+m][n]);
      fprintf(out,"\n");
    }
  }
  fprintf (out,"%10.5f%10.5f%10.5f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
  fclose (out);
}

void main (int argc, char *argv[])
{
  int     i,m,natoms,nmol;
  char    title[STRLEN];
  rvec    *x,*v;
  int     *renum;
  char    **atomname,**resname;
  matrix  box;
  t_block cgs;

  if (argc < 2) {
    fprintf(stderr,"Usage: %s confin [ confout ]\n",argv[0]);
    exit(1);
  }
  get_coordnum(argv[1],&natoms);
  nmol=natoms/3;
  snew(x,natoms);
  snew(renum,nmol);
  snew(v,natoms);
  snew(atomname,natoms);
  snew(resname,natoms);
  snew(cgs.index,nmol+1);
  snew(cgs.a,natoms);

  read_whole_conf(argv[1],title,&natoms,x,v,box,atomname,resname);

  if (argc > 2) {
    for(i=0; (i<nmol); i++)
      scanf("%*s%*s%*s%d",&(renum[i]));
    
    dump_ren(argv[2],title,natoms,nmol,x,v,box,atomname,resname,renum);
  }
  else {
    for(i=0; (i<nmol); i++) {
      for(m=0; (m<DIM); m++)
	printf("%8.3f",x[3*i][m]);
      printf("%5d\n",i);
    }
  }
}
