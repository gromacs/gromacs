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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
#include "stdio.h"
#include "smalloc.h"
#include "copyrite.h"
#include "statutil.h"
#include "macros.h"
#include "rdgroup.h"

#define MAXLINE 10000

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "merge_freq reads a frequencey file produced by g_hbond",
    "and sums some of the columns."
  };

  FILE    *fp; /* file pointers */
  char    line[MAXLINE];    /* read frame */
  char    **matrix,**mat2;
  int     *n;       /* binary 1 0 array */
  real    *ac;      /* autocorrelation function */
  int     i,j,k;       /* just counters */
  int     t,tau;     /* time and tau counter */
  int     nbonds;    /* the number of hydrogen bonds (cols ) */
  int     nframes;   /* the number of time frames ( rows ) */
  int     nbout;
  int     bond,nlast;
  bool    **bBits;
  int     nres,*rsize,hsize;
  atom_id **rindex=NULL,*hindex=NULL;
  char    **rgrpname,*hgrpname;
  t_block *bl;
  real t2,dt=1.0;
  t_filenm fnm[] = {
    efDAT, "-f", "freq",  FALSE,
    efDAT, "-o", "freq2", FALSE,
    efNDX, "-nh","hb",    FALSE,
    efNDX, "-nr","res",   FALSE
  };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,asize(desc),desc,0,NULL);
  
  fp=opt2FILE("-f",NFILE,fnm,"r");
  fscanf(fp,"%d%d",&nbonds,&nframes);
  snew(matrix,nframes);
  for(i=0; (i<nframes); i++) {
    snew(matrix[i],nbonds+2);
    fscanf(fp,"%s",matrix[i]);
  }
  fclose(fp);
  
  fprintf(stderr,"How many residues do you want ? ");
  scanf("%d",&nres);
  snew(rsize,nres);
  snew(rgrpname,nres);
  snew(rindex,nres);
  fprintf(stderr,"Select your favorite residues or whatever: \n");
  rd_index(opt2fn("-nr",NFILE,fnm),nres,rsize,rindex,rgrpname);
  fprintf(stderr,"Select the hbond list: \n");
  rd_index(opt2fn("-nh",NFILE,fnm),1,&hsize,&hindex,&hgrpname);
  
  snew(bBits,nres);
  for(i=0; (i<nres); i++)
    snew(bBits[i],nbonds);
  for(i=0; (i<hsize); i+=3) {
    for(k=0; (k<nres); k++) {
      for(j=0; (j<rsize[k]); j++) {
	if ((rindex[k][j]==hindex[i])   ||
	    (rindex[k][j]==hindex[i+1]) ||
	    (rindex[k][j]==hindex[i+2])) {
	  bBits[k][i/3]=TRUE;
	}
      }
    }
  }
  fp=opt2FILE("-o",NFILE,fnm,"w");
  fprintf(fp,"%d  %d\n",nres,nframes);
  for(i=0; (i<nframes); i++) {
    for(k=0; (k<nres); k++) {
      line[k]='0';
      for(j=0; (j<nbonds); j++)
	if (bBits[k][j])
	  if (matrix[i][j] == '1')
	    line[k]='1';
    }
    line[k]='\0';
    fprintf(fp,"%s\n",line);
  }
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
