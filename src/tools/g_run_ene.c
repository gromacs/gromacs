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
 * Great Red Owns Many ACres of Sand 
 */
#include "typedefs.h"
#include "smalloc.h"
#include "enerio.h"
#include "statutil.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "xvgr.h"

static int *select_it(int nre,char *nm[],int *nset)
{
  bool *bE;
  int  n,k,j,i;
  int  *set;
  
  printf("Select the terms you want from the following list\n");
  printf("End your selection with 0\n");
  for(k=0; (k<nre); ) {
    for(j=0; (j<4) && (k<nre); j++,k++) 
      printf(" %3d=%14s",k+1,nm[k]);
    printf("\n");
  }
  snew(bE,nre);
  do {
    scanf("%d",&n);
    if ((n>0) && (n<=nre))
      bE[n-1]=TRUE;
  } while (n != 0);

  snew(set,nre);
  for(i=(*nset)=0; (i<nre); i++)
    if (bE[i])
      set[(*nset)++]=i;
 
  sfree(bE);
  
  return set;
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_run_ene extracts energy components from an energy file.",
    "The user is prompted to interactively select the energy terms",
    "she wants."
  };
  t_manual man = { asize(desc),desc,0,NULL,NULL,0,NULL};

  FILE     *in,*out;
  t_energy *ee;
  int      teller=0,nre;
  real     t;
  int      *set,i,nset;
  char     **enm,**leg;
  t_filenm fnm[] = {
    { efENE, "-f", NULL, ffREAD },
    { efXVG, NULL, NULL, ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,
		    NFILE,fnm,TRUE,&man);
  
  in=ftp2FILE(efENE,NFILE,fnm,"r");
  rd_ener_nms(in,&nre,&enm);
  set=select_it(nre,enm,&nset);
    
  out=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Gromacs Energies","Time (ps)",
	       "E (kJ mol\\S-1\\N)");
  snew(leg,nset);
  for(i=0; (i<nset); i++)
    leg[i]=enm[set[i]];
  xvgr_legend(out,nset,leg);    
  sfree(leg);
  
  snew(ee,nre);
  while (rd_ener(in,&t,ee)) {
    if ((teller++ % 1) == 0)
      fprintf(stderr,"\rFrame: %d",teller-1);
    fprintf(out,"%10g ",t);
    for(i=0; (i<nset); i++)
      fprintf(out,"%16.10e  ",ee[set[i]].esum);
    fprintf(out,"\n");
  }
  fprintf(stderr,"\n");
  fclose(in);
  fclose(out);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
    
  thanx(stdout);
}
