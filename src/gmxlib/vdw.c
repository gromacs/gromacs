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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "vdw.h"
#include "macros.h"

int read_vdw(char *vdwdata,t_vdw **vdw) 
{
  FILE   *fp;
  t_vdw  *vdum;
  double dd;
  char   line[STRLEN],name[6];
  int    i;

  fp=libopen(vdwdata);
  *vdw=NULL;
  for (i=0; fgets(line,STRLEN,fp); i++) {
    srenew(*vdw,i+1);
    vdum=*vdw;
    strncpy(name,line,6);
    name[5]='\0';
    if ((int)strlen(name) < 2) 
      break; 
    sscanf(name,"%s",vdum[i].atomname);
    sscanf(line+5,"%lf",&dd);
    vdum[i].distance=dd;
  }
  fclose(fp);
  
  return i;
} /*read_vdw()*/

void write_vdw(char *vdwdata,t_vdw vdw[],int nvdw)
{
  FILE *fp;
  int i;

  fp=ffopen(vdwdata,"w");
  for(i=0;(i<nvdw);i++) {
    fprintf(fp,"%-5s%8.3f\n",vdw[i].atomname,vdw[i].distance);
  }
  fclose(fp);
}

real get_vdw(int nvdw,t_vdw vdw[],char *atom)
{
  int i,j,best,mlen,len;

  best=-1;
  mlen=0;
  for(i=0; (i<nvdw); i++) {
    len=min((int) strlen(atom),(int) strlen(vdw[i].atomname));
    for(j=0; (j<len); j++) 
      if (atom[j] != vdw[i].atomname[j])
	break;
    if (j > mlen) {
      mlen=j;
      best=i;
    }
  }

  if (best == -1) {
    fprintf(stderr,"Could not find atom %s\n",atom);
    return 0.0;
  }
  else
    return vdw[best].distance;
}
