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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_recomb_c = "$Id$";

#include "recomb.h"
#include "futil.h"
#include "wgms.h"
#include "smalloc.h"

real *read_gammaf(char *fn,int nframes)
{
  FILE   *in;
  real   *gf;
  double y;
  int    i;
  
  snew(gf,nframes);
  in=ffopen(fn,"r");
  for(i=0; (i<nframes); i++) {
    fscanf(in,"%lf",&y);
    gf[i]=y;
  }
  fclose(in);
  fprintf(stderr,"Succesfully read gamma\n");
  return gf;
}

void recombine(char *base,char *gammaf,int nskip,
	       int nframes,int nev,int natoms,
	       rvec *ev[],real *evprj[],
	       rvec yav[],atom_id all_index[])
{
  static char *format=
    "Recombined projection of Gamma trj (EV %d) in Cartesian Space\n";
  FILE *out;
  rvec *xxx,*evptr;
  real *gamma;
  real prj;
  char buf[256];
  int  i,j,n;
  real gt;
  
  gamma=read_gammaf(gammaf,nframes);
  snew(xxx,natoms);
  for(n=0; (n<nev); n++) {
    sprintf(buf,"%s%d",base,n+1);
    out=ffopen(buf,"w");
    fprintf(out,format,n+1);
    fprintf(stderr,format,n+1);
    evptr=ev[n];
    
    for(j=0; (j<nframes); j++) {
      if ((j % 50) == 0)
	fprintf(stderr,"\r frame %d",j);
      if ((nskip == 0) || ((j % nskip) == 0)) {
	gt=1.0/gamma[j];
	prj=evprj[n][j];
	for(i=0; (i<natoms); i++) {
	  xxx[i][XX]=(yav[i][XX]+prj*evptr[i][XX])*gt;
	  xxx[i][YY]=(yav[i][YY]+prj*evptr[i][YY])*gt;
	  xxx[i][ZZ]=(yav[i][ZZ]+prj*evptr[i][ZZ])*gt;
	}
	write_gms_ndx(out,natoms,all_index,xxx,NULL);
      }
    }
    fclose(out);
    fprintf(stderr,"\r");
  }
  fprintf(stderr,"\n");
  sfree(xxx);
  sfree(gamma);
}
