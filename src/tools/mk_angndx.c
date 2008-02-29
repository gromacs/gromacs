/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "statutil.h"
#include "macros.h"
#include "string2.h"
#include "futil.h"
#include "gmx_fatal.h"

static int calc_nftype(int FTYPE,t_idef *idef)
{
  int i,nf=0;

  for(i=0; (i<idef->ntypes); i++)
    if (idef->functype[i] == FTYPE) 
      nf++;

  return nf;
}

static void fill_ft_ind(int FTYPE,int ft_ind[],t_idef *idef,char *grpnames[])
{
  char buf[125];
  int  i,ft,ind=0;

  for(i=0; (i<idef->ntypes); i++) {
    ft=idef->functype[i];
    if (ft == FTYPE) {
      ft_ind[i]=ind;
      switch (FTYPE) {
      case F_G96ANGLES:
	sprintf(buf,"Theta=%.1f_%g",idef->iparams[i].harmonic.rA,
		idef->iparams[i].harmonic.krA);
	break;
      case F_ANGLES:
	sprintf(buf,"Theta=%.1f_%g",idef->iparams[i].harmonic.rA,
		idef->iparams[i].harmonic.krA);
	break;
      case F_PDIHS:
	sprintf(buf,"Phi=%.1f_%d_%g",idef->iparams[i].pdihs.phiA,
		idef->iparams[i].pdihs.mult,idef->iparams[i].pdihs.cpA);
	break;
      case F_IDIHS:
	sprintf(buf,"Xi=%.1f_%g",idef->iparams[i].harmonic.rA,
		idef->iparams[i].harmonic.krA);
	break;
      case F_RBDIHS:
	sprintf(buf,"RB-Dihs");
	break;
      default:
	gmx_fatal(FARGS,"kjdshfgkajhgajytgajtgasuyf");
      }
      grpnames[ind]=strdup(buf);
      ind++;
    }
    else
      ft_ind[i]=-1;
  }
}

static void fill_ang(int FTYPE,int fac,
		     int nr[],int *index[],int ft_ind[],t_idef *idef)
{
  int     i,j,ft,fft,nr_fac;
  t_iatom *ia;

  ia=idef->il[FTYPE].iatoms;
  for(i=0; (i<idef->il[FTYPE].nr); ) {
    ft  = idef->functype[ia[0]];
    fft = ft_ind[ia[0]];
    if (fft == -1)
      gmx_incons("Routine fill_ang");
    nr_fac=fac*nr[fft];
    for(j=0; (j<fac); j++)
      index[fft][nr_fac+j]=ia[j+1];
    nr[fft]++;
    ia += interaction_function[ft].nratoms+1;
    i  += interaction_function[ft].nratoms+1;
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "mk_angndx makes an index file for calculation of",
    "angle distributions etc. It uses a run input file ([TT].tpx[tt]) for the",
    "definitions of the angles, dihedrals etc."
  };
  static char *opt[] = { NULL, "angle", "g96-angle", "dihedral", "improper", "ryckaert-bellemans", "phi-psi", NULL };
  t_pargs pa[] = {
    { "-type", FALSE, etENUM, {opt},
      "Type of angle" }
  };
  
  FILE       *out;
  t_topology *top;
  int        i,j,nftype,nang;
  int        mult=0,FTYPE=0;
  bool       bPP;
  int        **index;
  int        *ft_ind;
  int        *nr;
  char       **grpnames;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL, ffREAD  },
    { efNDX, NULL, "angle", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  bPP  = FALSE;
  mult = 4;
  switch(opt[0][0]) {
  case 'a':
    mult=3;
    FTYPE=F_ANGLES;
    break;
  case 'g':
    mult=3;
    FTYPE=F_G96ANGLES;
    break;
  case 'd':
    FTYPE=F_PDIHS;
    break;
  case 'i':
    FTYPE=F_IDIHS;
    break;
  case 'r':
    FTYPE=F_RBDIHS;
    break;
  case 'p':
    bPP=TRUE;
    break;
  default:
    break;
  }

  top=read_top(ftp2fn(efTPX,NFILE,fnm));

  if (!bPP) {
    nftype=calc_nftype(FTYPE,&(top->idef));
    snew(grpnames,nftype);
    snew(ft_ind,top->idef.ntypes);
    fill_ft_ind(FTYPE,ft_ind,&top->idef,grpnames);
    nang=top->idef.il[FTYPE].nr;
    
    snew(nr,nftype);
    snew(index,nftype);
    for(i=0; (i<nftype); i++)
      snew(index[i],nang*mult);
    
    fill_ang(FTYPE,mult,nr,index,ft_ind,&(top->idef));

    out=ftp2FILE(efNDX,NFILE,fnm,"w");
    for(i=0; (i<nftype); i++) {
      fprintf(out,"[ %s ]\n",grpnames[i]);
      for(j=0; (j<nr[i]*mult); j++) {
	fprintf(out," %5d",index[i][j]+1);
	if ((j % 12) == 11)
	  fprintf(out,"\n");
      }
      fprintf(out,"\n");
    }
    fclose(out);
  }
  else {
    fprintf(stderr,"Sorry, maybe later...\n");
  }
  
  thanx(stderr);
  
  return 0;
}
