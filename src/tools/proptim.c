/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Grunge ROck MAChoS
 */
static char *SRCID_proptim_c = "$Id$";

#include "typedefs.h"
#include "maths.h"
#include "string2.h"
#include "hxprops.h"
#include "fatal.h"
#include "futil.h"
#include "smalloc.h"
#include "readev.h"
#include "macros.h"
#include "confio.h"
#include "copyrite.h"
#include "statutil.h"
#include "orise.h"
#include "pinput.h"
#include "recomb.h"

void calc_prj(int natoms,rvec xav[],rvec x[],rvec ev[],real eprj)
{
  int  i,m;
  
  for(i=0; (i<natoms); i++) {
    for(m=0; (m<DIM); m++)
      x[i][m]=xav[i][m]+ev[i][m]*eprj;
  }
}

void mkptrj(char *prop,int nSel,
	    int natoms,rvec xav[],int nframes,
	    int nev,rvec **EV,real **evprj,
	    int nca,atom_id ca_index[],atom_id bb_index[],
	    t_atom atom[],matrix box)
{
  FILE       *out;
  char       buf[256];
  double     propje,*pav,*pav2;
  int        i,j;
  rvec       *xxx;
  
  snew(pav,nev);
  snew(pav2,nev);
  snew(xxx,natoms);

  sprintf(buf,"%s.trj1",prop);
  out=ffopen(buf,"w");
  fprintf(out,"Projection of %s on EigenVectors\n",prop);
  for(j=0; (j<nframes); j++) {
    if ((j % 10) == 0)
      fprintf(stderr,"\rFrame %d",j);
    for(i=0; (i<nev); i++) {
      calc_prj(natoms,xav,xxx,EV[i],evprj[i][j]);
      switch (nSel) {
      case efhRAD:
	propje=radius(NULL,nca,ca_index,xxx);
	break;
      case efhLEN:
	propje=ahx_len(nca,ca_index,xxx,box);
	break;
      case efhTWIST:
	propje=twist(NULL,nca,ca_index,xxx);
	break;
      case efhCPHI:
	propje=ca_phi(nca,ca_index,xxx);
	break;
      case efhDIP:
	propje=dip(natoms,bb_index,xxx,atom);
	break;
      default:
	fatal_error(0,"Not implemented");
      }
      pav[i]+=propje;
      pav2[i]+=propje*propje;
      fprintf(out,"%8.3f",propje);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  fprintf(stderr,"\n");
  for(i=0; (i<nev); i++) {
    printf("ev %2d, average: %8.3f  rms: %8.3f\n",
	   i+1,pav[i]/nframes,sqrt(pav2[i]/nframes-sqr(pav[i]/nframes)));
  }
  sfree(pav);
  sfree(pav2);
  sfree(xxx);
}

void proptrj(char *fngro,char *fndat,t_topology *top,t_pinp *p)
{
  static char *ppp[efhNR] = { 
    "RAD", "TWIST", "RISE", "LEN", "NHX", "DIP", "RMS", "CPHI", 
    "RMSA", "PHI", "PSI", "HB3", "HB4", "HB5" 
  };
  FILE       *out;
  rvec       **EV;
  real       **evprj;
  atom_id    *index;
  int        natoms,nca,nSel,nframes,nev;
  rvec       *xav,*vav;
  atom_id    *ca_index,*bb_index;
  matrix     box;
  char       buf[256],*prop;
  double     x;
  int        i,j,d;
  
  nframes = p->nframes;
  nSel    = p->nSel;
  nev     = p->nev;
  
  prop=ppp[nSel];
  evprj=read_proj(nev,nframes,p->base);
  
  get_coordnum(fngro,&natoms);
  snew(xav,natoms);
  snew(vav,natoms);
  read_conf(fngro,buf,&natoms,xav,vav,box);
  fprintf(stderr,"Succesfully read average positions (%s)\n",buf);
  
  EV=read_ev(fndat,natoms);
  
  fprintf(stderr,"Succesfully read eigenvectors\n");

  snew(index,nev);
  for(i=0; (i<nev); i++)
    index[i]=i;
  snew(bb_index,natoms);
  for(i=0; (i<natoms); i++)
    bb_index[i]=i;
  snew(ca_index,natoms);
  for(i=nca=0; (i<natoms); i++)
    if ((strcmp("CA",*(top->atoms.atomname[i])) == 0))
      ca_index[nca++]=i;

  switch (p->funct) {
  case ptMC:
    switch(nSel) {
    case efhRAD:
      optim_radius(nev,xav,EV,evprj,natoms,nca,ca_index,p);
      break;
    case efhRISE:
      optim_rise(nev,xav,EV,evprj,natoms,nca,ca_index,p);
      break;
    default:
      break;
    }
    break;
  case ptREC:
    recombine(p->recomb,p->gamma,p->nskip,nframes,nev,natoms,
	      EV,evprj,xav,bb_index);
    break;
  case ptPTRJ:
    mkptrj(prop,nSel,natoms,xav,nframes,
	   nev,EV,evprj,nca,ca_index,bb_index,
	   top->atoms.atom,box);
    break;
  default:
    fatal_error(0,"I Don't Know What to Do");
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "proptrj"
  };
  t_manual man = { asize(desc),desc,0,NULL,NULL,0,NULL};

  t_filenm  fnm[] = {
    { efGRO, "-c", "aver",FALSE },
    { efDAT, "-d", "eigenvec", FALSE },
    { efTPB, NULL, NULL, FALSE },
    { efDAT, "-pi","pinp", FALSE },
    { efDAT, "-po","poutp", FALSE }
  };
#define NFILE asize(fnm)
  t_topology *top;
  t_pinp     *p;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,
		    NFILE,fnm,TRUE,&man);
		      
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
  init_debug("proptim.dbg",0);
  snew(p,1);
  read_inp(opt2fn("-pi",NFILE,fnm),opt2fn("-po",NFILE,fnm),p);
  
  proptrj(ftp2fn(efGRO,NFILE,fnm),ftp2fn(efDAT,NFILE,fnm),top,p);    
  
  thanx(stdout);
  
  return 0;
}

