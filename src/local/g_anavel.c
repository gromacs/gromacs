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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_com_c = "$Id$";

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "random.h"
#include "names.h"
#include "matio.h"
#include "physics.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "xvgr.h"
#include "string2.h"
#include "rdgroup.h"
#include "tpxio.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_anavel computes temperature profiles in a sample. The sample",
    "can be analysed radial, i.e. the temperature as a function of",
    "distance from the center, cylindrical, i.e. as a function of distance",
    "from the vector (0,0,1) through the center of the box, or otherwise",
    "(will be specified later)"
  };
  t_filenm fnm[] = {
    { efTRN,  "-f",  NULL, ffREAD },
    { efTPX,  "-s",  NULL, ffREAD },
    { efXPM,  "-o", "xcm", ffWRITE }
  };
#define NFILE asize(fnm)

  static int  mode = 0,   nlevels = 10;
  static real tmax = 300, xmax    = -1;
  t_pargs pa[] = {
    { "-mode",    FALSE, etINT,  {&mode},    "mode" },
    { "-nlevels", FALSE, etINT,  {&nlevels}, "number of levels" },
    { "-tmax",    FALSE, etREAL, {&tmax},    "max temperature in output" },
    { "-xmax",    FALSE, etREAL, {&xmax},    "max distance from center" }
  };
  
  FILE       *fp;
  int        *npts,nmax;
  int        status;
  int        i,j,idum,step,nframe=0,natoms,index;
  real       t,temp,rdum,hboxx,hboxy,scale,xnorm;
  real       **profile=NULL;
  real       *t_x=NULL,*t_y,hi=0;
  rvec       *x,*v;
  t_topology *top;
  int        d,m,n;
  matrix     box;
  atom_id    *sysindex;
  bool       bHaveV,bReadV;
  t_rgb      rgblo = { 0, 0, 1 },rgbhi = { 1, 0, 0 };
  
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,0,NULL);

  top    = read_top(ftp2fn(efTPX,NFILE,fnm));
		    
  natoms = read_first_x_v(&status,ftp2fn(efTRN,NFILE,fnm),&t,&x,&v,box);
  
  if (xmax > 0) {
    scale  = 5;
    nmax   = xmax*scale;
  }
  else {
    scale  = 5;
    nmax   = (0.5*sqrt(sqr(box[XX][XX])+sqr(box[YY][YY])))*scale; 
  }
  snew(npts,nmax+1);
  snew(t_y,nmax+1);
  for(i=0; (i<=nmax); i++) {
    npts[i] = 0;
    t_y[i]  = i/scale;
  }
  do {
    srenew(profile,++nframe);
    snew(profile[nframe-1],nmax+1);
    srenew(t_x,nframe);
    t_x[nframe-1] = t*1000;
    hboxx = box[XX][XX]/2;
    hboxy = box[YY][YY]/2;
    for(i=0; (i<natoms); i++) {
      /* determine position dependent on mode */
      switch (mode) {
      case 0:
	xnorm = sqrt(sqr(x[i][XX]-hboxx) + sqr(x[i][YY]-hboxy));
	break;
      default:
	fatal_error(0,"Unknown mode %d",mode);
      }
      index = xnorm*scale;
      if (index <= nmax) {
	temp = top->atoms.atom[i].m*iprod(v[i],v[i])/(2*BOLTZ);
	if (temp > hi)
	  hi = temp;
	npts[index]++;
	profile[nframe-1][index] += temp;
      }
    }
    for(i=0; (i<=nmax); i++) {
      if (npts[i] != 0) 
	profile[nframe-1][i] /= npts[i];
      npts[i] = 0;
    }
  } while (read_next_x_v(status,&t,natoms,x,v,box));
  close_trn(status);

  fp = ftp2FILE(efXPM,NFILE,fnm,"w");
  write_xpm(fp,"Temp. profile","T (a.u.)",
	    "t (fs)","R (nm)",
	    nframe,nmax+1,t_x,t_y,profile,0,tmax,
	    rgblo,rgbhi,&nlevels);
  
  thanx(stdout);
  
  return 0;
}

