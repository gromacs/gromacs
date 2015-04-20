/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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

#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "gromacs/commandline/pargs.h"
#include "random.h"
#include "names.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "copyrite.h"
#include "gromacs/fileio/tpxio.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_anavel[tt] computes temperature profiles in a sample. The sample",
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
  int        i,j,idum,step,nframe=0,index;
  real       temp,rdum,hboxx,hboxy,scale,xnorm=0;
  real       **profile=NULL;
  real       *t_x=NULL,*t_y,hi=0;
  t_topology *top;
  int        d,m,n;
  matrix     box;
  atom_id    *sysindex;
  gmx_bool       bHaveV,bReadV;
  t_rgb      rgblo = { 0, 0, 1 },rgbhi = { 1, 0, 0 };
  int        flags = TRX_READ_X | TRX_READ_V;
  t_trxframe fr;

  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,0,NULL);

  top    = read_top(ftp2fn(efTPX,NFILE,fnm));

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
	
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
    t_x[nframe-1] = fr.time*1000;
    hboxx = box[XX][XX]/2;
    hboxy = box[YY][YY]/2;
    for(i=0; (i<fr.natoms); i++) {
      /* determine position dependent on mode */
      switch (mode) {
      case 0:
	xnorm = sqrt(sqr(fr.x[i][XX]-hboxx) + sqr(fr.x[i][YY]-hboxy));
	break;
      default:
	gmx_fatal(FARGS,"Unknown mode %d",mode);
      }
      index = xnorm*scale;
      if (index <= nmax) {
	temp = top->atoms.atom[i].m*iprod(fr.v[i],fr.v[i])/(2*BOLTZ);
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
  } while (read_next_frame(status,&fr));
  close_trx(status);

  fp = ftp2FILE(efXPM,NFILE,fnm,"w");
  write_xpm(fp,0,"Temp. profile","T (a.u.)",
	    "t (fs)","R (nm)",
	    nframe,nmax+1,t_x,t_y,profile,0,tmax,
	    rgblo,rgbhi,&nlevels);
  
  gmx_thanx(stderr);
  
  return 0;
}

