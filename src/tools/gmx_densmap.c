/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "mshift.h"
#include "xvgr.h"
#include "princ.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "tpxio.h"
#include "gstat.h"
#include "matio.h"
#include "pbc.h"

int gmx_densmap(int argc,char *argv[])
{
  static char *desc[] = {
    "g_densmap computes 2D number-density maps.",
    "It can make planar and axial-radial density maps.",
    "The output [TT].xpm[tt] file can be visualized with for instance xv",
    "and can be converted to postscript with xpm2ps.",
    "[PAR]",
    "The default analysis is a 2-D number-density map for a selected",
    "group of atoms in the x-z plane. The grid spacing is set with the option",
    "[TT]-bin[tt]. When [TT]-nx[tt] or [TT]-nz[tt] is non-zero, the grid",
    "size is set by this option. Box size fluctuations are properly taken",
    "into account.",
    "[PAR]",
    "When options [TT]-amax[tt] and [TT]-rmax[tt] are set, an axial-radial",
    "number-density map is made. Three groups should be supplied, the centers",
    "of mass of the first two groups define the axis, the third defines the",
    "analysis group. The axial direction goes from -amax to +amax, where",
    "the center is defined as the midpoint between the centers of mass and",
    "the positive direction goes from the first to the second center of mass.",
    "The radial direction goes from 0 to rmax or from -rmax to +rmax",
    "when the [TT]-mirror[tt] option has been set.",
    "[PAR]",
    "When you do not want the scale in the output to go",
    "from zero to the maximum density, you can set the maximum",
    "with the option [TT]-dmax[tt]."
  };
  static int nx=0,nz=0;
  static real bin=0.02,dmax=0,amax=0,rmax=0;
  static bool bMirror=FALSE;

  t_pargs pa[] = {
    { "-bin", FALSE, etREAL, {&bin},
      "Grid size" },
    { "-nx", FALSE, etINT, {&nx},
      "Number of grid cells in x direction" },
    { "-nz", FALSE, etINT, {&nz},
      "Number of grid cells in z direction" },
    { "-amax", FALSE, etREAL, {&amax},
      "Maximum axial distance from the center"},
    { "-rmax", FALSE, etREAL, {&rmax},
      "Maximum radial distance" },
    { "-mirror", FALSE, etBOOL, {&bMirror},
      "Add the mirror image below the axial axis" },
    { "-dmax", FALSE, etREAL, {&dmax},
      "Maximum density (0 means calculate it)"},
  };
  bool       bRadial;
  FILE       *fp;
  int        status;
  t_topology top;
  rvec       *x,xcom[2],direction,center,dx;
  matrix     box;
  real       t,m,mtot;
  t_pbc      pbc;
  int        natoms;
  char       **grpname,title[256],buf[STRLEN];
  int        i,j,k,l,ngrps,anagrp,*gnx=NULL,nindex,nradial=0,nfr;
  atom_id    **ind=NULL,*index;
  real       **grid,maxgrid,mx,mz,boxx,boxz,*tickx,*tickz,invcellvol;
  real       invspa=0,invspz=0,axial,r,vol_old,vol;
  int        nlev=51;
  t_rgb rlo={1,1,1}, rhi={0,0,0};
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD }, 
    { efTPS, NULL,   NULL,       ffOPTRD }, 
    { efNDX, NULL,   NULL,       ffOPTRD }, 
    { efXPM, "-o",   "densmap",   ffWRITE } 
  }; 
#define NFILE asize(fnm) 
  int     npargs;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL); 
   
  bRadial = (amax>0 || rmax>0);
  if (bRadial) {
    if (amax<=0 || rmax<=0)
      gmx_fatal(FARGS,"Both amax and rmax should be larger than zero");
  }

  if (ftp2bSet(efTPS,NFILE,fnm) || !ftp2bSet(efNDX,NFILE,fnm))
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,bRadial);
  if (!bRadial) {
    ngrps = 1;
    fprintf(stderr,"\nSelect an analysis group\n");
  } else {
    ngrps = 3;
    fprintf(stderr,
	    "\nSelect two groups to define the axis and an analysis group\n");
  }
  snew(gnx,ngrps);
  snew(grpname,ngrps);
  snew(ind,ngrps);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,gnx,ind,grpname);
  anagrp = ngrps - 1;
  nindex = gnx[anagrp];
  index = ind[anagrp];
  if (bRadial) {
    if ((gnx[0]>1 || gnx[1]>1) && !ftp2bSet(efTPS,NFILE,fnm))
      gmx_fatal(FARGS,"No run input file was supplied (option -s), this is required for the center of mass calculation");
  }
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 

  if (!bRadial) {
    if (nx == 0)
      nx = (int)(box[XX][XX]/bin + 0.5);
    if (nz == 0)
      nz = (int)(box[ZZ][ZZ]/bin + 0.5);
  } else {
    nx = (int)(ceil(2*amax/bin));
    nradial = (int)(ceil(rmax/bin));
    invspa = nx/(2*amax);
    invspz = nradial/rmax;
    if (bMirror)
      nz = 2*nradial;
    else
      nz = nradial;
  }
  
  snew(grid,nx);
  for(i=0; i<nx; i++)
    snew(grid[i],nz);

  boxx = 0;
  boxz = 0;
  nfr = 0;
  do {
    if (!bRadial) {
      boxx += box[XX][XX];
      boxz += box[ZZ][ZZ];
      invcellvol = nx*nz/det(box);
      for(i=0; i<nindex; i++) {
	j = index[i];
	mx = x[j][XX]/box[XX][XX];
	if (mx >= 1)
	  mx -= 1;
	if (mx < 0)
	  mx += 1;
	mz = x[j][ZZ]/box[ZZ][ZZ];
	if (mz >= 1)
	  mz -= 1;
	if (mz < 0)
	  mz += 1;
	grid[(int)(mx*nx)][(int)(mz*nz)] += invcellvol;
      }
    } else {
      set_pbc(&pbc,box);
      for(i=0; i<2; i++) {
	if (gnx[i] == 1) {
	  /* One atom, just copy the coordinates */
	  copy_rvec(x[ind[i][0]],xcom[i]);
	} else {
	  /* Calculate the center of mass */
	  clear_rvec(xcom[i]);
	  mtot = 0;
	  for(j=0; j<gnx[i]; j++) {
	    k = ind[i][j];
	    m = top.atoms.atom[k].m;
	    for(l=0; l<DIM; l++)
	      xcom[i][l] += m*x[k][l];
	    mtot += m;
	  }
	  svmul(1/mtot,xcom[i],xcom[i]);
	}
      }
      pbc_dx(&pbc,xcom[1],xcom[0],direction);
      for(i=0; i<DIM; i++)
	center[i] = xcom[0][i] + 0.5*direction[i];
      unitv(direction,direction);
      for(i=0; i<nindex; i++) {
	j = index[i];
	pbc_dx(&pbc,x[j],center,dx);
	axial = iprod(dx,direction);
	r = sqrt(norm2(dx) - axial*axial);
	if (axial>=-amax && axial<amax && r<rmax) {
	  if (bMirror)
	    r += rmax;
	  grid[(int)((axial + amax)*invspa)][(int)(r*invspz)] += 1;
	}
      }
    }
    nfr++;
  } while(read_next_x(status,&t,natoms,x,box));
  close_trj(status);

  /* normalize gridpoints */
  maxgrid = 0;
  if (!bRadial) {
    for (i=0; i<nx; i++) {
      for (j=0; j<nz; j++) {
	grid[i][j] /= nfr;
	if (grid[i][j] > maxgrid)
	  maxgrid = grid[i][j];
      }
    }
  } else {
    for (i=0; i<nx; i++) {
      vol_old = 0;
      for (j=0; j<nradial; j++) {
	vol = M_PI*(j+1)*(j+1)/(invspz*invspz*invspa);
	if (bMirror)
	  k = j + nradial;
	else
	  k = j;
	grid[i][k] /= nfr*(vol - vol_old);
	if (bMirror)
	  grid[i][nradial-1-j] = grid[i][k];
	vol_old = vol;
	if (grid[i][k] > maxgrid)
	  maxgrid = grid[i][k];
      }
    }
  }
  fprintf(stdout,"\n  The maximum density is %f nm^-3\n",maxgrid);
  if (dmax > 0)
    maxgrid = dmax;

  snew(tickx,nx+1);
  snew(tickz,nz+1);
  if (!bRadial) {
    /* normalize box-axes */
    boxx /= nfr;
    boxz /= nfr;
    for (i=0; i<=nx; i++)
      tickx[i] = i*boxx/nx;
    for (i=0; i<=nz; i++)
      tickz[i] = i*boxz/nz;
  } else {
    for (i=0; i<=nx; i++)
      tickx[i] = i/invspa - amax;
    if (bMirror) {
      for (i=0; i<=nz; i++)
	tickz[i] = i/invspz - rmax;
    } else {
      for (i=0; i<=nz; i++)
	tickz[i] = i/invspz;
    }
  }
  
  sprintf(buf,"%s number density map",grpname[anagrp]);
  fp = ffopen(ftp2fn(efXPM,NFILE,fnm),"w");
  write_xpm(fp,MAT_SPATIAL_X | MAT_SPATIAL_Y,buf,"(nm^-3)",
	    bRadial ? "axial (nm)" : "x (nm)",bRadial ? "r (nm)" : "z (nm)",
	    nx,nz,tickx,tickz,grid,0,maxgrid,rlo,rhi,&nlev);     
  ffclose(fp);
  
  thanx(stderr);

  do_view(opt2fn("-o",NFILE,fnm),NULL);

  return 0;
}
