/*
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
#include "gmx_ana.h"


int gmx_densmap(int argc,char *argv[])
{
  const char *desc[] = {
    "g_densmap computes 2D number-density maps.",
    "It can make planar and axial-radial density maps.",
    "The output [TT].xpm[tt] file can be visualized with for instance xv",
    "and can be converted to postscript with xpm2ps.",
	"Optionally, output can be in text form to a .dat file.",
    "[PAR]",
    "The default analysis is a 2-D number-density map for a selected",
    "group of atoms in the x-y plane.",
    "The averaging direction can be changed with the option [TT]-aver[tt].",
    "When [TT]-xmin[tt] and/or [TT]-xmax[tt] are set only atoms that are",
    "within the limit(s) in the averaging direction are taken into account.",
    "The grid spacing is set with the option [TT]-bin[tt].",
    "When [TT]-n1[tt] or [TT]-n2[tt] is non-zero, the grid",
    "size is set by this option.",
    "Box size fluctuations are properly taken into account.",
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
    "The normalization of the output is set with the [TT]-unit[tt] option.",
    "The default produces a true number density. Unit [TT]nm-2[tt] leaves out",
    "the normalization for the averaging or the angular direction.",
    "Option [TT]count[tt] produces the count for each grid cell.",
    "When you do not want the scale in the output to go",
    "from zero to the maximum density, you can set the maximum",
    "with the option [TT]-dmax[tt]."
  };
  static int n1=0,n2=0;
  static real xmin=-1,xmax=-1,bin=0.02,dmin=0,dmax=0,amax=0,rmax=0;
  static gmx_bool bMirror=FALSE, bSums=FALSE;
  static const char *eaver[]={ NULL, "z", "y", "x", NULL };
  static const char *eunit[]={ NULL, "nm-3", "nm-2", "count", NULL };

  t_pargs pa[] = {
    { "-bin", FALSE, etREAL, {&bin},
      "Grid size (nm)" },
    { "-aver", FALSE, etENUM, {eaver},
      "The direction to average over" },
    { "-xmin", FALSE, etREAL, {&xmin},
      "Minimum coordinate for averaging" },
    { "-xmax", FALSE, etREAL, {&xmax},
      "Maximum coordinate for averaging" },
    { "-n1", FALSE, etINT, {&n1},
      "Number of grid cells in the first direction" },
    { "-n2", FALSE, etINT, {&n2},
      "Number of grid cells in the second direction" },
    { "-amax", FALSE, etREAL, {&amax},
      "Maximum axial distance from the center"},
    { "-rmax", FALSE, etREAL, {&rmax},
      "Maximum radial distance" },
    { "-mirror", FALSE, etBOOL, {&bMirror},
      "Add the mirror image below the axial axis" },
    { "-sums", FALSE, etBOOL, {&bSums},
      "Print density sums (1D map) to stdout" },
    { "-unit", FALSE, etENUM, {eunit},
      "Unit for the output" },
    { "-dmin", FALSE, etREAL, {&dmin},
      "Minimum density in output"},
    { "-dmax", FALSE, etREAL, {&dmax},
      "Maximum density in output (0 means calculate it)"},
  };
  gmx_bool       bXmin,bXmax,bRadial;
  FILE       *fp;
  t_trxstatus *status;
  t_topology top;
  int        ePBC=-1;
  rvec       *x,xcom[2],direction,center,dx;
  matrix     box;
  real       t,m,mtot;
  t_pbc      pbc;
  int        cav=0,c1=0,c2=0,natoms;
  char       **grpname,title[256],buf[STRLEN];
  const char *unit;
  int        i,j,k,l,ngrps,anagrp,*gnx=NULL,nindex,nradial=0,nfr,nmpower;
  atom_id    **ind=NULL,*index;
  real       **grid,maxgrid,m1,m2,box1,box2,*tickx,*tickz,invcellvol;
  real       invspa=0,invspz=0,axial,r,vol_old,vol,rowsum;
  int        nlev=51;
  t_rgb rlo={1,1,1}, rhi={0,0,0};
  output_env_t oenv;
  const char *label[]={ "x (nm)", "y (nm)", "z (nm)" };
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD }, 
    { efTPS, NULL,   NULL,       ffOPTRD }, 
    { efNDX, NULL,   NULL,       ffOPTRD }, 
    { efDAT, "-od",  "densmap",   ffOPTWR }, 
    { efXPM, "-o",   "densmap",   ffWRITE } 
  }; 
#define NFILE asize(fnm) 
  int     npargs;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL,&oenv); 
   
  bXmin = opt2parg_bSet("-xmin",npargs,pa);
  bXmax = opt2parg_bSet("-xmax",npargs,pa);
  bRadial = (amax>0 || rmax>0);
  if (bRadial) {
    if (amax<=0 || rmax<=0)
      gmx_fatal(FARGS,"Both amax and rmax should be larger than zero");
  }

  if (strcmp(eunit[0],"nm-3") == 0) {
    nmpower = -3;
    unit = "(nm^-3)";
  } else if (strcmp(eunit[0],"nm-2") == 0) {
    nmpower = -2;
    unit = "(nm^-2)";
  } else {
    nmpower = 0;
    unit = "count";
  }
  
  if (ftp2bSet(efTPS,NFILE,fnm) || !ftp2bSet(efNDX,NFILE,fnm))
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x,NULL,box,
		  bRadial);
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
  
  switch (eaver[0][0]) {
  case 'x': cav = XX; c1 = YY; c2 = ZZ; break;
  case 'y': cav = YY; c1 = XX; c2 = ZZ; break;
  case 'z': cav = ZZ; c1 = XX; c2 = YY; break;
  }

  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 

  if (!bRadial) {
    if (n1 == 0)
      n1 = (int)(box[c1][c1]/bin + 0.5);
    if (n2 == 0)
      n2 = (int)(box[c2][c2]/bin + 0.5);
  } else {
    n1 = (int)(2*amax/bin + 0.5);
    nradial = (int)(rmax/bin + 0.5);
    invspa = n1/(2*amax);
    invspz = nradial/rmax;
    if (bMirror)
      n2 = 2*nradial;
    else
      n2 = nradial;
  }
  
  snew(grid,n1);
  for(i=0; i<n1; i++)
    snew(grid[i],n2);

  box1 = 0;
  box2 = 0;
  nfr = 0;
  do {
    if (!bRadial) {
      box1 += box[c1][c1];
      box2 += box[c2][c2];
      invcellvol = n1*n2;
      if (nmpower == -3)
	invcellvol /= det(box);
      else if (nmpower == -2)
	invcellvol /= box[c1][c1]*box[c2][c2];
      for(i=0; i<nindex; i++) {
	j = index[i];
	if ((!bXmin || x[j][cav] >= xmin) &&
	    (!bXmax || x[j][cav] <= xmax)) {
	  m1 = x[j][c1]/box[c1][c1];
	  if (m1 >= 1)
	    m1 -= 1;
	  if (m1 < 0)
	    m1 += 1;
	  m2 = x[j][c2]/box[c2][c2];
	  if (m2 >= 1)
	    m2 -= 1;
	  if (m2 < 0)
	    m2 += 1;
	  grid[(int)(m1*n1)][(int)(m2*n2)] += invcellvol;
	}
      }
    } else {
      set_pbc(&pbc,ePBC,box);
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
  } while(read_next_x(oenv,status,&t,natoms,x,box));
  close_trj(status);

  /* normalize gridpoints */
  maxgrid = 0;
  if (!bRadial) {
    for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) {
	grid[i][j] /= nfr;
	if (grid[i][j] > maxgrid)
	  maxgrid = grid[i][j];
      }
    }
  } else {
    for (i=0; i<n1; i++) {
      vol_old = 0;
      for (j=0; j<nradial; j++) {
	switch (nmpower) {
	case -3: vol = M_PI*(j+1)*(j+1)/(invspz*invspz*invspa); break;
	case -2: vol =            (j+1)/(invspz*invspa);        break;
	default: vol =             j+1;                         break;
	}
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
  fprintf(stdout,"\n  The maximum density is %f %s\n",maxgrid,unit);
  if (dmax > 0)
    maxgrid = dmax;

  snew(tickx,n1+1);
  snew(tickz,n2+1);
  if (!bRadial) {
    /* normalize box-axes */
    box1 /= nfr;
    box2 /= nfr;
    for (i=0; i<=n1; i++)
      tickx[i] = i*box1/n1;
    for (i=0; i<=n2; i++)
      tickz[i] = i*box2/n2;
  } else {
    for (i=0; i<=n1; i++)
      tickx[i] = i/invspa - amax;
    if (bMirror) {
      for (i=0; i<=n2; i++)
	tickz[i] = i/invspz - rmax;
    } else {
      for (i=0; i<=n2; i++)
	tickz[i] = i/invspz;
    }
  }
  
  if (bSums)
  {
    for (i=0;i<n1;++i)
    {
	fprintf(stdout,"Density sums:\n");
    rowsum=0;
    for (j=0;j<n2;++j)
	  rowsum+=grid[i][j];
	fprintf(stdout,"%g\t",rowsum);
    }
	fprintf(stdout,"\n");
  }
  
  sprintf(buf,"%s number density",grpname[anagrp]);
  if (!bRadial && (bXmin || bXmax)) {
    if (!bXmax)
      sprintf(buf+strlen(buf),", %c > %g nm",eaver[0][0],xmin);
    else if (!bXmin)
      sprintf(buf+strlen(buf),", %c < %g nm",eaver[0][0],xmax);
    else
      sprintf(buf+strlen(buf),", %c: %g - %g nm",eaver[0][0],xmin,xmax);
  }
  if (ftp2bSet(efDAT,NFILE,fnm))
  {
  fp = ffopen(ftp2fn(efDAT,NFILE,fnm),"w");
  /*optional text form output:  first row is tickz; first col is tickx */
  fprintf(fp,"0\t");
  for(j=0;j<n2;++j)
	fprintf(fp,"%g\t",tickz[j]);
  fprintf(fp,"\n");
  
  for (i=0;i<n1;++i)
  {
    fprintf(fp,"%g\t",tickx[i]);
    for (j=0;j<n2;++j)
	  fprintf(fp,"%g\t",grid[i][j]);
	fprintf(fp,"\n");
  }
  ffclose(fp);
  }
  else
  {
  fp = ffopen(ftp2fn(efXPM,NFILE,fnm),"w");
  write_xpm(fp,MAT_SPATIAL_X | MAT_SPATIAL_Y,buf,unit,
	    bRadial ? "axial (nm)" : label[c1],bRadial ? "r (nm)" : label[c2],
	    n1,n2,tickx,tickz,grid,dmin,maxgrid,rlo,rhi,&nlev);     
  ffclose(fp);
  }
  
  thanx(stderr);

  do_view(oenv,opt2fn("-o",NFILE,fnm),NULL);

  return 0;
}
