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

real calc_gyro(rvec x[],int gnx,atom_id index[],t_atom atom[],real tm,
	       rvec gvec,rvec d,bool bQ,bool bRot,bool bMOI,matrix trans)
{
  int    i,ii,m;
  real   gyro,dx2,m0,Itot;
  rvec   comp;

  if (bRot) {
    principal_comp(gnx,index,atom,x,trans,d);
    Itot = norm(d);
    if (bMOI)
      return Itot;
    for(m=0; (m<DIM); m++)
      d[m]=sqrt(d[m]/tm);
#ifdef DEBUG
    pr_rvecs(stderr,0,"trans",trans,DIM);
#endif
    /* rotate_atoms(gnx,index,x,trans); */
  }
  clear_rvec(comp);
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    if (bQ)
      m0=fabs(atom[ii].q);
    else
      m0=atom[ii].m;
    for(m=0; (m<DIM); m++) {
      dx2=x[ii][m]*x[ii][m];
      comp[m]+=dx2*m0;
    }
  }
  gyro=comp[XX]+comp[YY]+comp[ZZ];
  
  for(m=0; (m<DIM); m++)
    gvec[m]=sqrt((gyro-comp[m])/tm);
  
  return sqrt(gyro/tm);
}

int gmx_gyrate(int argc,char *argv[])
{
  static char *desc[] = {
    "g_gyrate computes the radius of gyration of a group of atoms",
    "and the radii of gyration about the x, y and z axes,",
    "as a function of time. The atoms are explicitly mass weighted.",
    "With the [TT]-nmol[tt] option the radius of gyration will be calculated",
    "for multiple molecules by splitting the analysis group in equally",
    "sized parts."
  };
  static int  nmol=1;
  static bool bQ=FALSE,bRot=FALSE,bMOI=FALSE;
  t_pargs pa[] = {
    { "-nmol", FALSE, etINT, {&nmol},
      "The number of molecules to analyze" },
    { "-q", FALSE, etBOOL, {&bQ},
      "Use absolute value of the charge of an atom as weighting factor instead of mass" },
    { "-p", FALSE, etBOOL, {&bRot},
      "Calculate the radii of gyration about the principal axes." },
    { "-moi", FALSE, etBOOL, {&bMOI},
      "Calculate the moments of inertia (defined by  the principal axes)." }
  };
  FILE       *out;
  int        status;
  t_topology top;
  rvec       *x,*x_s;
  rvec       xcm,gvec,gvec1;
  matrix     box,trans;
  bool       bACF;
  real       **moi_trans=NULL;
  int        max_moi=0,delta_moi=100;
  rvec       d,d1;         /* eigenvalues of inertia tensor */
  real       t,t0,tm,gyro;
  int        natoms;
  char       *grpname,title[256];
  int        i,j,m,gnx,nam,mol;
  atom_id    *index;
  char       *leg[]  = { "Rg", "RgX", "RgY", "RgZ" }; 
  char       *legI[] = { "Itot", "I1", "I2", "I3" }; 
#define NLEG asize(leg) 
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD }, 
    { efTPS, NULL,   NULL,       ffREAD },
    { efNDX, NULL,   NULL,       ffOPTRD },
    { efXVG, NULL,   "gyrate",   ffWRITE }, 
    { efXVG, "-acf", "moi-acf",  ffOPTWR }
  }; 
#define NFILE asize(fnm) 
  int     npargs;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL); 
  bACF = opt2bSet("-acf",NFILE,fnm);
  if (bACF && nmol!=1)
    gmx_fatal(FARGS,"Can only do acf with nmol=1");
  bRot = bRot || bMOI || bACF;
  if (bRot) {
    printf("Will rotate system along principal axes\n"); 
    snew(moi_trans,DIM);
  }
  if (bMOI) {
    printf("Will print moments of inertia\n");
    bQ = FALSE;
  }
  if (bQ) 
    printf("Will print radius normalised by charge\n"); 
    
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,TRUE);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (nmol > gnx || gnx % nmol != 0) {
    gmx_fatal(FARGS,"The number of atoms in the group (%d) is not a multiple of nmol (%d)",gnx,nmol);
  }
  nam = gnx/nmol;

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 
  snew(x_s,natoms); 

  j  = 0; 
  t0 = t;
  if (bQ) 
    out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
		 "Radius of Charge","Time (ps)","Rg (nm)"); 
  else if (bMOI)
    out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
		 "Moments of inertia","Time (ps)","I (a.m.u. nm\\S2\\N)"); 
  else 
    out=xvgropen(ftp2fn(efXVG,NFILE,fnm), 
		 "Radius of gyration","Time (ps)","Rg (nm)"); 
  if (bMOI) 
    xvgr_legend(out,NLEG,legI);
  else {
    if (bRot)
      if (bPrintXvgrCodes())
	fprintf(out,"@ subtitle \"Axes are principal component axes\"\n");
    xvgr_legend(out,NLEG,leg);
  }
  do {
    rm_pbc(&top.idef,natoms,box,x,x_s);
    gyro = 0;
    clear_rvec(gvec);
    clear_rvec(d);
    for(mol=0; mol<nmol; mol++) {
      tm    = sub_xcm(x_s,nam,index+mol*nam,top.atoms.atom,xcm,bQ);
      gyro += calc_gyro(x_s,nam,index+mol*nam,top.atoms.atom,
			tm,gvec1,d1,bQ,bRot,bMOI,trans);
      rvec_inc(gvec,gvec1);
      rvec_inc(d,d1);
    }
    if (nmol > 0) {
      gyro /= nmol;
      svmul(1.0/nmol,gvec,gvec);
      svmul(1.0/nmol,d,d);
    }

    if (bRot) {
      if (j >= max_moi) {
	max_moi += delta_moi;
	for(m=0; (m<DIM); m++)
	  srenew(moi_trans[m],max_moi*DIM);
      }
      for(m=0; (m<DIM); m++)
	copy_rvec(trans[m],moi_trans[m]+DIM*j);
      fprintf(out,"%10g  %10g  %10g  %10g  %10g\n",
	      t,gyro,d[XX],d[YY],d[ZZ]); }
    else {
      fprintf(out,"%10g  %10g  %10g  %10g  %10g\n",
	      t,gyro,gvec[XX],gvec[YY],gvec[ZZ]); }
    j++;
  } while(read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  
  fclose(out);

  if (bACF) {
    int mode = eacVector;
  
    do_autocorr(opt2fn("-acf",NFILE,fnm),
		"Moment of inertia vector ACF",
		j,3,moi_trans,(t-t0)/j,mode,FALSE);
    do_view(opt2fn("-acf",NFILE,fnm),"-nxy");
  }
  
  do_view(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}
