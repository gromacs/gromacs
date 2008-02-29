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
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "gstat.h"
#include "tpxio.h"

static void make_dist_leg(FILE *fp,int gnx,atom_id index[],t_atoms *atoms)
{
  char **leg;
  int  i;

  snew(leg,gnx/2);
  for(i=0; i<gnx; i+=2) {
    snew(leg[i/2],256);
    sprintf(leg[i/2],"%s %s%d - %s %s%d",
	    *(atoms->atomname[index[i]]),
	    *(atoms->resname[atoms->atom[index[i]].resnr]),
            atoms->atom[index[i]].resnr+1,
	    *(atoms->atomname[index[i+1]]),
	    *(atoms->resname[atoms->atom[index[i+1]].resnr]),
            atoms->atom[index[i+1]].resnr+1);
  }
  xvgr_legend(fp,gnx/2,leg);
  for(i=0; i<gnx/2; i++)
    sfree(leg[i]);
  sfree(leg);
}

static void do_bonds(FILE *log,char *fn,char *fbond,char *fdist,
		     int gnx,atom_id index[],
		     real blen,real tol,bool bAver,
		     t_topology *top,bool bAverDist)
{
#define MAXTAB 1000
  FILE   *out,*outd=NULL;
  int    *btab=NULL;
  real   b0=0,b1,db=0;
  real   bond,bav;
  t_lsq  b_one,*b_all=NULL;
  /*real   mean, mean2, sqrdev2, sigma2; 
    int    counter;*/
  rvec   *x;
  rvec   dx;
  int    status,natoms;
  matrix box;
  real   t,fac;
  int    bind,i,nframes,i0,i1;
  t_pbc  pbc;
  
  if (!bAver) {
    snew(b_all,gnx/2);
    for(i=0; (i<gnx/2); i++)
      init_lsq(&(b_all[i]));
  }
  else {
    init_lsq(&b_one);
    snew(btab,MAXTAB+1);
  }
  
  natoms=read_first_x(&status,fn,&t,&x,box);
  if (natoms == 0) 
    gmx_fatal(FARGS,"No atoms in trajectory!");
  
  if (fdist) {
    outd = xvgropen(fdist,bAverDist ? "Average distance" : "Distances",
		    "Time (ps)","Distance (nm)");
    if (!bAverDist) 
      make_dist_leg(outd,gnx,index,&(top->atoms));
  }
  
  nframes=0;
  do {
    set_pbc(&pbc,box);
    if (fdist)
      fprintf(outd," %8.4f",t);
    nframes++; /* count frames */
    bav = 0.0;
    for(i=0; (i<gnx); i+=2) {
      pbc_dx(&pbc,x[index[i]],x[index[i+1]],dx);
      bond   = norm(dx);
      if (bAverDist)
	bav += bond;
      else if (fdist)
	fprintf(outd," %.3f",bond);
      if (bAver) {
	add_lsq(&b_one,t,bond);
	if (db == 0) {
	  if (blen == -1) {
	    b0 = 0;
	    b1 = 0.2;
	    db = (b1-b0)/MAXTAB;
	  }
	  else {
	    b0   = (1.0-tol)*blen;
	    b1   = (1.0+tol)*blen;
	    db   = (2.0*(b1-b0))/MAXTAB;
	  }
	}
	bind = (int)((bond-b0)/db+0.5);
	if ((bind >= 0) && (bind <= MAXTAB))
	  btab[bind]++;
	else {
	  /*
	    printf("bond: %4d-%4d bond=%10.5e, dx=(%10.5e,%10.5e,%10.5e)\n",
	    index[i],index[i+1],bond,dx[XX],dx[YY],dx[ZZ]);
	  */
	}
      }
      else {
	add_lsq(&(b_all[i/2]),t,bond);
      }
    }
    if (bAverDist)
      fprintf(outd," %.5f",bav*2.0/gnx);
    if (fdist)
      fprintf(outd,"\n");
  } while (read_next_x(status,&t,natoms,x,box));
  close_trj(status);

  if (fdist)
    ffclose(outd);
  
  /*
    mean = mean / counter;
    mean2 = mean2 / counter;
    sqrdev2 = (mean2 - mean*mean);
    sigma2 = sqrdev2*counter / (counter - 1);
  */
  /* For definitions see "Weet wat je meet" */
  if (bAver) {
    printf("\n");
    printf("Total number of samples               : %d\n",
	    npoints_lsq(&b_one));
    printf("Mean                                  : %g\n",
	    aver_lsq(&b_one));
    printf("Standard deviation of the distribution: %g\n",
	    sigma_lsq(&b_one));
    printf("Standard deviation of the mean        : %g\n",
	    error_lsq(&b_one));
	    
    out=xvgropen(fbond,"Bond Stretching Distribution",
		 "Bond Length (nm)","");
    
    for(i0=0;      ((i0 < MAXTAB) && (btab[i0]==0)); i0++)
      ;
    i0=max(0,i0-1);
    for(i1=MAXTAB; ((i1 > 0)      && (btab[i1]==0)); i1--)
      ;
    i1=min(MAXTAB,i1+1);
    
    if (i0 >= i1)
      gmx_fatal(FARGS,"No distribution... (i0 = %d, i1 = %d)? ? ! ! ? !",i0,i1);
    
    fac=2.0/(nframes*gnx*db);
    for(i=i0; (i<=i1); i++)
      fprintf(out,"%8.5f  %8.5f\n",b0+i*db,btab[i]*fac);
    fclose(out);
  }
  else {
    fprintf(log,"%5s  %5s  %8s  %8s\n","i","j","b_aver","sigma");
    for(i=0; (i<gnx/2); i++) {
      fprintf(log,"%5u  %5u  %8.5f  %8.5f\n",1+index[2*i],1+index[2*i+1],
	      aver_lsq(&b_all[i]),sigma_lsq(&b_all[i]));
    }
  }
}

int gmx_bond(int argc,char *argv[])
{
  static char *desc[] = {
    "g_bond makes a distribution of bond lengths. If all is well a",
    "gaussian distribution should be made when using a harmonic potential.",
    "bonds are read from a single group in the index file in order i1-j1",
    "i2-j2 thru in-jn.[PAR]",
    "[TT]-tol[tt] gives the half-width of the distribution as a fraction",
    "of the bondlength ([TT]-blen[tt]). That means, for a bond of 0.2",
    "a tol of 0.1 gives a distribution from 0.18 to 0.22.[PAR]",
    "Option [TT]-d[tt] plots all the distances as a function of time.",
    "This requires a structure file for the atom and residue names in",
    "the output. If however the option [TT]-averdist[tt] is given (as well",
    "or separately) the average bond length is plotted instead."
  };
  static char *bugs[] = {
    "It should be possible to get bond information from the topology."
  };
  static real blen=-1.0,tol=0.1;
  static bool bAver=TRUE,bAverDist=TRUE;
  t_pargs pa[] = {
    { "-blen", FALSE, etREAL, {&blen}, 
      "Bond length. By default length of first bond" },
    { "-tol",  FALSE, etREAL, {&tol}, 
      "Half width of distribution as fraction of blen" },
    { "-aver", FALSE, etBOOL, {&bAver},
      "Average bond length distributions" },
    { "-averdist", FALSE, etBOOL, {&bAverDist},
      "Average distances (turns on -d)" }
  };
  FILE      *fp;
  char      *grpname,*fdist;
  int       gnx;
  atom_id   *index;
  char      title[STRLEN];
  t_topology top;
  rvec      *x;
  matrix    box;

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD  },
    { efNDX, NULL, NULL, ffREAD  },
    { efTPS, NULL, NULL, ffOPTRD },
    { efXVG, "-o", "bonds", ffWRITE },
    { efLOG, NULL, "bonds", ffOPTWR },
    { efXVG, "-d", "distance", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE ,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  if (bAverDist)
    fdist = opt2fn("-d",NFILE,fnm);
  else {
    fdist = opt2fn_null("-d",NFILE,fnm);
    if (fdist)
      read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,FALSE);
  }
  
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  if ( !even(gnx) )
    fprintf(stderr,"WARNING: odd number of atoms (%d) in group!\n",gnx);
  fprintf(stderr,"Will gather information on %d bonds\n",gnx/2);
  
  if (!bAver)
    fp = ftp2FILE(efLOG,NFILE,fnm,"w");
  else
    fp = NULL;
  
  do_bonds(fp,ftp2fn(efTRX,NFILE,fnm),opt2fn("-o",NFILE,fnm),fdist,gnx,index,
	   blen,tol,bAver,&top,bAverDist);
  
  do_view(opt2fn("-o",NFILE,fnm),"-nxy");
  do_view(opt2fn_null("-d",NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}
 
