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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_g_bond_c = "$Id$";

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
#include "fatal.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "gstat.h"

void do_bonds(FILE *log,char *fn,char *outf,int gnx,atom_id index[],
	      real blen,real tol,bool bAver)
{
#define MAXTAB 1000
  FILE   *out;
  int    *btab;
  real   b0=0,b1,db=0;
  real   bond;
  t_lsq  b_one,*b_all;
  /*real   mean, mean2, sqrdev2, sigma2; 
    int    counter;*/
  rvec   *x;
  rvec   dx;
  int    status,natoms;
  matrix box;
  real   t,fac;
  int    bind,i,j,i0,i1;
  
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
  init_pbc(box,FALSE);
  if (natoms == 0) 
    fatal_error(0,"No atoms in trajectory!");
    
  j=0;
  do {
    j++; /* count frames */
    for(i=0; (i<gnx); i+=2) {
      pbc_dx(x[index[i]],x[index[i+1]],dx);
      bond   = norm(dx);
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
	bind = ((bond-b0)/db);
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
  } while (read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  
  /*
    mean = mean / counter;
    mean2 = mean2 / counter;
    sqrdev2 = (mean2 - mean*mean);
    sigma2 = sqrdev2*counter / (counter - 1);
  */
  /* For definitions see "Weet wat je meet" */
  if (bAver) {
    fprintf(stderr,"\n");
    fprintf(stderr,"Total number of samples               : %d\n",
	    npoints_lsq(&b_one));
    fprintf(stderr,"Mean                                  : %g\n",
	    aver_lsq(&b_one));
    fprintf(stderr,"Standard deviation of the distribution: %g\n",
	    sigma_lsq(&b_one));
    fprintf(stderr,"Standard deviation of the mean        : %g\n",
	    error_lsq(&b_one));
	    
    out=xvgropen(outf,"Bond Stretching Distribution",
		 "Bond Length (nm)","Time (ps)");
    
    for(i0=0;      ((i0 < MAXTAB) && (btab[i0]==0)); i0++)
      ;
    i0=max(0,i0-1);
    for(i1=MAXTAB; ((i1 > 0)      && (btab[i1]==0)); i1--)
      ;
    i1=min(MAXTAB,i1+1);
    
    if (i0 >= i1)
      fatal_error(0,"No distribution... (i0 = %d, i1 = %d)? ? ! ! ? !",i0,i1);
    
    fac=((2.0/j)/gnx);
    for(i=i0; (i<=i1); i++)
      fprintf(out,"%16.10e  %16.10e  %10d\n",b0+i*db,btab[i]*fac,btab[i]);
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_bond makes a distribution of bond lengths. If all is well a",
    "gaussian distribution should be made when using a harmonic potential.",
    "bonds are read from a single group in the index file in order i1-j1",
    "i2-j2 thru in-jn.[PAR]",
    "[TT]-tol[tt] gives the half-width of the distribution as a fraction",
    "of the bondlength ([TT]-blen[tt]). That means, for a bond of 0.2",
    "a tol of 0.1 gives a distribution from 0.18 to 0.22"
  };
  static char *bugs[] = {
    "It should be possible to get bond information from the topology."
  };
  static real blen=-1.0,tol=0.1;
  static bool bAver=TRUE;
  t_pargs pa[] = {
    { "-blen", FALSE, etREAL, &blen, 
      "Bond length. By default length of first bond" },
    { "-tol",  FALSE, etREAL, &tol, 
      "Half width of distribution as fraction of blen" },
    { "-aver", FALSE, etBOOL, &bAver,
      "Sum up distributions" }
  };
  FILE      *fp;
  char      *grpname;
  int       gnx;
  atom_id   *index;
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD  },
    { efNDX, NULL, NULL, ffREAD  },
    { efXVG, NULL, "bonds", ffWRITE },
    { efLOG, NULL, "bonds", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  if ( !even(gnx) )
    fprintf(stderr,"WARNING: odd number of atoms (%d) in group!\n",gnx);
  fprintf(stderr,"Will gather information on %d bonds\n",gnx/2);

  if (!bAver)
    fp = ftp2FILE(efLOG,NFILE,fnm,"w");
  else
    fp = NULL;
  
  do_bonds(fp,ftp2fn(efTRX,NFILE,fnm),ftp2fn(efXVG,NFILE,fnm),gnx,index,
	   blen,tol,bAver);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
    
  thanx(stdout);
  
  return 0;
}
