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
static char *SRCID_g_diffuse_c = "$Id$";

#include <stdio.h> 
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "maths.h"
#include "futil.h"
#include "vec.h"
#include "rdgroup.h"
#include "copyrite.h"
#include "typedefs.h"
#include "xvgr.h"
#include "gstat.h"

void prepare_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		  matrix box)
{
  int  i,m,ind;
  rvec hbox;
  
  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m]=0.5*box[m][m];
  for(i=0; (i<gnx); i++) {
    ind=index[i];
    for(m=0; (m<DIM); m++) {
      while(xcur[ind][m]-xprev[ind][m] <= hbox[m])
	xcur[ind][m] += box[m][m];
      while(xcur[ind][m]-xprev[ind][m] >  hbox[m])
	xcur[ind][m] -= box[m][m];
    }      
  }
}

static real calc_msd(rvec xref[],int nx,atom_id index[],rvec x[])
{
  int  i,m,ind;
  rvec dx;
  real g;
  
  g = 0;
  for(i=0; (i<nx); i++) {
    ind = index[i];
    rvec_sub(xref[ind],x[ind],dx);
    g += iprod(dx,dx);
  }
  return g/nx;
}

static void analyse_msd(char *xvg,char *grpname,
			int nframes,int nstart,real **msd,real xtime[],
			real bfit,real efit)
{
  FILE  *fp;
  char  buf[256];
  int   i,j,nj,minj,legind;
  real  *D,*Dsq,t0,msd0,D2,DeltaD;
  char  **leg;
  
  snew(D,nstart+1);
  snew(Dsq,nstart+1);
  snew(leg,nstart+1);
  
  sprintf(buf,"Diffusion constant for %s",grpname);
  fp=xvgropen(xvg,buf,"Time (ps)","D (cm\\S2\\Ns\\S-1\\N)");
  minj = nframes;
  for(i=0; (i<nstart); i++) {
    D[i]   = 0;
    Dsq[i] = 0;
    nj     = 0;
    for(j=0; (j<nframes); j++) {
      if (j > 0) {
	msd[i][j] *= 1000.0/(6*(xtime[j] - xtime[0]));
	if (msd[i][j] == 0.0)
	  break;
      }
      if ((xtime[j] >= bfit) && (xtime[j] <= efit)) {
	D[i]   += msd[i][j];
	Dsq[i] += sqr(msd[i][j]);
	nj++;
      }
    }
    legind = (nstart > 1) ? i+1 : 0;
    if (nj > 0) {
      D[i] /= (real) nj;
      sprintf(buf,"D%d = %.2f, N=%d",i+1,D[i],nj);
      leg[legind] = strdup(buf);
    }
    else {
      leg[legind] = strdup("No data!");
    }
    minj=min(minj,j);
  }
  /* Compute average D */
  if (nstart > 1) {
    D[nstart] = 0.0;
    D2        = 0.0;
    for(i=0; (i<nstart); i++) {
      D[nstart] += D[i];
      D2 += sqr(D[i]);
    }
    D[nstart] /= nstart;
    D2        /= nstart;
    DeltaD = sqrt(D2 - sqr(D[nstart]));
    fprintf(stderr,"Dav = %.2f, RMS = %.2f\n",D[nstart],DeltaD);
    fprintf(fp,"# Dav = %.2f, RMS = %.2f nstart = %d\n",D[nstart],DeltaD,nstart);
  
    sprintf(buf,"D\\sav\\N = %.2f",D[nstart]);
    leg[0] = strdup(buf);
  }
  xvgr_legend(fp,((nstart > 1) ? nstart+1 : 1),leg);
  
  for(j=0; (j<minj); j++) {
    fprintf(fp,"%10g",xtime[j]-xtime[0]);
    if (nstart > 1) {
      msd0 = 0;
      for(i=0; (i<nstart); i++)
	msd0 += msd[i][j];
      fprintf(fp,"  %10g",msd0/nstart);
    }
    for(i=0; (i<nstart); i++)
      fprintf(fp,"  %10g",msd[i][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}

static void do_msd(char *trx,char *xvg,
		   int nframes,int nstart,real dt,
		   int nx,atom_id index[],char *grpname,
		   real bfit,real efit)
{
  real   **msd,*xtime,t,t0,msdt;
  int    *n_offs;
  int    i,n,status,iframe,istart,natoms;
  bool   bEOF;
  rvec   *x[2],**xref=NULL;
  int    cur=0;
#define prev 1-cur
  matrix box;
  
  if (nstart <= 0)
    fatal_error(0,"nstart should be > 0 (iso %d)",nstart);
  if (nframes <= 0)
    fatal_error(0,"nframes should be > 0 (iso %d)",nframes);
    
  natoms=read_first_x(&status,trx,&t0,&(x[cur]),box);
  t = t0;
  snew(x[prev],natoms);
  memcpy(x[prev],x[cur],natoms*sizeof(x[cur][0]));
  
  snew(n_offs,nstart);
  snew(msd,nstart);
  snew(xref,nstart);
  for(i=0; (i<nstart); i++) {
    n_offs[i] = -1;
    snew(msd[i],nframes);
    snew(xref[i],natoms);
  }
  snew(xtime,nframes);
  
  /* Check index */
  for(i=0; (i<nx); i++)
    if (index[i] >= natoms)
      fatal_error(0,"Index[%d] = %d should be >=0 and < %d",i,index[i],natoms);

  /* Index is OK if we got here */
  iframe = 0;
  istart = 0;
  do {
    xtime[iframe] = t;
    /* Check for new starting point */
    if (istart < nstart) {
      if ((t >= (t0+istart*dt)) && (n_offs[istart] == -1)) {
	fprintf(stderr,"  New starting point\n");
	memcpy(xref[istart],x[cur],natoms*sizeof(x[cur][0]));
	n_offs[istart]=iframe;
	istart++;
      }
    }
    for(n=0; (n<istart); n++) {
      /* Compute the MSD for all starting points */
      msdt = calc_msd(xref[n],nx,index,x[cur]);
      msd[n][iframe-n_offs[n]] = msdt;
    }
    
    /* Read new data into new array and remove PBC */
    cur  = prev;
    bEOF = !read_next_x(status,&t,natoms,x[cur],box);
    prepare_data(nx,index,x[cur],x[prev],box);
    
    iframe++;
  } while ((iframe < nframes) && !bEOF);
  close_trj(status);

  analyse_msd(xvg,grpname,iframe,istart,msd,xtime,bfit,efit);
      
  for(i=0; (i<nstart); i++) {
    sfree(xref[i]);
    sfree(msd[i]);
  }
  sfree(xref);
  sfree(msd);
  sfree(xtime);
  sfree(n_offs);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_diffuse computes the diffusion constant using the mean square",
    "displacement of atoms from",
    "their initial positions (using the Einstein relation).[PAR]",
  };
  static char *bugs[] = {
    "When the number of starting points is too high, average and standard deviation may be meaningless"
  };
  static int  nrframes = 10;
  static int  nstart   = 1;
  static real dt       = 10.0;
  static real bfit=0,efit=10000;
  t_pargs pa[] = {
    { "-nframes", FALSE, etINT, &nrframes,
      "Number of frames in your trajectory" },
    { "-nstart",  FALSE, etINT, &nstart,
      "Number of starting points" },
    { "-tstart",  FALSE, etREAL, &dt,
      "Time between starting points" },
    { "-beginfit",FALSE, etREAL, &bfit,
      "Begin fitting to a straight line at this time (ps)" },
    { "-endfit",  FALSE, etREAL, &efit,
      "End fitting here" }
  };
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD },
    { efXVG, NULL, "diff.xvg",  ffWRITE },
  };
#define NFILE asize(fnm)
  int     nx;
  atom_id *index=NULL;
  char    *grpname=NULL;
  
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&nx,&index,&grpname);

  do_msd(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efXVG,NFILE,fnm),
	 nrframes,nstart,dt,nx,index,grpname,bfit,efit);
	
  if (bDoView())
    xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
	 
  thanx(stderr);
  
  return 0;
}
