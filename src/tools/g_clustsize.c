/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
 */
static char *SRCID_g_rdf_c = "$Id$";
#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "statutil.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "fftgrid.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "shift_util.h"
#include "pme.h"
#include "gstat.h"
#include "matio.h"

static void clust_size(char *ndx,char *trx,char *xpm,char *ncl,char *acl, char *mcl,
		       real cut,int nskip,int nlevels)
{
  FILE    *fp,*gp,*hp;
  atom_id *index=NULL;
  int     nindex,natoms,status;
  rvec    *x=NULL,dx;
  matrix  box;
  char    *gname;
  real    t,dx2,cut2,**cs_dist=NULL,*t_x=NULL,*t_y,mid,cav;
  int     i,j,k,ai,aj,ak,ci,cj,nframe,nclust,n_x,n_y,max_size=0;
  int     *clust_index,*clust_size,max_clust_size,nav;
  t_rgb   rlo = { 1.0, 1.0, 1.0 },
		  rmid = { 1.0, 1.0, 0.0 },rhi = { 0.0, 0.0, 1.0 };
  
  fp = xvgropen(ncl,"Number of clusters","Time (ps)","N");
  gp = xvgropen(acl,"Average cluster size","Time (ps)","#molecules");
  hp = xvgropen(mcl,"Max cluster size","Time (ps)","#molecules");
  rd_index(ndx,1,&nindex,&index,&gname);
  natoms = read_first_x(&status,trx,&t,&x,box);
  snew(clust_index,nindex);
  snew(clust_size,nindex);
  cut2   = cut*cut;
  nframe = 0;
  n_x    = 0;
  snew(t_y,nindex);
  for(i=0; (i<nindex); i++) 
    t_y[i] = i+1;
  do {
    if ((nskip == 0) || ((nskip > 0) && ((nframe % nskip) == 0))) {
      init_pbc(box,FALSE);
      max_clust_size = 1;
      for(i=0; (i<nindex); i++) {
	clust_index[i] = i;
	clust_size[i]  = 1;
      }
      
      for(i=0; (i<nindex); i++) {
	ai = index[i];
	ci = clust_index[i];
	for(j=i+1; (j<nindex); j++) {
	  cj = clust_index[j];
	  if (ci != cj) {
	    aj = index[j];
	    pbc_dx(x[ai],x[aj],dx);
	    dx2 = iprod(dx,dx);
	    if (dx2 < cut2) {
	      /* Merge clusters */
	      for(k=j; (k<nindex); k++) {
		if (clust_index[k] == cj) {
		  clust_size[cj]--;
		  clust_index[k] = ci;
		  clust_size[i]++;
		}
	      }
	    }
	  }
	}
      }
      n_x++;
      srenew(t_x,n_x);
      t_x[n_x-1] = t;
      srenew(cs_dist,n_x);
      snew(cs_dist[n_x-1],nindex);
      nclust = 0;
      cav    = 0;
      nav    = 0;
      for(i=0; (i<nindex); i++) {
	ci = clust_size[i];
	if (ci > max_clust_size) max_clust_size = ci;
	if (ci > 0) {
	  nclust++;
	  cs_dist[n_x-1][ci-1] += 100.0*ci/(real)nindex;
	  max_size = max(max_size,ci);
	  if (ci > 1) {
	    cav += ci;
	    nav++;
	  }
	}
      }
      fprintf(fp,"%10.3e  %10d\n",t,nclust);
      if (nav > 0)
	fprintf(gp,"%10.3e  %10.3f\n",t,cav/nav);
      fprintf(hp, "%10.3e  %10d\n", t, max_clust_size);
    }
    nframe++;
  } while (read_next_x(status,&t,natoms,x,box));
  close_trx(status);
  fclose(fp);
  fclose(gp);
  fclose(hp);

  /* Look for the smallest entry that is not zero */
  mid = 100.0;
  for(i=0; (i<n_x); i++)
    for(j=0; (j<max_size); j++) 
      if ((cs_dist[i][j] > 0) && (cs_dist[i][j] < mid))
	mid = cs_dist[i][j];
      
  fp = ffopen(xpm,"w");
  write_xpm3(fp,"Cluster size distribution","Fraction","Time (ps)","Size",
	     n_x,max_size-1,t_x,t_y,cs_dist,0,mid,100.0,
	     rlo,rmid,rhi,&nlevels);
  fclose(fp);
	        
  sfree(clust_index);
  sfree(clust_size);
  sfree(index);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "This program computes the size distributions of molecular/atomic clusters in",
    "the gas phase. The output is given in the form of a XPM file.",
    "The total number of clusters is written to a XVG file."
  };
  static real cutoff   = 0.35;
  static int  nskip    = 0;
  static int  nlevels  = 20;
  t_pargs pa[] = {
    { "-cut",      FALSE, etREAL, {&cutoff},
      "Largest distance (nm) to be considered in a cluster" },
    { "-nskip",    FALSE, etINT,  {&nskip},
      "Number of frames to skip between writing" },
    { "-nlevels",    FALSE, etINT,  {&nlevels},
      "Number of levels of grey in xpm output" }
  };
#define NPA asize(pa)
  char       *fnTPS,*fnNDX;
  bool       bSQ,bRDF;
  
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,    ffREAD  },
    { efNDX, NULL,  NULL,    ffOPTRD },
    { efXPM, "-o", "csize",  ffWRITE },
    { efXVG, "-nc","nclust", ffWRITE },
    { efXVG, "-mc","maxclust", ffWRITE },
    { efXVG, "-ac","avclust", ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  
  clust_size(fnNDX,ftp2fn(efTRX,NFILE,fnm),ftp2fn(efXPM,NFILE,fnm),
	     opt2fn("-nc",NFILE,fnm),opt2fn("-ac",NFILE,fnm),opt2fn("-mc",NFILE,fnm),
	     cutoff,nskip,nlevels);

  thanx(stderr);
  
  return 0;
}
