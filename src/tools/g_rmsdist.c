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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_rmsdist_c = "$Id$";

#include <math.h>
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "string2.h"
#include "vec.h"
#include "macros.h"
#include "rdgroup.h"
#include "pbc.h"
#include "xvgr.h"
#include "futil.h"
#include "matio.h"

static void calc_dist(int nind,atom_id index[],
		      rvec x[],matrix box,real **d)
{
  int     i,j;
  real    *xi;
  rvec    dx;

  for(i=0; (i<nind-1); i++) {
    xi=x[index[i]];
    for(j=i+1; (j<nind); j++) {
      pbc_dx(xi,x[index[j]],dx);
      d[i][j]=norm(dx);
    }
  }
}

static void calc_dist_tot(int nind,atom_id index[], rvec x[],
			  matrix box,real **d, real **dtot, real **dtot2,
			  bool bNMR, real **dtot1_3, real **dtot1_6)
{
  int     i,j;
  real    *xi;
  real    temp, temp2, temp1_3;
  rvec    dx;

  for(i=0; (i<nind-1); i++) {
    xi=x[index[i]];
    for(j=i+1; (j<nind); j++) {
      pbc_dx(xi,x[index[j]],dx);
      temp2=dx[XX]*dx[XX]+dx[YY]*dx[YY]+dx[ZZ]*dx[ZZ];
      temp =sqrt(temp2);
      d[i][j]=temp;
      dtot[i][j]+=temp; 
      dtot2[i][j]+=temp2;
      if (bNMR) {
	temp1_3 = 1.0/(temp*temp2);
	dtot1_3[i][j] += temp1_3;
	dtot1_6[i][j] += temp1_3*temp1_3;
      }
    }
  }
}

static void calc_nmr(int nind, int nframes, real **dtot1_3, real **dtot1_6,
		     real *max1_3, real *max1_6)
{
  int     i,j;
  real    temp1_3,temp1_6;


  for(i=0; (i<nind-1); i++) {
    for(j=i+1; (j<nind); j++) {
      temp1_3 = pow(dtot1_3[i][j]/nframes,-1.0/3.0);
      temp1_6 = pow(dtot1_6[i][j]/nframes,-1.0/6.0);
      if (temp1_3 > *max1_3) *max1_3 = temp1_3;
      if (temp1_6 > *max1_6) *max1_6 = temp1_6;
      dtot1_3[i][j] = temp1_3;
      dtot1_6[i][j] = temp1_6;
      dtot1_3[j][i] = temp1_3;
      dtot1_6[j][i] = temp1_6;
    }
  }
}

static void calc_rms(int nind, int nframes, 
		     real **dtot, real **dtot2, 
		     real **rmsmat,  real *rmsmax,
		     real **rmscmat, real *rmscmax, 
		     real **meanmat, real *meanmax)
{
  int     i,j;
  real    mean, mean2, rms, rmsc;
/*  
 * N.B. dtot and dtot2 contain the total distance and the total squared
 * distance respectively, BUT they return RMS and the scaled RMS resp.
 *
 */

  *rmsmax=-1000;
  *rmscmax=-1000;
  *meanmax=-1000;

  for(i=0; (i<nind-1); i++) {
    for(j=i+1; (j<nind); j++) {
      mean =dtot[i][j]/nframes;   
      mean2=dtot2[i][j]/nframes;  
      rms=sqrt(max(0,mean2-mean*mean));
      rmsc=rms/mean;
      if (mean > *meanmax) *meanmax=mean;
      if (rms  > *rmsmax ) *rmsmax =rms;
      if (rmsc > *rmscmax) *rmscmax=rmsc;
      meanmat[i][j]=meanmat[j][i]=mean;
      rmsmat[i][j] =rmsmat[j][i] =rms;           
      rmscmat[i][j]=rmscmat[j][i]=rmsc;      
    }
  }
}

real rms_diff(int isize,real **d,real **d_r)
{
  int  i,j;
  real r,r2;
  
  r2=0.0;
  for(i=0; (i<isize-1); i++)
    for(j=i+1; (j<isize); j++) {
      r=d[i][j]-d_r[i][j];
      r2+=r*r;
    }    
  r2/=(isize*(isize-1))/2;
  
  return sqrt(r2);
}

void main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_rmsdist computes the root mean square deviation of atom distances,",
    "which has the advantage that no fit is needed like in standard RMS",
    "deviation as computed by g_rms.",
    "The reference structure is taken from the structure file.",
    "The rmsd at time t is calculated as the rms",
    "of the differences in distance between atom-pairs in the reference",
    "structure and the structure at time t.[PAR]",
    "g_rmsdist can also produce matrices of the rms distances, rms distances",
    "scaled with the mean distance and the mean distances and matrices with",
    "NMR averaged distances (1/r^3 and 1/r^6 averaging)."
  };
  
  int          natom,i,teller;
  real         t;

  t_topology   top;
  matrix       box;
  rvec         *x;
  FILE         *fp;

  int      status,isize;
  atom_id  *index;
  char     *grpname;
  real     **d_r,**d,**dtot,**dtot2,**mean,**rms,**rmsc,*resnr;
  real     **dtot1_3,**dtot1_6;
  real     rmsnow,meanmax,rmsmax,rmscmax;
  real     max1_3,max1_6;
  t_rgb    rlo,rhi;
  char     buf[255];
  bool bRMS, bScale, bMean, bNMR3, bNMR6;
  
  static int  nlevels=40;
  static real scalemax=-1.0;
  t_pargs pa[] = {
    { "-nlevels",   FALSE, etINT,  &nlevels,
      "Discretize rms in # levels" },
    { "-max",   FALSE, etREAL, &scalemax,
      "Maximum level in matrices" }
  };
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD },
    { efTPS, NULL,   NULL,       ffREAD },
    { efNDX, NULL,   NULL,       ffOPTRD },
    { efXVG, NULL,   "distrmsd",  ffWRITE },
    { efXPM, "-rms", "rmsdist",  ffOPTWR },
    { efXPM, "-scl", "rmsscale", ffOPTWR },
    { efXPM, "-mean","rmsmean",  ffOPTWR },
    { efXPM, "-nmr3","nmr3",     ffOPTWR },
    { efXPM, "-nmr6","nmr6",     ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  bRMS  =opt2bSet("-rms", NFILE,fnm);
  bScale=opt2bSet("-scl", NFILE,fnm);
  bMean =opt2bSet("-mean",NFILE,fnm);
  bNMR3 =opt2bSet("-nmr3",NFILE,fnm);
  bNMR6 =opt2bSet("-nmr6",NFILE,fnm);
  
  /* don't read mass-database as masses (and top) are not used */
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),buf,&top,&x,NULL,box,FALSE);
  
  get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&isize,&index,&grpname);
  
  snew(d,isize);
  snew(dtot,isize);
  snew(dtot2,isize);
  if (bNMR3 || bNMR6) {
    snew(dtot1_3,isize);
    snew(dtot1_6,isize);
  }
  snew(mean,isize);
  snew(rms,isize);
  snew(rmsc,isize);
  snew(d_r,isize);
  snew(resnr,isize);
  for(i=0; (i<isize); i++) {
    snew(d[i],isize);
    snew(dtot[i],isize);
    snew(dtot2[i],isize);
    if (bNMR3 || bNMR6) {
      snew(dtot1_3[i],isize);
      snew(dtot1_6[i],isize);
    }
    snew(mean[i],isize);
    snew(rms[i],isize);
    snew(rmsc[i],isize);
    snew(d_r[i],isize);
    resnr[i]=i+1;
  }

  /*set box type*/
  init_pbc(box,FALSE);
  calc_dist(isize,index,x,box,d_r);
  sfree(x);

  /*open output files*/
  fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"RMS Deviation","Time (ps)","RMSD (nm)");
  fprintf(fp,"@ subtitle \"of distances between %s atoms\"\n",grpname);
  
  /*do a first step*/
  natom=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  
  do {
    calc_dist_tot(isize,index,x,box,d,dtot,dtot2,
		  (bNMR3 || bNMR6),dtot1_3,dtot1_6);
    
    rmsnow=rms_diff(isize,d,d_r);
    fprintf(fp,"%g  %g\n",t,rmsnow);
  } while (read_next_x(status,&t,natom,x,box));
  fprintf(stderr, "\n");

  fclose(fp);
  close_trj(status);

  teller = nframes_read();
  calc_rms(isize,teller,dtot,dtot2,mean,&meanmax,rms,&rmsmax,rmsc,&rmscmax);
  fprintf(stderr,"rmsmax = %g, rmscmax = %g\n",rmsmax,rmscmax);
  
  if (scalemax > -1.0) {
    rmsmax=scalemax;
    rmscmax=scalemax;
    meanmax=scalemax;
  }
  
  if (bNMR3 || bNMR6) {
    max1_3=0;
    max1_6=0;
    calc_nmr(isize,teller,dtot1_3,dtot1_6,&max1_3,&max1_6);
  }
  
  rlo.r=1.0, rlo.g=1.0, rlo.b=1.0;
  rhi.r=0.0, rhi.g=0.0, rhi.b=0.0;

  if ( bRMS )
    write_xpm(opt2FILE("-rms",NFILE,fnm,"w"),
	      "RMS of distance","RMS (nm)","Residue Index","Residue Index",
	      isize,isize,resnr,resnr,rms,0.0,rmsmax,rlo,rhi,&nlevels);
  
  if ( bScale )
    write_xpm(opt2FILE("-scl",NFILE,fnm,"w"),
	      "Relative RMS","RMS","Residue Index","Residue Index",
	      isize,isize,resnr,resnr,rmsc,0.0,rmscmax,rlo,rhi,&nlevels);
  
  if ( bMean )
    write_xpm(opt2FILE("-mean",NFILE,fnm,"w"),
	      "Mean Distance","Distance (nm)","Residue Index","Residue Index",
	      isize,isize,resnr,resnr,mean,0.0,meanmax,rlo,rhi,&nlevels);
  
  if (bNMR3)
      write_xpm(opt2FILE("-nmr3",NFILE,fnm,"w"),"1/r^3 averaged distances",
		"Distance (nm)","Residue Index","Residue Index",
		isize,isize,resnr,resnr,dtot1_3,0.0,max1_3,rlo,rhi,&nlevels);
  if (bNMR6)
      write_xpm(opt2FILE("-nmr6",NFILE,fnm,"w"),"1/r^6 averaged distances",
		"Distance (nm)","Residue Index","Residue Index",
		isize,isize,resnr,resnr,dtot1_6,0.0,max1_6,rlo,rhi,&nlevels);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
 
  thanx(stdout);
}
