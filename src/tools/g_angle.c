/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_g_angle_c = "$Id$";

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "physics.h"
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "rdgroup.h"
#include "macros.h"
#include "fatal.h"
#include "xvgr.h"
#include "gstat.h"

#define MAX_FRAMES 10000 /* slightly arbitrary. */

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_angle computes the angle distribution for a number of angles",
    "or dihedrals. This way you can check whether your simulation",
    "is correct. With option -ov you can plot the average angle of",
    "a group of angles as a function of time.[PAR]",
    "With the -of option g_angle also calculates the fraction of trans",
    "dihedrals (only for dihedrals) as function of time, but this is",
    "probably only fun for a selected few.[PAR]",
    "With option -oc a dihedral correlation function is calculated.[PAR]",
    "It should be noted that the indexfile should contain",
    "atom-triples for angles or atom-quadruplets for dihedrals.",
    "If this is not the case, the program will crash."
  };
  static char *opt=NULL;
  static bool bALL=FALSE,bFour=TRUE,bNumber=FALSE;
  static int  resolution=5,nrestart=1;
  static int  nframes = 10;
  t_pargs pa[] = {
    { "-type", FALSE, etSTR, &opt,
      "Select either A (angles), D (dihedrals), I (impropers), R (Ryckaert-Bellemans)" },
    { "-all",    FALSE,  etBOOL, &bALL,
      "Plot all angles separately in the averages file, in the order of appearance in the index file. This way the first graph is the average, the rest are the individual angles." },
    { "-resolution", FALSE, etINT, &resolution,
      "The number of points per degree in the distribution" },
    { "-nframes",   FALSE, etINT,  &nframes,
      "Number of frames in your trajectory" },
    { "-nrestart", FALSE, etINT, &nrestart,
      "Without FFT this is the number of points for calculation of ACF, when set to 1 all points are taken into account" },
    { "-number", FALSE,  etBOOL, &bNumber,
      "Use Chandler correlation function (N[trans] = 1, N[gauche] = 0) rather than cosine correlation function. Trans is defined as phi < -60 || phi > 60." },
    { "-fft",    FALSE,  etBOOL, &bFour,
      "Use FFT for correlation function" }
  };
  static char *bugs[] = {
    "Counting transitions only works for dihedrals with multiplicity 3"
  };
  
  FILE       *out;
  real       tmp,dt;
  int        status,isize;
  atom_id    *index;
  char       *grpname;
  real       maxang=0,Jc,S2;
  unsigned long mode;
  int        maxangstat=0,mult,*angstat;
  int        i,j,total,nangles,natoms,nat2,first,last,angind;
  bool       bAver,bRb=FALSE,
    bFrac,          /* calculate fraction too?  */
    bTrans,         /* worry about transtions too? */
    bCorr;          /* correlation function ? */    
  real       t,aa,fraction;       /* fraction trans dihedrals */
  double     tfrac;
  char       title[256];
  real       **dih;          /* mega array with all dih. angles at all times*/
  char       buf[80];       
  real       *time,*trans_frac,*aver_angle;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL,  ffREAD  },
    { efTPB, NULL, NULL,  ffREAD  },
    { efNDX, NULL, NULL,  ffREAD  },
    { efXVG, "-od", "angdist",  ffWRITE },
    { efXVG, "-ov", "angaver",  ffOPTWR },
    { efXVG, "-of", "dihfrac",  ffOPTWR },
    { efXVG, "-ot", "dihtrans", ffOPTWR },
    { efXVG, "-oh", "trhisto",  ffOPTWR },
    { efXVG, "-oc", "dihcorr",  ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,asize(bugs),bugs);
		    
  if (!opt)
    usage(argv[0],opt);
    
  mult   = 4;
  maxang = 360.0;
  switch(opt[0]) {
  case 'A':
    mult   = 3;
    maxang = 180.0;
    break;
  case 'D':
    break;
  case 'I':
    break;
  case 'R':
    bRb = TRUE;
    break;
  default:
    usage(argv[0],opt);
  }

  /* Calculate bin size */
  maxangstat=resolution*maxang;
    
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
  nangles=isize/mult;
  if ((isize % mult) != 0) 
    fatal_error(0,"number of index elements not multiple of %d, "
		"these can not be %s\n",
		mult,(mult==3) ? "angle triplets" : "dihedral quadruplets");
  

  /* Check whether specific analyses have to be performed */
  bCorr=opt2bSet("-oc",NFILE,fnm);
  bAver=opt2bSet("-ov",NFILE,fnm);
  bTrans=opt2bSet("-ot",NFILE,fnm);
  bFrac=opt2bSet("-of",NFILE,fnm);

  if (bNumber && !bCorr)
    bCorr=TRUE;
    
  if (bFrac && bRb==FALSE) {
    fprintf(stderr,"Warning:"
	    " calculating fractions as defined in this program\n"
	    "makes sense for Ryckaert Bellemans dihs. only. Ignoring -of\n\n"); 
    bFrac = FALSE;
  }
  
  if ( (bTrans || bFrac || bCorr) && mult==3)
    fatal_error(0,"Can only do transition, fraction or correlation\n"
		"on dihedrals. Select -d\n");
  
  /* 
   * We need to know the nr of frames so we can allocate memory for an array 
   * with all dihedral angles at all timesteps. Works for me.
   */
  if (bTrans || bCorr  || bALL) {
    fprintf(stderr,"Reading at most %d frames. Hang on.\n",nframes);
    
    snew(dih,nangles);
    for (i=0; (i<nangles); i++)
      snew(dih[i],nframes);
  }

  /* Allocate time array */
  snew(time,nframes);
  snew(trans_frac,nframes);
  snew(aver_angle,nframes);
  snew(angstat,maxangstat+1);

  read_ang_dih(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPB,NFILE,fnm),(mult == 3),
	       bALL || bCorr || bTrans,bRb,maxangstat,angstat,
	       &nframes,time,isize,index,trans_frac,aver_angle,dih);
	       
  dt=(time[nframes-1]-time[0])/(nframes-1);
  
  if (bAver) {
    sprintf(title,"Average Angle: %s",grpname);
    out=xvgropen(opt2fn("-ov",NFILE,fnm),
		 title,"Time (ps)","Angle (degrees)");
    for(i=0; (i<nframes); i++) {
      fprintf(out,"%10.5f  %8.3f",time[i],aver_angle[i]);
      if (bALL)
	for(j=0; (j<nangles); j++)
	  fprintf(out,"  %8.3f",dih[j][i]);
      fprintf(out,"\n");
    }	
    fclose(out);
  }
  
  if (bFrac) {
    sprintf(title,"Trans fraction: %s",grpname);
    out=xvgropen(opt2fn("-of",NFILE,fnm),
		  title,"Time (ps)","Fraction");
    tfrac = 0.0;
    for(i=0; (i<nframes); i++) {
      fprintf(out,"%10.5f  %10.3f\n",time[i],trans_frac[i]);
      tfrac += trans_frac[i];
    }
    fclose(out);
    
    tfrac/=nframes;
    fprintf(stderr,"Average trans fraction: %g\n",tfrac);
  }
  sfree(trans_frac);
  
  if (bTrans) 
    ana_dih_trans(opt2fn("-ot",NFILE,fnm),opt2fn("-oh",NFILE,fnm),
		  dih,nframes,nangles,grpname,time[0],dt,bRb);
		  
  if (bCorr) {
    /* Autocorrelation function */
    if (nframes < 2)
      fprintf(stderr,"Not enough frames for correlation function\n");
    else {
      
      if (bNumber) {
	real dval,sixty=DEG2RAD*60;
	bool bTest;

	for(i=0; (i<nangles); i++)
	  for(j=0; (j<nframes); j++) {
	    dval = dih[i][j];
	    if (bRb)
	      bTest=(dval > -sixty) && (dval < sixty);
	    else
	      bTest=(dval < -sixty) || (dval > sixty);
	    if (bTest)
	      dih[i][j] = dval-tfrac;
	    else
	      dih[i][j] = -tfrac;
	  }
      }
      if (bNumber)
	mode = eacNormal;
      else
	mode = eacCos;
      do_autocorr(opt2fn("-oc",NFILE,fnm),"Dihedral Autocorrelation Function",
		  nframes,nangles,dih,dt,mode,FALSE,NULL,NULL);
    }
  }

  
  /* Determine the non-zero part of the distribution */
  for(first=0; (first <= maxangstat) && (angstat[first] == 0); first++)
    ;
  for(last=maxangstat; (last >= 0) && (angstat[last] == 0) ; last--)
    ;

  fprintf(stderr,"Found points in the range from %d to %d (max %d)\n",
	  first,last,maxangstat);
    
  if (mult == 3)
    sprintf(title,"Angle Distribution: %s",grpname);
  else {
    sprintf(title,"Dihedral Distribution: %s",grpname);
    
    calc_distribution_props(maxangstat,angstat,-180.0,0,NULL,&S2);
    fprintf(stderr,"Order parameter S^2 = %g\n",S2);
  }
  out=xvgropen(opt2fn("-od",NFILE,fnm),title,"Degrees","");
  for(i=first; (i<=last); i++) 
    fprintf(out,"%10g  %10d\n",
	    ((i*maxang)/maxangstat)+180.0-maxang,angstat[i]);
  
  fclose(out);

  xvgr_file(opt2fn("-od",NFILE,fnm),NULL);
  if (bAver)
    xvgr_file(opt2fn("-ov",NFILE,fnm),"-nxy");
    
  thanx(stdout);
    
  return 0;
}
