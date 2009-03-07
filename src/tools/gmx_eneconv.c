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
#include <string.h>
#include <math.h>

#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "statutil.h"
#include "disre.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "enxio.h"
#include "vec.h"

#define TIME_EXPLICIT 0
#define TIME_CONTINUE 1
#define TIME_LAST     2
#ifndef FLT_MAX
#define FLT_MAX 1e36
#endif

static int *select_it(int nre,gmx_enxnm_t *nm,int *nset)
{
  bool *bE;
  int  n,k,j,i;
  int  *set;
  bool bVerbose = TRUE;
  
  if ((getenv("VERBOSE")) != NULL)
    bVerbose = FALSE;
  
  fprintf(stderr,"Select the terms you want to scale from the following list\n");
  fprintf(stderr,"End your selection with 0\n");

  if ( bVerbose ) {
    for(k=0; (k<nre); ) {
      for(j=0; (j<4) && (k<nre); j++,k++) 
	fprintf(stderr," %3d=%14s",k+1,nm[k].name);
      fprintf(stderr,"\n");
    }
  }

  snew(bE,nre);
  do {
    if(1 != scanf("%d",&n))
    {
      gmx_fatal(FARGS,"Cannot read energy term");
    }
    if ((n>0) && (n<=nre))
      bE[n-1]=TRUE;
  } while (n != 0);

  snew(set,nre);
  for(i=(*nset)=0; (i<nre); i++)
    if (bE[i])
      set[(*nset)++]=i;
 
  sfree(bE);
  
  return set;
}

static bool same_time(real t1,real t2)
{
  const real tol=1e-5;

  return (fabs(t1-t2) < tol);
}


bool bRgt(double a,double b)
{
  double tol = 1e-6;
  
  if ( a > (b - tol*(a+b)) )
    return TRUE;
  else
    return FALSE;
}

static void sort_files(char **fnms,real *settime,int nfile)
{
    int i,j,minidx;
    real timeswap;
    char *chptr;

    for(i=0;i<nfile;i++) {
	minidx=i;
	for(j=i+1;j<nfile;j++) {
	    if(settime[j]<settime[minidx])
		minidx=j;
	}
	if(minidx!=i) {
	    timeswap=settime[i];
	    settime[i]=settime[minidx];
	    settime[minidx]=timeswap;
	    chptr=fnms[i];
	    fnms[i]=fnms[minidx];
	    fnms[minidx]=chptr;
	}
    }
}


static int scan_ene_files(char **fnms, int nfiles,
			  real *readtime, real *timestep, int *nremax)
{
  /* Check number of energy terms and start time of all files */
  int        f,i,in,nre,nremin=0,nresav=0;
  real       t1,t2;
  char       inputstring[STRLEN];
  gmx_enxnm_t *enm=NULL;
  t_enxframe *fr;
  
  snew(fr,1);

  for(f=0; f<nfiles; f++) {
    in = open_enx(fnms[f],"r");
    do_enxnms(in,&nre,&enm);

    if (f == 0) {
      nresav  = nre;
      nremin  = nre;
      *nremax = nre;
      do_enx(in,fr);
      t1 = fr->t;
      do_enx(in,fr);
      t2 = fr->t;
      *timestep=t2-t1;
      readtime[f]=t1;
      close_enx(in);
    } else {
      nremin  = min(nremin,fr->nre);
      *nremax = max(*nremax,fr->nre);
      if (nre != nresav) {
	fprintf(stderr,
		"Energy files don't match, different number of energies:\n"
		" %s: %d\n %s: %d\n",fnms[f-1],nresav,fnms[f],fr->nre);
	fprintf(stderr,
		"\nContinue conversion using only the first %d terms (n/y)?\n"
		"(you should be sure that the energy terms match)\n",nremin);
	if(NULL==fgets(inputstring,STRLEN-1,stdin))
        { 
	      gmx_fatal(FARGS,"Error reading user input");
	}
	if (inputstring[0]!='y' && inputstring[0]!='Y') {
	  fprintf(stderr,"Will not convert\n");
	  exit(0);
	}
	nresav = fr->nre;
      }
      do_enx(in,fr);
      readtime[f] = fr->t;
      close_enx(in);
    }
    fprintf(stderr,"\n");
    free_enxnms(nre,enm);
  }

  free_enxframe(fr);
  sfree(fr);

  return nremin;
}


static void edit_files(char **fnms,int nfiles,real *readtime, 
		       real *settime,int *cont_type,bool bSetTime,bool bSort)
{
  int i;
  bool ok;
  char inputstring[STRLEN],*chptr;
  
  if(bSetTime) {
    if(nfiles==1)
      fprintf(stderr,"\n\nEnter the new start time:\n\n");
    else
      fprintf(stderr,"\n\nEnter the new start time for each file.\n"
	      "There are two special options, both disables sorting:\n\n"
	      "c (continue) - The start time is taken from the end\n"
	      "of the previous file. Use it when your continuation run\n"
	      "restarts with t=0 and there is no overlap.\n\n"
	      "l (last) - The time in this file will be changed the\n"
	      "same amount as in the previous. Use it when the time in the\n"
		"new run continues from the end of the previous one,\n"
	      "since this takes possible overlap into account.\n\n");
    
    fprintf(stderr,"          File             Current start       New start\n"
	    "---------------------------------------------------------\n");
    
    for(i=0;i<nfiles;i++) {
      fprintf(stderr,"%25s   %10.3f             ",fnms[i],readtime[i]);
      ok=FALSE;
      do {
	if(NULL==fgets(inputstring,STRLEN-1,stdin))
	{
	    gmx_fatal(FARGS,"Error reading user input");
	}
	inputstring[strlen(inputstring)-1]=0;
	
	if(inputstring[0]=='c' || inputstring[0]=='C') {
	  cont_type[i]=TIME_CONTINUE;
	  bSort=FALSE;
	  ok=TRUE;
	  settime[i]=FLT_MAX;
	}
	else if(inputstring[0]=='l' ||
		inputstring[0]=='L') {
	  cont_type[i]=TIME_LAST;
	  bSort=FALSE;
	  ok=TRUE;
	  settime[i]=FLT_MAX;			  
	}
	else {
	  settime[i]=strtod(inputstring,&chptr);
	  if(chptr==inputstring) {
	    fprintf(stderr,"Try that again: ");
	  }
	  else {
	    cont_type[i]=TIME_EXPLICIT;
	    ok=TRUE;
	  }
	}
      } while (!ok);
    }
    if(cont_type[0]!=TIME_EXPLICIT) {
      cont_type[0]=TIME_EXPLICIT;
      settime[0]=0;
    }
  }
  else 
    for(i=0;i<nfiles;i++)
      settime[i]=readtime[i];
  
  if(bSort && (nfiles>1)) 
    sort_files(fnms,settime,nfiles);
  else
    fprintf(stderr,"Sorting disabled.\n");
  
  
  /* Write out the new order and start times */
  fprintf(stderr,"\nSummary of files and start times used:\n\n"
	  "          File                Start time\n"
	  "-----------------------------------------\n");
  for(i=0;i<nfiles;i++)
    switch(cont_type[i]) {
    case TIME_EXPLICIT:
      fprintf(stderr,"%25s   %10.3f\n",fnms[i],settime[i]);
      break;
    case TIME_CONTINUE:
      fprintf(stderr,"%25s        Continue from end of last file\n",fnms[i]);
      break;	      
    case TIME_LAST:
      fprintf(stderr,"%25s        Change by same amount as last file\n",fnms[i]);
      break;
    }
  fprintf(stderr,"\n");
  
  settime[nfiles]=FLT_MAX;
  cont_type[nfiles]=TIME_EXPLICIT;
  readtime[nfiles]=FLT_MAX;
}


static void copy_ee(t_energy *src, t_energy *dst, int nre)
{
  int i;

  for(i=0;i<nre;i++) {
    dst[i].e=src[i].e;
    dst[i].esum=src[i].esum;  
    dst[i].eav=src[i].eav;
  }
}


static void remove_last_eeframe(t_energy *lastee, gmx_step_t laststep,
				t_energy *ee, int nre)
{
    int i;
    gmx_step_t p=laststep+1;
    double sigmacorr;
    
    for(i=0;i<nre;i++) {
	lastee[i].esum-=ee[i].e;
	sigmacorr=lastee[i].esum-(p-1)*ee[i].e;
	lastee[i].eav-=(sigmacorr*sigmacorr)/((p-1)*p);
    }
}



static void update_ee(t_energy *lastee,gmx_step_t laststep,
		      t_energy *startee,gmx_step_t startstep,
		      t_energy *ee, int step,
		      t_energy *outee, int nre)
{
  int i; 
  double sigmacorr,nom,denom;
  double prestart_esum;
  double prestart_sigma;
  
  for(i=0;i<nre;i++) {
	  outee[i].e=ee[i].e;
      /* add statistics from earlier file if present */
      if(laststep>0) {
	  outee[i].esum=lastee[i].esum+ee[i].esum;
	  nom=(lastee[i].esum*(step+1)-ee[i].esum*(laststep));
	  denom=((step+1.0)*(laststep)*(step+1.0+laststep));	
	  sigmacorr=nom*nom/denom;
	  outee[i].eav=lastee[i].eav+ee[i].eav+sigmacorr;
      }  
      else {
	  /* otherwise just copy to output */
	  outee[i].esum=ee[i].esum;
	  outee[i].eav=ee[i].eav; 
      }
      
      /* if we didnt start to write at the first frame
       * we must compensate the statistics for this
       * there are laststep frames in the earlier file
       * and step+1 frames in this one.
       */
    if(startstep>0) {
      gmx_step_t q=laststep+step; 
      gmx_step_t p=startstep+1;
      prestart_esum=startee[i].esum-startee[i].e;
      sigmacorr=prestart_esum-(p-1)*startee[i].e;
      prestart_sigma=startee[i].eav-
	sigmacorr*sigmacorr/(p*(p-1));
      sigmacorr=prestart_esum/(p-1)-
	outee[i].esum/(q);
      outee[i].esum-=prestart_esum;
      if (q-p+1 > 0)
	outee[i].eav=outee[i].eav-prestart_sigma-
	  sigmacorr*sigmacorr*((p-1)*q)/(q-p+1);
    }
 
    if((outee[i].eav/(laststep+step+1))<(GMX_REAL_EPS))
      outee[i].eav=0;
  }
}


static void update_last_ee(t_energy *lastee, gmx_step_t laststep,
			   t_energy *ee,gmx_step_t step,int nre)
{
    t_energy *tmp;
    snew(tmp,nre);
    update_ee(lastee,laststep,NULL,0,ee,step,tmp,nre);
    copy_ee(tmp,lastee,nre);
    sfree(tmp);
}

static void update_ee_sum(int nre,
			  gmx_step_t *ee_sum_step,gmx_step_t *ee_sum_nsum,
			  t_energy *ee_sum,
			  t_enxframe *fr,int out_step)
{
  gmx_step_t ns,fr_nsum;
  int i;
 
  ns = *ee_sum_nsum;

  fr_nsum = fr->nsum;
  if (fr_nsum == 0) {
    fr_nsum = 1;
  }

  if (ns == 0) {
    if (fr_nsum == 1) {
      for(i=0;i<nre;i++) {
	ee_sum[i].esum = fr->ener[i].e;
	ee_sum[i].eav  = 0;
      }
    } else {
      for(i=0;i<nre;i++) {
	ee_sum[i].esum = fr->ener[i].esum;
	ee_sum[i].eav  = fr->ener[i].eav;
      }
    }
    ns = fr_nsum;
  } else if (out_step - *ee_sum_step == ns + fr_nsum) {
    if (fr_nsum == 1) {
      for(i=0;i<nre;i++) {
	ee_sum[i].eav  +=
	  dsqr(ee_sum[i].esum/ns - (ee_sum[i].esum + fr->ener[i].e)/(ns + 1))*
	  ns*(ns + 1);
	ee_sum[i].esum += fr->ener[i].e;
      }
    } else {
      for(i=0; i<fr->nre; i++) {
	ee_sum[i].eav  +=
	  fr->ener[i].eav +
	  dsqr(ee_sum[i].esum/ns - (ee_sum[i].esum + fr->ener[i].esum)/(ns + fr->nsum))*
	  ns*(ns + fr->nsum)/(double)fr->nsum;
	ee_sum[i].esum += fr->ener[i].esum;
      }
    }
    ns += fr_nsum;
  } else {
    if (fr->nsum != 0) {
      fprintf(stderr,"\nWARNING: missing energy sums at time %f\n",fr->t);
    }
    ns = 0;
  }

  *ee_sum_step = out_step;
  *ee_sum_nsum = ns;
}
 
int gmx_eneconv(int argc,char *argv[])
{
  static char *desc[] = {
    "With [IT]multiple files[it] specified for the [TT]-f[tt] option:[BR]",
    "Concatenates several energy files in sorted order.",
    "In case of double time frames the one",
    "in the later file is used. By specifying [TT]-settime[tt] you will be",
    "asked for the start time of each file. The input files are taken",
    "from the command line,",
    "such that the command [TT]eneconv -o fixed.edr *.edr[tt] should do",
    "the trick. [PAR]",
    "With [IT]one file[it] specified for [TT]-f[tt]:[BR]",
    "Reads one energy file and writes another, applying the [TT]-dt[tt],",
    "[TT]-offset[tt], [TT]-t0[tt] and [TT]-settime[tt] options and",
    "converting to a different format if necessary (indicated by file",
    "extentions).[PAR]",
    "[TT]-settime[tt] is applied first, then [TT]-dt[tt]/[TT]-offset[tt]",
    "followed by [TT]-b[tt] and [TT]-e[tt] to select which frames to write."
  };
  static char *bugs[] = {
    "When combining trajectories the sigma and E^2 (necessary for statistics) are not updated correctly. Only the actual energy is correct. One thus has to compute statistics in another way."
  };
  int        in,out=0;
  gmx_enxnm_t *enm=NULL;
  t_enxframe *fr,*fro;
  gmx_step_t ee_sum_step=0,ee_sum_nsum;
  t_energy   *ee_sum;
  gmx_step_t laststep,startstep,startstep_file=0;
  int        noutfr;
  int        nre,nremax,this_nre,nfile,f,i,j,kkk,nset,*set=NULL;
  real       t; 
  char       **fnms;
  real       *readtime,*settime,timestep,t1,tadjust;
  char       inputstring[STRLEN],*chptr,buf[22],buf2[22];
  bool       ok;
  int        *cont_type;
  bool       bNewFile,bFirst,bNewOutput;
  
  t_filenm fnm[] = {
    { efEDR, "-f", NULL,    ffRDMULT },
    { efEDR, "-o", "fixed", ffWRITE  },
  };

#define NFILE asize(fnm)  
  bool   bWrite;
  static real  delta_t=0.0, toffset=0,scalefac=1;
  static bool  bSetTime=FALSE;
  static bool  bSort=TRUE,bError=TRUE;
  static real  begin=-1;
  static real  end=-1;
  
  t_pargs pa[] = {
    { "-b",        FALSE, etREAL, {&begin},
      "First time to use"},
    { "-e",        FALSE, etREAL, {&end},
      "Last time to use"},
    { "-dt",       FALSE, etREAL, {&delta_t},
      "Only write out frame when t MOD dt = offset" },
    { "-offset",   FALSE, etREAL, {&toffset},
      "Time offset for -dt option" }, 
    { "-settime",  FALSE, etBOOL, {&bSetTime}, 
      "Change starting time interactively" },
    { "-sort",     FALSE, etBOOL, {&bSort},
      "Sort energy files (not frames)"},
    { "-scalefac", FALSE, etREAL, {&scalefac},
      "Multiply energy component by this factor" },
    { "-error",    FALSE, etBOOL, {&bError},
      "Stop on errors in the file" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_BE_NICE ,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  tadjust  = 0;
  nremax   = 0;
  nset     = 0;
  timestep = 0.0;
  snew(fnms,argc);
  nfile=0;
  laststep=startstep=0;
  
  nfile = opt2fns(&fnms,"-f",NFILE,fnm);
  
  if (!nfile)
    gmx_fatal(FARGS,"No input files!");
  
  snew(settime,nfile+1);
  snew(readtime,nfile+1);
  snew(cont_type,nfile+1);
  
  nre=scan_ene_files(fnms,nfile,readtime,&timestep,&nremax);   
  edit_files(fnms,nfile,readtime,settime,cont_type,bSetTime,bSort);     

  ee_sum_nsum = 0;
  snew(ee_sum,nremax);

  snew(fr,1);
  snew(fro,1);
  fro->t = -1;
  fro->nre = nre;
  snew(fro->ener,nremax);

  noutfr=0;
  bFirst=TRUE;

  t = -1e20;
  for(f=0;f<nfile;f++) {
    bNewFile=TRUE;
    bNewOutput=TRUE;
    in=open_enx(fnms[f],"r");
    do_enxnms(in,&this_nre,&enm);
    if (f == 0) {
      if (scalefac != 1)
	set = select_it(nre,enm,&nset);
      
      /* write names to the output file */
      out=open_enx(opt2fn("-o",NFILE,fnm),"w");  
      do_enxnms(out,&nre,&enm);
    }
    
    /* start reading from the next file */
    while ((t <= (settime[f+1] + GMX_REAL_EPS)) &&
	  do_enx(in,fr)) {
      if(bNewFile) {
	startstep_file = fr->step;
	tadjust = settime[f] - fr->t;	  
	if(cont_type[f+1]==TIME_LAST) {
	  settime[f+1]   = readtime[f+1]-readtime[f]+settime[f];
	  cont_type[f+1] = TIME_EXPLICIT;
	}
	bNewFile = FALSE;
      }
      
      if (debug) {
	fprintf(debug,"fr->step %s, ee_sum_nsum %s\n",
		gmx_step_str(fr->step,buf),gmx_step_str(ee_sum_nsum,buf2));
      }

      if (tadjust + fr->t <= t) {
	/* Skip this frame, since we already have it / past it */
	continue;
      }

      fro->step = laststep + fr->step - startstep_file;
      t = tadjust + fr->t;

      /*bWrite = ((begin<0 || (begin>=0 && (t >= begin-GMX_REAL_EPS))) && 
		(end  <0 || (end  >=0 && (t <= end  +GMX_REAL_EPS))) &&
		(t < settime[i+1]-GMX_REAL_EPS));*/
      bWrite = ((begin<0 || (begin>=0 && (t >= begin-GMX_REAL_EPS))) && 
		(end  <0 || (end  >=0 && (t <= end  +GMX_REAL_EPS))) &&
		(t <= settime[f+1]+0.5*timestep));
      
      if (bError)      
	if ((end > 0) && (t > end+GMX_REAL_EPS)) {
	  f = nfile;
	  break;
	}
      
      if (t >= begin-GMX_REAL_EPS) {
	if (bFirst) {
	  bFirst = FALSE;
	  startstep = fr->step;	
	}
	update_ee_sum(nre,&ee_sum_step,&ee_sum_nsum,ee_sum,fr,fro->step);
      }	  
      
      /* determine if we should write it */
      if (bWrite && (delta_t==0 || bRmod(t,toffset,delta_t))) {
	fro->t = t;
	if(bNewOutput) {
	  bNewOutput=FALSE;
	  fprintf(stderr,"\nContinue writing frames from t=%g, step=%s\n",
		  t,gmx_step_str(fro->step,buf));
	}

	/* Copy the energies */
	for(i=0; i<nre; i++) {
	  fro->ener[i].e = fr->ener[i].e;
	}

	if (ee_sum_nsum <= 1) {
	  fro->nsum = 0;
	} else {
	  fro->nsum = ee_sum_nsum;
	  /* Copy the energy sums */
	  for(i=0; i<nre; i++) {
	    fro->ener[i].esum = ee_sum[i].esum;
	    fro->ener[i].eav  = ee_sum[i].eav;
	  }
	}
	/* We wrote the sums, so reset the sum count */
	ee_sum_nsum = 0;
	
	if (scalefac != 1) {
	  for(kkk=0; kkk<nset; kkk++) {
	    fro->ener[set[kkk]].e    *= scalefac;
	    if (fro->nsum > 0) {
	      fro->ener[set[kkk]].eav  *= scalefac*scalefac;
	      fro->ener[set[kkk]].esum *= scalefac;
	    }
	  }
	}
	/* Copy restraint stuff */
	fro->ndisre       = fr->ndisre;
	fro->disre_rm3tav = fr->disre_rm3tav;
	fro->disre_rt     = fr->disre_rt;
	fro->nblock       = fr->nblock;
	fro->nr           = fr->nr;
	fro->block        = fr->block;
	
	do_enx(out,fro);
	if (noutfr % 1000 == 0)
	  fprintf(stderr,"Writing frame time %g    ",fro->t);
	noutfr++;
      }
    }

    if (nfile > 1) {
      laststep = fro->step;
      printf("laststep=%s step=%s\n",
	     gmx_step_str(laststep,buf),gmx_step_str(fro->step,buf2));
    }
    
    /* set the next time from the last in previous file */
    if (cont_type[f+1]==TIME_CONTINUE) {
	settime[f+1] = fro->t;
	/* in this case we have already written the last frame of
	 * previous file, so update begin to avoid doubling it
	 * with the start of the next file
	 */
	begin = fro->t+0.5*timestep;
	/* cont_type[f+1]==TIME_EXPLICIT; */
    }
    
    if ((fro->t < end) && (f < nfile-1) &&
	(fro->t < settime[f+1]-1.5*timestep)) 
      fprintf(stderr,
	      "\nWARNING: There might be a gap around t=%g\n",t);
    
    /* move energies to lastee */
    close_enx(in);
    free_enxnms(this_nre,enm);

    fprintf(stderr,"\n");
  }
  if (noutfr == 0)
    fprintf(stderr,"No frames written.\n");
  else {
    fprintf(stderr,"Last frame written was at step %s, time %f\n",
	    gmx_step_str(fro->step,buf),fro->t);
    fprintf(stderr,"Wrote %d frames\n",noutfr);
  }

  thanx(stderr);
  return 0;
}
