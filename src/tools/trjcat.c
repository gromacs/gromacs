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
static char *SRCID_trjcat_c = "$Id$";
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "macros.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "copyrite.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "statutil.h"
#include "futil.h"
#include "pdbio.h"
#include "confio.h"
#include "names.h"
#include "rdgroup.h"
#include "vec.h"
#include "xtcio.h"
#include "do_fit.h"
#include "rmpbc.h"
#include "wgms.h"
#include "magic.h"
#include "binio.h"
#include "pbc.h"

#define TIME_EXPLICIT 0
#define TIME_CONTINUE 1
#define TIME_LAST     2
#ifndef FLT_MAX
#define FLT_MAX 1e36
#endif

static real get_timestep(char *fnm)
{
  /* read first two frames in trajectory 'fnm' to determine timestep */
  
  int        status;
  real       t0, dt;
  t_trxframe fr;
  bool ok;
  
  ok=read_first_frame(&status,fnm,&fr,TRX_NEED_X);
  if(!ok || !fr.bTime)
    fatal_error(0,"\nCouldn't read time from first frame.");
  t0=fr.time;
    
  ok=read_next_frame(status,&fr);
  if(!ok || !fr.bTime) 
    fatal_error(0,"\nCouldn't read time from second frame.");
  dt=fr.time-t0;

  close_trj(status);
  
  return dt;
}

static void scan_trj_files(char **fnms,int nfiles,real *readtime,atom_id imax)
{
  /* Check start time of all files */
  int i,flags,status,natoms=0;
  real t;
  t_trxframe fr;
  bool ok;
  
  flags = TRX_NEED_X;
  
  for(i=0;i<nfiles;i++) {
    ok=read_first_frame(&status,fnms[i],&fr,flags);
    
    if(!ok) 
      fatal_error(0,"\nCouldn't read frame from file.");
    if(!fr.bTime)
      fatal_error(0,"\nCouldn't find a time in the frame.");
    readtime[i]=fr.time;
    
    if(i==0) {
      natoms=fr.natoms;
    }
    else {
      if (imax==NO_ATID) {
	if(natoms!=fr.natoms) 
	  fatal_error(0,"\nDifferent numbers of atoms (%d/%d) in files",
		      natoms,fr.natoms);
      } else {
	if(fr.natoms <= imax)
	  fatal_error(0,"\nNot enough atoms (%d) for index group (%d)",
		      fr.natoms,imax);
      }
    }
    close_trj(status);
  }
  fprintf(stderr,"\n");
  
  sfree(fr.x);
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


static void edit_files(char **fnms,int nfiles,real *readtime, real
		       *settime, int *cont_type, bool bSetTime,bool bSort)
{
    int i;
    bool ok;
    char inputstring[STRLEN],*chptr;
    
    if(bSetTime) {
	fprintf(stderr,"\n\nEnter the new start time (%s) for each file.\n"
		"There are two special options, both disable sorting:\n\n"
		"c (continue) - The start time is taken from the end\n"
		"of the previous file. Use it when your continuation run\n"
		"restarts with t=0.\n\n"
		"l (last) - The time in this file will be changed the\n"
		"same amount as in the previous. Use it when the time in the\n"
		"new run continues from the end of the previous one,\n"
		"since this takes possible overlap into account.\n\n",
		time_unit() );
	
	  fprintf(stderr,
	  "          File             Current start (%s)  New start (%s)\n"
		  "---------------------------------------------------------\n",
		  time_unit(), time_unit() );
	  
	  for(i=0;i<nfiles;i++) {
	      fprintf(stderr,"%25s   %10.3f %s          ",
		      fnms[i],convert_time(readtime[i]), time_unit());
	      ok=FALSE;
	      do {
		  fgets(inputstring,STRLEN-1,stdin);
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
		    settime[i]=strtod(inputstring,&chptr)/time_factor();
		    if(chptr==inputstring) {
		      fprintf(stderr,"'%s' not recognized as a floating point number, 'c' or 'l'. "
			      "Try again: ",inputstring);
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
    
    if(!bSort) 
	fprintf(stderr,"Sorting disabled.\n");
    else 
	sort_files(fnms,settime,nfiles);
    
    
    /* Write out the new order and start times */
    fprintf(stderr,"\nSummary of files and start times used:\n\n"
	    "          File                Start time\n"
	    "-----------------------------------------\n");
    for(i=0;i<nfiles;i++)
	switch(cont_type[i]) {
	case TIME_EXPLICIT:
	  fprintf(stderr,"%25s   %10.3f %s",
		  fnms[i],convert_time(settime[i]),time_unit());
	  if ( i>0 && 
	       cont_type[i-1]==TIME_EXPLICIT && settime[i]==settime[i-1] )
	    fprintf(stderr," WARNING: same Start time as previous");
	  fprintf(stderr,"\n");
	  break;
	case TIME_CONTINUE:
	  fprintf(stderr,"%25s        Continue from last file\n",fnms[i]);
	  break;	      
	case TIME_LAST:
	  fprintf(stderr,"%25s        Change by same amount as last file\n",
		  fnms[i]);
	  break;
	}
    fprintf(stderr,"\n");

    settime[nfiles]=FLT_MAX;
    cont_type[nfiles]=TIME_EXPLICIT;
    readtime[nfiles]=FLT_MAX;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
      "trjcat concatenates several input trajectory files in sorted order. ",
      "In case of double time frames the one in the later file is used. ",
      "By specifying [TT]-settime[tt] you will be asked for the start time ",
      "of each file. The input files are taken from the command line, ",
      "such that a command like [TT]trjcat -o fixed.trr *.trr[tt] should do ",
      "the trick. Using [TT]-cat[tt] you can simply paste several files ",
      "together without removal of frames with identical time stamps."
  };
  static bool  bVels=TRUE;
  static int   prec=3;
  static bool  bCat=FALSE;
  static bool  bSort=TRUE;
  static bool  bSetTime=FALSE;
  static real  begin=-1;
  static real  end=-1;
  static real  dt=0;
  
  t_pargs pa[] = {
    { "-b",       FALSE, etTIME, {&begin},
      "First time to use (%t)"},
    { "-e",       FALSE, etTIME, {&end},
      "Last time to use (%t)"},
    { "-dt",      FALSE, etTIME, {&dt},
      "Only write frame when t MOD dt = first time (%t)" },
    { "-prec",    FALSE, etINT,  {&prec},
      "Precision for .xtc and .gro writing in number of decimal places" },
    { "-vel",     FALSE, etBOOL, {&bVels},
      "Read and write velocities if possible" },
    { "-settime", FALSE, etBOOL, {&bSetTime}, 
      "Change starting time interactively" },
    { "-sort",    FALSE, etBOOL, {&bSort},
      "Sort trajectory files (not frames)" },
    { "-cat",     FALSE, etBOOL, {&bCat},
      "do not discard double time frames" }
  };
      
  int         status,ftp,ftpin,i,frame,frame_out,step,trjout=0;
  rvec        *x,*v;
  real        xtcpr,t_corr;
  t_trxframe  fr,frout;
  char        **fnms;
  int         trxout=-1;
  bool        bNewFile,bIndex;
  char        *in_file,*out_file;
  int         flags,earliersteps,nfile,*cont_type,last_ok_step;
  real        *readtime,*settime,first_time=0,last_ok_t=-1,timestep;
  int         isize;
  atom_id     *index=NULL,imax;
  char        *grpname;

  t_filenm fnm[] = {
      { efTRX, "-o", "trajout", ffWRITE },
      { efNDX, NULL,  NULL,     ffOPTRD }
  };
  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS|PCA_BE_NICE|PCA_TIME_UNIT,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

  bIndex=ftp2bSet(efNDX,NFILE,fnm);
  
  imax=NO_ATID;
  if (bIndex) {
    printf("Select group for output\n");
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
    /* scan index */
    imax=index[0];
    for(i=1; i<isize; i++)
      imax = max(imax, index[i]);
  }

  /* prec is in nr of decimal places, xtcprec is a multiplication factor: */
  xtcpr=1;
  for (i=0; i<prec; i++)
      xtcpr*=10;
  
  nfile=0;
      
  snew(fnms,argc);
  for(i=1; (i<argc); i++) 
      fnms[nfile++]=argv[i];
  
  if(!nfile)
      fatal_error(0,"No input files!");
  
  timestep = get_timestep(fnms[0]);
  snew(readtime,nfile+1);
  scan_trj_files(fnms,nfile,readtime,imax);
  
  snew(settime,nfile+1);
  snew(cont_type,nfile+1);
  edit_files(fnms,nfile,readtime,settime,cont_type,bSetTime,bSort);
  
  earliersteps=0;    
  out_file=opt2fn("-o",NFILE,fnm);
  ftp=fn2ftp(out_file);
  
  /* Not checking input format, could be dangerous :-) */
 
  flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
 
  frame=-1;
  frame_out=-1;
  /* the default is not to change the time at all,
   * but this is overridden by the edit_files routine
   */
  t_corr=0;
     
  /* Lets stitch up some files */
  for(i=0;i<nfile;i++) {
      /* Open next file */
  
      read_first_frame(&status,fnms[i],&fr,flags);
      if(!fr.bTime)
	fatal_error(0,"Couldn't find a time in the frame.");
      
      if(cont_type[i]==TIME_EXPLICIT)
	t_corr=settime[i]-fr.time;
      /* t_corr is the amount we want to change the time.
       * If the user has chosen not to change the time for
       * this part of the trajectory t_corr remains at 
       * the value it had in the last part, changing this
       * by the same amount.
       * If no value was given for the first trajectory part
       * we let the time start at zero, see the edit_files routine.
       */
      
      bNewFile=TRUE;
      if(i==0) {
	switch(ftp) {
	case efXTC:
	case efTRR:
	case efTRJ:
	  trxout = open_trx(out_file,"w");
	  break;
	default:
	  fatal_error(0,"This fileformat doesn't work here yet.");
	  
	  break;
	}
      }
      
      do {
	/* copy the input frame to the output frame */
	frout=fr;
	/* set the new time by adding the correct calculated above */
	frout.time += t_corr; 
	/* quit if we have reached the end of what should be written */
	if((end > 0) && (frout.time > end+GMX_REAL_EPS)) {
	  i=nfile;
	  break;
	}
	/* determine if we should write this frame */
	if( ( bCat || 
	      ( !bCat && frout.time < settime[i+1]-0.5*timestep ) ) && 
	    frout.time >= begin) {
	  frame++;
	  if (frame_out == -1)
	    first_time = frout.time;
	  if (dt==0 || bRmod(frout.time-first_time,dt)) {
	    frame_out++;
	    last_ok_t=frout.time;
	    if(bNewFile) {
	      fprintf(stderr,"\nContinue writing frames from %s t=%g %s, "
		      "frame=%d      \n",
		      fnms[i],convert_time(frout.time),time_unit(),frame);
	      bNewFile=FALSE;
	    }
	    
	    switch(ftp) {
	    case efTRJ:
	    case efTRR:
	    case efXTC:
	      if (bIndex)
		write_trxframe_indexed(trxout,&frout,isize,index);
	      else
		write_trxframe(trxout,&frout);
	      break;
	    default:
	      fatal_error(0,"This fileformat doesn't work here yet.");
	    }
	    if ( ((frame % 10) == 0) || (frame < 10) )
	      fprintf(stderr," ->  frame %6d time %8.3f %s     \r",
		      frame_out,convert_time(frout.time),time_unit());
	  }
	}
      } while( ( bCat || ( !bCat && frout.time<settime[i+1]) ) && 
	       read_next_frame(status,&fr));
      
      close_trj(status);
      
      earliersteps+=step;	  
      
      /* set the next time from the last frame in previous file */
      if(cont_type[i+1]==TIME_CONTINUE) {
	begin=frout.time+0.5*timestep;
	settime[i+1]=frout.time;
	cont_type[i+1]=TIME_EXPLICIT;	  
      }
      else if(cont_type[i+1]==TIME_LAST)
	begin=frout.time+0.5*timestep;
      /* Or, if the time in the next part should be changed by the
       * same amount, start at half a timestep from the last time
       * so we dont repeat frames.
       */
      
      /* Or, if time is set explicitly, we check for overlap/gap */
      if(cont_type[i+1]==TIME_EXPLICIT) 
	if((i<(nfile-1)) &&
	   (frout.time<(settime[i+1]-1.5*timestep))) 
	  fprintf(stderr, "WARNING: Frames around t=%f %s have a different "
		  "spacing than the rest,\n"
		  "might be a gap or overlap that couldn't be corrected "
		  "automatically.\n",convert_time(frout.time),time_unit());
  }
  if (trxout >= 0)
    close_trx(trxout);
     
  fprintf(stderr,"\nLast frame written was %d, time %f %s\n",
	  frame,convert_time(last_ok_t),time_unit()); 

  thanx(stderr);
  
  return 0;
}
