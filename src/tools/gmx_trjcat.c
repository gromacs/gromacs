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
#include "index.h"
#include "vec.h"
#include "xtcio.h"
#include "do_fit.h"
#include "rmpbc.h"
#include "wgms.h"
#include "magic.h"
#include "pbc.h"
#include "xvgr.h"

#define TIME_EXPLICIT 0
#define TIME_CONTINUE 1
#define TIME_LAST     2
#ifndef FLT_MAX
#define FLT_MAX 1e36
#endif
#define FLAGS (TRX_READ_X | TRX_READ_V | TRX_READ_F)

static void scan_trj_files(char **fnms,int nfiles,
			   real *readtime,real *timestep,atom_id imax)
{
  /* Check start time of all files */
  int i,status,natoms=0;
  real t;
  t_trxframe fr;
  bool ok;
  
  for(i=0;i<nfiles;i++) {
    ok=read_first_frame(&status,fnms[i],&fr,FLAGS);
    
    if(!ok) 
      gmx_fatal(FARGS,"\nCouldn't read frame from file.");
    if(fr.bTime)
      readtime[i]=fr.time;
    else {
      readtime[i]=0;
      fprintf(stderr,"\nWARNING: Couldn't find a time in the frame.\n");
    }
    
    if(i==0) {
      natoms=fr.natoms;
    }
    else {
      if (imax==NO_ATID) {
	if(natoms!=fr.natoms) 
	  gmx_fatal(FARGS,"\nDifferent numbers of atoms (%d/%d) in files",
		      natoms,fr.natoms);
      } else {
	if(fr.natoms <= imax)
	  gmx_fatal(FARGS,"\nNot enough atoms (%d) for index group (%d)",
		      fr.natoms,imax);
      }
    }
    ok=read_next_frame(status,&fr);
    if(ok && fr.bTime) {
      timestep[i] = fr.time - readtime[i];
    } else {
      timestep[i] = 0;
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


static void edit_files(char **fnms,int nfiles,real *readtime, real *timestep,
		       real *settime, int *cont_type, bool bSetTime,bool bSort)
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
		    settime[i]=strtod(inputstring,&chptr)*time_invfactor();
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
	    "          File                Start time       Time step\n"
	    "---------------------------------------------------------\n");
    for(i=0;i<nfiles;i++)
	switch(cont_type[i]) {
	case TIME_EXPLICIT:
	  fprintf(stderr,"%25s   %10.3f %s   %10.3f %s",
		  fnms[i],
		  convert_time(settime[i]),time_unit(),
		  convert_time(timestep[i]),time_unit());
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

static void do_demux(int nset,char *fnms[],int nfile_out,char *fnms_out[],
		     int nval,real **value,real *time,real dt_remd,
		     int isize,atom_id index[])
{
  int        i,j,k,natoms,nnn;
  int        *fp_in,*fp_out;
  bool       bCont,*bSet;
  real       t;
  t_trxframe *trx;
  
  snew(fp_in,nset);
  snew(trx,nset);
  snew(bSet,nset);
  natoms = -1;
  t = -1;
  for(i=0; (i<nset); i++) {
    nnn = read_first_frame(&(fp_in[i]),fnms[i],&(trx[i]),TRX_NEED_X);
    if (natoms == -1)
      natoms = nnn;
    else if (natoms != nnn) 
      gmx_fatal(FARGS,"Trajectory file %s has %d atoms while previous trajs had %d atoms",fnms[i],nnn,natoms);
    if (t == -1)
      t = trx[i].time;
    else if (t != trx[i].time) 
      gmx_fatal(FARGS,"Trajectory file %s has time %f while previous trajs had time %f",fnms[i],trx[i].time,t);
  }
   
  snew(fp_out,nfile_out);
  for(i=0; (i<nfile_out); i++)
    fp_out[i] = open_trx(fnms_out[i],"w");
  k = 0;
  if (gmx_nint(time[k] - t) != 0) 
    gmx_fatal(FARGS,"First time in demuxing table does not match trajectories");
    
  do {
    while ((k+1 < nval) && ((trx[0].time - time[k+1]) > dt_remd*0.1))
      k++;
    if (debug)
      fprintf(debug,"trx[0].time = %g, time[k] = %g\n",trx[0].time,time[k]);
    for(i=0; (i<nfile_out); i++) 
      bSet[i] = FALSE;
    for(i=0; (i<nfile_out); i++) {
      j = gmx_nint(value[i][k]);
      range_check(j,0,nset);
      if (bSet[j])
	gmx_fatal(FARGS,"Demuxing the same replica %d twice at time %f",
		  j,trx[0].time);
      bSet[j] = TRUE;
      if (index)
	write_trxframe_indexed(fp_out[j],&trx[i],isize,index);
      else
	write_trxframe(fp_out[j],&trx[i]);
    }
    
    bCont = (k < nval);
    for(i=0; (i<nset); i++) 
      bCont = bCont && read_next_frame(fp_in[i],&trx[i]);
  } while (bCont);
  
  for(i=0; (i<nset); i++) 
    close_trx(fp_in[i]);
  for(i=0; (i<nfile_out); i++) 
    close_trx(fp_out[i]);
}

int gmx_trjcat(int argc,char *argv[])
{
  static char *desc[] = {
      "trjcat concatenates several input trajectory files in sorted order. ",
      "In case of double time frames the one in the later file is used. ",
      "By specifying [TT]-settime[tt] you will be asked for the start time ",
      "of each file. The input files are taken from the command line, ",
      "such that a command like [TT]trjcat -o fixed.trr *.trr[tt] should do ",
      "the trick. Using [TT]-cat[tt] you can simply paste several files ",
      "together without removal of frames with identical time stamps.[PAR]",
      "One important option is inferred when the output file is amongst the",
      "input files. In that case that particular file will be appended to",
      "which implies you do not need to store double the amount of data.",
      "Obviously the file to append to has to be the one with lowest starting",
      "time since one can only append at the end of a file.[PAR]",
      "If the [TT]-demux[tt] option is given, the N trajectories that are",
      "read, are written in another order as specified in the xvg file."
      "The xvg file should contain something like:[PAR]",
      "0  0  1  2  3  4  5[BR]",
      "2  1  0  2  3  5  4[BR]",
      "Where the first number is the time, and subsequent numbers point to",
      "trajectory indices.",
      "The frames corresponding to the numbers present at the first line",
      "are collected into the output trajectory. If the number of frames in",
      "the trajectory does not match that in the xvg file then the program",
      "tries to be smart. Beware."
  };
  static bool  bVels=TRUE;
  static int   prec=3;
  static bool  bCat=FALSE;
  static bool  bSort=TRUE;
  static bool  bKeepLast=FALSE;
  static bool  bSetTime=FALSE;
  static bool  bDeMux;
  static int   nrepl=1;
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
    { "-nrepl",   FALSE, etINT,  {&nrepl},
      "Number of replicas assumed to be numbered consecutively. Only used with the -demux option" },
    { "-sort",    FALSE, etBOOL, {&bSort},
      "Sort trajectory files (not frames)" },
    { "-keeplast",FALSE, etBOOL, {&bKeepLast},
      "keep overlapping frames at end of trajectory" },
    { "-cat",     FALSE, etBOOL, {&bCat},
      "do not discard double time frames" }
  };
#define npargs asize(pa)
  int         status,ftpin,i,frame,frame_out,step=0,trjout=0;
  rvec        *x,*v;
  real        xtcpr,t_corr;
  t_trxframe  fr,frout;
  char        **fnms,**fnms_out,*in_file,*out_file;
  int         n_append;
  int         trxout=-1;
  bool        bNewFile,bIndex,bWrite;
  int         earliersteps,nfile_in,nfile_out,*cont_type,last_ok_step;
  real        *readtime,*timest,*settime;
  real        first_time=0,lasttime=NOTSET,last_ok_t=-1,timestep;
  int         isize;
  atom_id     *index=NULL,imax;
  char        *grpname;
  real        **val=NULL,*t=NULL,dt_remd;
  int         n,nset;
  t_filenm fnm[] = {
      { efTRX, "-f",     NULL,      ffRDMULT },
      { efTRX, "-o",     "trajout", ffWRMULT },
      { efNDX, "-n",     "index",   ffOPTRD  },
      { efXVG, "-demux", "remd",    ffOPTRD  }
  };
  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_BE_NICE|PCA_TIME_UNIT,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

  bIndex = ftp2bSet(efNDX,NFILE,fnm);
  bDeMux = ftp2bSet(efXVG,NFILE,fnm) || (nrepl > 1);
  bSort  = bSort && !bDeMux;
  
  imax=NO_ATID;
  if (bIndex) {
    printf("Select group for output\n");
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
    /* scan index */
    imax=index[0];
    for(i=1; i<isize; i++)
      imax = max(imax, index[i]);
  }
  if (bDeMux) {
    nset    = 0;
    dt_remd = 0;
    val=read_xvg_time(opt2fn("-demux",NFILE,fnm),TRUE,
		      opt2parg_bSet("-b",npargs,pa),begin,
		      opt2parg_bSet("-e",npargs,pa),end,
		      1,&nset,&n,&dt_remd,&t);
    printf("Read %d sets of %d points, dt = %g\n\n",nset,n,dt_remd);
  }
  /* prec is in nr of decimal places, xtcprec is a multiplication factor: */
  xtcpr=1;
  for (i=0; i<prec; i++)
    xtcpr*=10;
  
  nfile_in = opt2fns(&fnms,"-f",NFILE,fnm);
  if (!nfile_in)
    gmx_fatal(FARGS,"No input files!");
  else if ((nfile_in == 1) && (nrepl > 1)) {
    /* We'll expand filenames ourselves */
    
  } else if ((nfile_in > 1) && (nrepl != 1) && (nrepl != nfile_in)) {
    gmx_fatal(FARGS,"You passed %d input file names while nrepl is %d",
	      nfile_in,nrepl);
  }
  if (bDeMux && (nfile_in != nset)) 
    gmx_fatal(FARGS,"You have specified %d files and %d entries in the demux table",nfile_in,nset);
    
  nfile_out = opt2fns(&fnms_out,"-o",NFILE,fnm);
  if (!nfile_out)
    gmx_fatal(FARGS,"No output files!");
  if ((nfile_out > 1) && !bDeMux) 
    gmx_fatal(FARGS,"Don't know what to do with more than 1 output file if  not demultiplexing");
  else if (bDeMux && (nfile_out != nset) && (nfile_out != 1))
    gmx_fatal(FARGS,"Number of output files should be 1 or %d (#input files), not %d",nset,nfile_out);

  if (bDeMux) 
    do_demux(nfile_in,fnms,nfile_out,fnms_out,n,val,t,dt_remd,isize,index);
  else {
    snew(readtime,nfile_in+1);
    snew(timest,nfile_in+1);
    scan_trj_files(fnms,nfile_in,readtime,timest,imax);
    
    snew(settime,nfile_in+1);
    snew(cont_type,nfile_in+1);
    edit_files(fnms,nfile_in,readtime,timest,settime,cont_type,bSetTime,bSort);
  
    /* Check whether the output file is amongst the input files 
     * This has to be done after sorting etc.
     */
    out_file = fnms_out[0];
    n_append = -1;
    for(i=0; ((i<nfile_in) && (n_append==-1)); i++) {
      if (strcmp(fnms[i],out_file) == 0) {
	n_append = i;
      }
    }
    if (n_append == 0)
      fprintf(stderr,"Will append to %s rather than creating a new file\n",
	      out_file);
    else if (n_append != -1)
      gmx_fatal(FARGS,"Can only append to the first file which is %s (not %s)",
		fnms[0],out_file);
    
    earliersteps=0;    
    
    /* Not checking input format, could be dangerous :-) */
    /* Not checking output format, equally dangerous :-) */
    
    frame=-1;
    frame_out=-1;
    /* the default is not to change the time at all,
     * but this is overridden by the edit_files routine
     */
    t_corr=0;
    
    if (n_append == -1) {
      trxout = open_trx(out_file,"w");
      memset(&frout,0,sizeof(frout));
    }
    else {
      /* Read file to find what is the last frame in it */
      if (!read_first_frame(&status,out_file,&fr,FLAGS))
	gmx_fatal(FARGS,"Reading first frame from %s",out_file);
      while (read_next_frame(status,&fr))
	;
      close_trj(status);
      lasttime = fr.time;
      bKeepLast = TRUE;
      trxout = open_trx(out_file,"a");
      frout = fr;
    }
    
    /* Lets stitch up some files */
    timestep = timest[0];
    for(i=n_append+1; (i<nfile_in); i++) {
      /* Open next file */
      /* set the next time from the last frame in previous file */
      if (i > 0) {
	if(cont_type[i]==TIME_CONTINUE) {
	  begin =frout.time;
	  begin += 0.5*timestep;
	  settime[i]=frout.time;
	  cont_type[i]=TIME_EXPLICIT;	  
	}
	else if(cont_type[i]==TIME_LAST)
	  begin=frout.time;
	begin += 0.5*timestep;
	/* Or, if the time in the next part should be changed by the
	 * same amount, start at half a timestep from the last time
	 * so we dont repeat frames.
	 */
      
	/* Or, if time is set explicitly, we check for overlap/gap */
	if(cont_type[i]==TIME_EXPLICIT) 
	  if( ( i < nfile_in ) &&
	      ( frout.time < settime[i]-1.5*timestep ) ) 
	    fprintf(stderr, "WARNING: Frames around t=%f %s have a different "
		    "spacing than the rest,\n"
		    "might be a gap or overlap that couldn't be corrected "
		    "automatically.\n",convert_time(frout.time),time_unit());
      }
      
      /* if we don't have a timestep in the current file, use the old one */
      if ( timest[i] != 0 )
	timestep = timest[i];
      
      read_first_frame(&status,fnms[i],&fr,FLAGS);
      if(!fr.bTime) {
	fr.time=0;
	fprintf(stderr,"\nWARNING: Couldn't find a time in the frame.\n");
      }
      
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
      
      printf("\n");
      if (lasttime != NOTSET)
	printf("lasttime %g\n", lasttime);
      
      do {
	/* copy the input frame to the output frame */
	frout=fr;
	/* set the new time by adding the correct calculated above */
	frout.time += t_corr; 
	/* quit if we have reached the end of what should be written */
	if((end > 0) && (frout.time > end+GMX_REAL_EPS)) {
	  i=nfile_in;
	  break;
	}
	
	/* determine if we should write this frame (dt is handled elsewhere) */
	if (bCat) /* write all frames of all files */ 
	  bWrite = TRUE;
	else if ( bKeepLast ) /* write till last frame of this traj
				 and skip first frame(s) of next traj */
	  bWrite = ( frout.time > lasttime+0.5*timestep );
	else /* write till first frame of next traj */
	  bWrite = ( frout.time < settime[i+1]-0.5*timestep );
	
	if( bWrite && (frout.time >= begin) ) {
	  frame++;
	  if (frame_out == -1)
	    first_time = frout.time;
	  lasttime = frout.time;
	  if (dt==0 || bRmod(frout.time,first_time,dt)) {
	    frame_out++;
	    last_ok_t=frout.time;
	    if(bNewFile) {
	      fprintf(stderr,"\nContinue writing frames from %s t=%g %s, "
		      "frame=%d      \n",
		      fnms[i],convert_time(frout.time),time_unit(),frame);
	      bNewFile=FALSE;
	    }
	    
	    if (bIndex)
	      write_trxframe_indexed(trxout,&frout,isize,index);
	    else
	      write_trxframe(trxout,&frout);
	    if ( ((frame % 10) == 0) || (frame < 10) )
	      fprintf(stderr," ->  frame %6d time %8.3f %s     \r",
		      frame_out,convert_time(frout.time),time_unit());
	  }
	}
      } while( read_next_frame(status,&fr));
      
      close_trj(status);
      
      earliersteps+=step;	  
    }
    if (trxout >= 0)
      close_trx(trxout);
    
    fprintf(stderr,"\nLast frame written was %d, time %f %s\n",
	    frame,convert_time(last_ok_t),time_unit()); 
  }
  thanx(stderr);
  
  return 0;
}
