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
 * Great Red Oystrich Makes All Chemists Sane
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

static void scan_trj_files(char **fnms,int nfiles,real *readtime, real *timestep)
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
      fatal_error(0,"Couldn't read frame from file.");
    if(!fr.bTime)
      fatal_error(0,"Couldn't find a time in the frame.");
    readtime[i]=fr.time;
    
    if(i==0) {
      natoms=fr.natoms;
      t=fr.time;
      
      ok=read_next_frame(status,&fr);
      if(!ok || !fr.bTime) 
	fatal_error(0,"Couldn't read time from second frame.");
      
      *timestep=fr.time-t;
    }
    else {
      if(natoms!=fr.natoms) 
	fatal_error(0,"Different number of atoms in files");
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
	fprintf(stderr,"\n\nEnter the new start time for each file.\n"
		"There are two special options, both disable sorting:\n\n"
		"c (continue) - The start time is taken from the end\n"
		"of the previous file. Use it when your continuation run\n"
		"restarts with t=0.\n\n"
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
		    settime[i]=strtod(inputstring,&chptr);
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
	    fprintf(stderr,"%25s   %10.3f\n",fnms[i],settime[i]);
	    break;
	case TIME_CONTINUE:
	    fprintf(stderr,"%25s        Continue from last file\n",fnms[i]);
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
      "trjcat concatenates several input trajectory files in sorted order. ",
      "In case of double time frames the one in the later file is used. ",
      "By specifying [TT]-settime[tt] you will be asked for the start time ",
      "of each file. The input files are taken from the command line, ",
      "such that a command like [TT]trjconv -o fixed.trr *.trr[tt] should do ",
      "the trick."
  };
  static bool  bVels=TRUE;
  static int   prec=3;
  static bool  bSort=TRUE;
  static bool  bSetTime=FALSE;
  static real  begin=-1;
  static real  end=-1;
  
  t_pargs pa[] = {
      { "-b",        FALSE, etREAL, {&begin},
	"First time to use"},
      { "-e",        FALSE, etREAL, {&end},
	"Last time to use"},
    { "-prec", FALSE,  etINT,  {&prec},
      "Precision for .xtc and .gro writing in number of decimal places" },
    { "-vel", FALSE, etBOOL, {&bVels},
      "Read and write velocities if possible" },
    { "-settime",  FALSE, etBOOL, {&bSetTime}, 
      "Change starting time interactively" },
    { "-sort",     FALSE, etBOOL, {&bSort},
      "Sort trajectory files (not frames)"}
  };
      
  int          status,ftp,ftpin,i,frame,step,trjout=0;
  rvec         *x,*v;
  real         xtcpr,t_corr;
  t_trxframe   fr,frout;
  char         **fnms;
  int          trxout=-1;
  bool         bNewFile;
  char         *in_file,*out_file;
  int          flags,earliersteps,nfile,*cont_type,last_ok_step;
  real         *readtime,*settime,tstart,last_ok_t=-1,timestep;

  t_filenm fnm[] = {
      { efTRX, "-o", "trajout", ffWRITE }
  };
  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

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
  
  snew(settime,nfile+1);
  snew(readtime,nfile+1);
  snew(cont_type,nfile+1);
  
  scan_trj_files(fnms,nfile,readtime,&timestep);
  edit_files(fnms,nfile,readtime,settime,cont_type,bSetTime,bSort);
  
  earliersteps=0;    
  out_file=opt2fn("-o",NFILE,fnm);
  ftp=fn2ftp(out_file);

  /* Not checking input format, could be dangerous :-) */
 
  flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
 
  frame=-1;
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
	  if(frout.time < settime[i+1]-0.5*timestep && frout.time >= begin) {
	    frame++;
	    last_ok_t=frout.time;
	    if(bNewFile) {
	      fprintf(stderr,"\nContinue writing frames from t=%g, frame=%d      \n",frout.time,frame);
	      bNewFile=FALSE;
	    }
	    
	    switch(ftp) {
	    case efTRJ:
	    case efTRR:
	      case efXTC:
		write_trxframe(trxout,&frout);
		break;
	    default:
	      fatal_error(0,"This fileformat doesn't work here yet.");
	      }
	    if ( ((frame % 10) == 0) || (frame < 10) )
	      fprintf(stderr," ->  frame %6d time %8.3f      \r",
		      frame,frout.time);
	  }
      } while((frout.time<settime[i+1]) && read_next_frame(status,&fr));
      
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
	  fprintf(stderr,
		  "WARNING: Frames around t=%f have a different spacing than the rest,\n"
		  "might be a gap or overlap that couldn't be corrected automatically.\n",frout.time);   
  }
  if (trxout >= 0)
    close_trx(trxout);
     
  fprintf(stderr,"\nLast frame written was %d, time %f\n",frame,last_ok_t); 

  thanx(stdout);
  
  return 0;
}
