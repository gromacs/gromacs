/*
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
 * Grunge ROck MAChoS
 */
static char *SRCID_trjconv_c = "";

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

static void scan_trj_files(char **fnms,int nfiles,real *readtime, real *timestep)
{
    /* Check start time of all files */
    int i,status,oldnatoms=0,natoms,step;
    real t,t2;
    rvec *x;
    matrix box;

    
    for(i=0;i<nfiles;i++) {
      natoms=read_first_x(&status,fnms[i],&t,&x,box);
      if(i==0) {
	oldnatoms=natoms;
	read_next_x(status,&t2,natoms,x,  box);
        *timestep=t2-t;
	sfree(x);
      }
      else {
	if(natoms!=oldnatoms) 
	  fatal_error(0,"Different number of atoms in files");
      }
      readtime[i]=t;
      close_trj(status);
    }
    fprintf(stderr,"\n");  
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


bool bRmod(double a,double b)
{
  int iq;
  double tol = 1e-6;
  
  iq = ((1.0+tol)*a)/b;
  
  if (fabs(a-b*iq) <= tol*a)
    return TRUE;
  else
    return FALSE;
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
      
  FILE         *out=NULL;
  int          status,ftp,ftpin,i,frame,natoms=0,step,trjout=0;
  rvec         *x,*v;
  real         xtcpr,t,t0=-1,t1;
  matrix       box;
  t_topology   top;
  char         **fnms;
  bool         bNewFile,bHaveX=FALSE,bHaveV=FALSE;
  char         *in_file,*out_file;
  int          xdr=0,earliersteps,nfile,*cont_type,last_ok_step;
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
  
  t=-1;
  earliersteps=0;    
  out_file=opt2fn("-o",NFILE,fnm);
  ftp=fn2ftp(out_file);

  bVels= ((ftp==efTRR) ||(ftp==efTRJ) || (ftp==efGRO));
  /* Not checking input format, could be dangerous :-) */
  if (!bVels) {
    bHaveX=TRUE;
    bHaveV=FALSE;
    v=NULL;
  }
  
 
  frame=-1;
  /* Lets stitch up some files */
  for(i=0;i<nfile;i++) {
      /* Open next file */
      if(bVels)
	  natoms=read_first_x_or_v(&status,fnms[i],&tstart,&x,&v,box);
      else
	  natoms=read_first_x(&status,fnms[i],&tstart,&x,box);
      if(cont_type[i]==TIME_EXPLICIT)
	  t0=settime[i]-tstart;
      t1=tstart;
   
      bNewFile=TRUE;
      if(i==0) {
	  switch(ftp) {
	  case efXTC:
	      xdr = open_xtc(out_file,"w");
	      break;
	  case efG87:
	      out=ffopen(out_file,"w");
	      break;
	  case efTRR:
	  case efTRJ:
	      out=NULL;
	      trjout = open_tpx(out_file,"w");
	      break;
	  case efGRO:
	  case efG96:
	  case efPDB:
		  out=ffopen(out_file,"w");
	      break;
	  }
      }
      
      do {
	  /* set the new time */
	  t=t0+t1;
	  if((end>0) && (t>(end+GMX_REAL_EPS))) {
	      i=nfile;
	      break;
	  }
	  /* determine if we should write it */
	  if((t<(settime[i+1]-0.5*timestep)) && (t>=begin)) {
	      frame++;
	      last_ok_t=t;
	      if(bNewFile) {
		  fprintf(stderr,"\nContinue writing frames from t=%g, frame=%d      \n",t,frame);
		  bNewFile=FALSE;
	      }
	      
	      switch(ftp) {
	      case efTRJ:
	      case efTRR:
		  fwrite_trn(trjout,frame,t,0,box,
			     natoms,x,v,NULL);
		  break;
	      case efXTC:
		  write_xtc(xdr,natoms,frame,t,box,x,xtcpr);
		  break;
	      default:
		  fatal_error(0,"DHE, ftp=%d\n",ftp);
	      }
	      
	      fprintf(stderr,"\rWriting frame %d, time %f         ",frame,t);
	  }
      } while((t<settime[i+1]) &&
	      ((bVels && read_next_x_or_v(status,&t1,natoms,x,v,box)) ||
	       (read_next_x(status,&t1,natoms,x,  box))));
      
      close_trj(status);
      fprintf(stderr,"\n");
      
      earliersteps+=step;	  
      
      /* set the next time from the last in previous file */
      if(cont_type[i+1]==TIME_CONTINUE) {
	  begin=t+0.5*timestep;
	  settime[i+1]=t;
	  cont_type[i+1]=TIME_EXPLICIT;	  
      }
      else if(cont_type[i+1]==TIME_LAST)
	  begin=t+0.5*timestep;
      
      if(cont_type[i+1]==TIME_EXPLICIT) 
	if((i<(nfile-1)) &&
	   (t<(settime[i+1]-1.5*timestep))) 
	  fprintf(stderr,
		  "WARNING: Frames around t=%f have a different spacing than the rest,\n"
		  "might be a gap or overlap that couldn't be corrected automatically.\n",t);   
  }
  
  if (ftp == efXTC)
      close_xtc(xdr);
  else if (out != NULL)
      fclose(out);
  fprintf(stderr,"\nLast frame written was %d, time %f\n",frame,last_ok_t); 

  thanx(stdout);
  
  return 0;
}
