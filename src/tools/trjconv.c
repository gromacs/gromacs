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
 * Grunge ROck MAChoS
 */
static char *SRCID_trjconv_c = "$Id$";

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
#define TIME_DITO     2

static void center_x(rvec x[],matrix box,int n,atom_id index[])
{
  int i,m,ai;
  rvec cmin,cmax,dx;
  
  if (n==0) return;
  copy_rvec(x[index[0]],cmin);
  copy_rvec(x[index[0]],cmax);
  for(i=0; (i<n); i++) {
    ai=index[i];
    for(m=0; (m<DIM); m++) {
      if (x[ai][m] < cmin[m])
	cmin[m]=x[ai][m];
      else if (x[ai][m] > cmax[m])
	cmax[m]=x[ai][m];
    }
  }
  for(m=0; (m<DIM); m++) {
    dx[m]=-(box[m][m]-(cmin[m]+cmax[m]))*0.5;
  }
  for(i=0; (i<n); i++) {
    ai=index[i];
    rvec_dec(x[ai],dx);
  }
}


static void scan_trj_files(char **fnms,int nfiles,real *readtime, real *timestep)
{
    /* Check number of energy terms and start time of all files */
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
		"There are two special options, both disables sorting:\n\n"
		"c (continue) - The start time is taken from the end\n"
		"of the previous file. Use it when your continuation run\n"
		"restarts with t=0.\n\n"
		"d (dito) - The time in this file will be changed the\n"
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
		  else if(inputstring[0]=='d' ||
			  inputstring[0]=='D') {
		    cont_type[i]=TIME_DITO;
		    bSort=FALSE;
		    ok=TRUE;
		    settime[i]=FLT_MAX;			  
		  }
		  else {
		    settime[i]=strtod(inputstring,&chptr);
		    if(chptr==inputstring) {
		      fprintf(stderr,"'%s' not recognized as a floating point number, 'c' or 'd'. "
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
	case TIME_DITO:
	    fprintf(stderr,"%25s        Change by same amount as last file\n",fnms[i]);
	    break;
	}
    fprintf(stderr,"\n");

    settime[nfiles]=FLT_MAX;
    cont_type[nfiles]=TIME_EXPLICIT;
    readtime[nfiles]=FLT_MAX;
}



void check_trn(char *fn)
{
  if ((fn2ftp(fn) != efTRJ)  && (fn2ftp(fn) != efTRR))
    fatal_error(0,"%s is not a trj file, exiting\n",fn);
}

#ifndef _win_
void do_trunc(char *fn, real t0)
{
  int          in;
  FILE         *fp;
  bool         bStop,bOK;
  t_trnheader  sh;
  long         fpos;
  char         yesno[256];
  int          j;
  real         t=0;
  
  if (t0 == -1)
    fatal_error(0,"You forgot to set the truncation time");
  
  /* Check whether this is a .trj file */
  check_trn(fn);
  
  in   = open_trn(fn,"r");
  fp   = fio_getfp(in);
  if (fp == NULL) {
    fprintf(stderr,"Sorry, can not trunc %s, truncation of this filetype is not supported\n",fn);
    close_trn(in);
  } else {
    j    = 0;
    fpos = fio_ftell(in);
    bStop= FALSE;
    while (!bStop && fread_trnheader(in,&sh,&bOK)) {
      fread_htrn(in,&sh,NULL,NULL,NULL,NULL);
      fpos=ftell(fp);
      t=sh.t;
      if (t>=t0) {
	fseek(fp,fpos,SEEK_SET);
	bStop=TRUE;
      }
    }
    if (bStop) {
      fprintf(stderr,"Do you REALLY want to truncate this trajectory (%s) at:\n"
	      "frame %d, time %g, bytes %ld ??? (type YES if so)\n",
	      fn,j,t,fpos);
      scanf("%s",yesno);
      if (strcmp(yesno,"YES") == 0) {
	fprintf(stderr,"Once again, I'm gonna DO this...\n");
	close_trn(in);
	truncate(fn,fpos);
      }
      else {
	fprintf(stderr,"Ok, I'll forget about it\n");
      }
    }
    else {
      fprintf(stderr,"Already at end of file (t=%g)...\n",t);
      close_trn(in);
    }
  }
}
#endif

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
    "trjconv can convert trajectory files in many ways:[BR]",
    "[BB] 1.[bb] from one format to another[BR]",
    "[BB] 2.[bb] select a subset of atoms[BR]",
    "[BB] 3.[bb] remove periodicity from molecules[BR]",
    "[BB] 4.[bb] keep multimeric molecules together[BR]",
    "[BB] 5.[bb] center atoms in the box[BR]",
    "[BB] 6.[bb] fit atoms to reference structure[BR]",
    "[BB] 7.[bb] remove duplicate frames[BR]",
    "[BB] 8.[bb] reduce the number of frames[BR]",
    "[BB] 9.[bb] change the timestamps of the frames (e.g. t0 and delta-t)",
    "[BB]10.[bb] concatenate several files into one[BR]",
    "[PAR]",
    "Currently seven formats are supported for input and output:",
    "[TT].xtc[tt], [TT].trr[tt], [TT].trj[tt], [TT].gro[tt], [TT].g96[tt],",
    "[TT].pdb[tt] and [TT].g87[tt].",
    "The file formats are detected from the file extension.",
    "For [TT].gro[tt] and [TT].xtc[tt] files the output precision ",
    "can be given as a number of ",
    "decimal places. Note that velocities are only supported in ",
    "[TT].trr[tt], [TT].trj[tt], [TT].gro[tt] and [TT].g96[tt] files.[PAR]",
    "The option [TT]-app[tt] can be used to",
    "append output to an existing trajectory file.",
    "No checks are performed to ensure integrity",
    "of the resulting combined trajectory file.",
    "[TT].pdb[tt] files with all frames concatenated can be viewed with",
    "[TT]rasmol -nmrpdb[tt].[PAR]",
    "It is possible to select part of your trajectory and write it out",
    "to a new trajectory file in order to save disk space, e.g. for leaving",
    "out the water from a trajectory of a protein in water.",
    "[BB]ALWAYS[bb] put the original trajectory on tape!",
    "We recommend to use the portable [TT].xtc[tt] format for your analysis",
    "to save disk space and to have portable files.[PAR]",
    "There are two options for fitting the trajectory to a reference",
    "either for essential dynamics analysis or for whatever.",
    "The first option is just plain fitting to a reference structure",
    "in the structure file, the second option is a progressive fit",
    "in which the first timeframe is fitted to the reference structure ",
    "in the structure file to obtain and each subsequent timeframe is ",
    "fitted to the previously fitted structure. This way a continuous",
    "trajectory is generated, which might not be the case when using the",
    "regular fit method, e.g. when your protein undergoes large",
    "conformational transitions.[PAR]",
    "The option [TT]-pbc[tt] sets the type of periodic boundary condition",
    "treatment. [TT]whole[tt] makes broken molecules whole (a run input",
    "file is required). [TT]-pbc[tt] is changed form [TT]none[tt] to",
    "[TT]whole[tt] when [TT]-fit[tt] or [TT]-pfit[tt] is set.",
    "[TT]inbox[tt] puts all the atoms in the box.",
    "[TT]nojump[tt] checks if atoms jump across the box and then puts",
    "them back. This has the effect that all molecules",
    "will remain whole (provided they were whole in the initial",
    "conformation), note that this ensures a continuous trajectory but",
    "molecules may diffuse out of the box. The starting configuration",
    "for this procedure is taken from the structure file, if one is",
    "supplied, otherwise it is the first frame.",
    "Use [TT]-center[tt] to put the system in the center of the box.",
    "This is especially useful for multimeric proteins, since this",
    "procedure will ensure the subunits stay together in the trajectory",
    "(due to PBC, they might be separated), providing they were together",
    "in the initial conformation.[PAR]",
    "With the option [TT]-dt[tt] it is possible to reduce the number of ",
    "frames in the output. This option relies on the accuracy of the times",
    "in your input trajectory, so if these are inaccurate use the",
    "[TT]-timestep[tt]",
    "option to modify the time (this can be done simultaneously).[PAR]",
    "Using [TT]-trunc[tt] trjconv can truncate [TT].trj[tt] in place, i.e.",
    "without copying the file. This is useful when a run has crashed",
    "during disk I/O (one more disk full), or when two contiguous",
    "trajectories must be concatenated without have double frames.[PAR]",
    "Also the option [TT]-checkdouble[tt] may be used to remove all",
    "duplicate frames from such a concatenated trajectory, this is done",
    "by ignoring all frames with a time smaller than or equal to the previous",
    "frame.[PAR]",
    "The option [TT]-dump[tt] can be used to extract a frame at or near",
    "one specific time from your trajectory. [PAR]"
    "When the option [TT]-f[tt] is not given, trjconv concatenates several",
    "input files in sorted order. In case of double time frames the one",
    "in the later file is used. By specifying [TT]-settime[tt] you will be",
    "asked for the start time of each file. The input files are taken",
    "from the command line, such that a command like",
    "[TT]trjconv -o fixed.trr *.trr[tt] should do the trick. The only other",
    "options with any influence in this case is -sort.[PAR]"
  };
  
  static char *pbc_opt[] = { NULL, "none", "whole", "inbox", "nojump", NULL };

  static bool  bAppend=FALSE,bSeparate=FALSE,bVels=TRUE;
  static bool  bCenter=FALSE,bFit=FALSE,bPFit=FALSE,bBox=TRUE;
  static bool  bCheckDouble=FALSE;
  static int   skip_nr=1,prec=3;
  static real  tzero=0.0,delta_t=0.0,timestep=0.0,ttrunc=-1,tdump=-1;
  static rvec  newbox = {0,0,0}, shift = {0,0,0};
  static char  *exec_command=NULL;
  static bool  bSort=TRUE;
  static bool  bSetTime=FALSE;

  t_pargs pa[] = {
    { "-pbc", FALSE,  etENUM, {pbc_opt},
      "PBC treatment" },
    { "-center", FALSE,  etBOOL, {&bCenter},
      "Center atoms in box" },
    { "-box", FALSE, etRVEC, {&newbox},
      "Size for new cubic box (default: read from input)" },
    { "-shift", FALSE, etRVEC, {&shift},
      "All coordinates will be shifted by framenr*shift" },
    { "-fit", FALSE,  etBOOL, {&bFit},
      "Fit molecule to ref structure in the structure file" },
    { "-pfit", FALSE,  etBOOL, {&bPFit},
      "Progressive fit, to the previous fitted structure" },
    { "-prec", FALSE,  etINT,  {&prec},
      "Precision for .xtc and .gro writing in number of decimal places" },
    { "-vel", FALSE, etBOOL, {&bVels},
      "Read and write velocities if possible" },
    { "-skip", FALSE,  etINT, {&skip_nr},
      "Only write every nr-th frame" },
    { "-dt", FALSE,  etREAL, {&delta_t},
      "Only write frame when t MOD dt = first time" },
    { "-t0", FALSE,  etREAL, {&tzero},
      "Starting time for trajectory"
      "(default: don't change)"},
#ifndef _win_
    { "-trunc", FALSE, etREAL, {&ttrunc},
      "Truncate input trj file after this amount of ps" },
#endif
    { "-dump", FALSE, etREAL, {&tdump},
      "Dump frame nearest specified time" },
    { "-g87box", FALSE,  etBOOL, {&bBox},
      "Write a box for .g87" },
    { "-exec", FALSE,  etSTR, {&exec_command},
      "Execute command for every output frame with the frame number "
      "as argument" },
    { "-timestep", FALSE,  etREAL, {&timestep},
      "Change time step between frames" },
    { "-app", FALSE,  etBOOL, {&bAppend},
      "Append output"},
    { "-sep", FALSE,  etBOOL, {&bSeparate},
      "Write each frame to a separate .gro or .pdb file"},
    { "-checkdouble", FALSE, etBOOL, {&bCheckDouble},
      "Only write frames with time larger than previous frame" },
    { "-settime",  FALSE, etBOOL, {&bSetTime}, 
      "change starting time interactively" },
    { "-sort",     FALSE, etBOOL, {&bSort},
      "sort energy files (not frames)"}
  };
      
  FILE         *out=NULL;
  int          trjout=0;
  int          status,ftp,ftpin,file_nr=0;
  rvec         *x,*xn,*xout,*v,*vn,*vout;
  rvec         *xp,x_shift;
  real         xtcpr, lambda,*w_rls=NULL;
  matrix       box;
  int          m,i,d,frame,outframe,natoms=0,nout,nre,step;
#define SKIP 10
  t_topology   top;
  t_atoms      *atoms=NULL,useatoms;
  atom_id      *index;
  char         *grpname,**fnms;
  int          ifit,irms;
  atom_id      *ind_fit,*ind_rms;
  char         *gn_fit,*gn_rms;
  real         t,t1,pt,tshift=0,t0=-1,dt=0.001;
  bool         bPBC,bInBox,bNoJump,bNewFile;
  bool         bCopy,bDoIt,bIndex,bTDump,bChangeTime=FALSE,bTPS=FALSE,bDTset=FALSE;
  bool         bExec,bTimeStep=FALSE,bDumpFrame=FALSE,bToldYouOnce=FALSE;
  bool         bHaveNextFrame,bHaveX=FALSE,bHaveV=FALSE,bSetBox;
  char         *grpnm;
  char         *top_file,*in_file,*out_file,out_file2[256],*charpt;
  char         top_title[256],title[256],command[256],filemode[5];
  int          xdr=0;
  int          earliersteps,nfile,*cont_type,last_ok_step;
  real         *readtime,*settime,tstart,begin=-1,last_ok_t;

  t_filenm fnm[] = {
    { efTRX, "-f",  NULL, ffREAD },
    { efTRX, "-o", "trajout", ffWRITE },
    { efTPS, NULL,  NULL, ffOPTRD },
    { efNDX, NULL,  NULL, ffOPTRD }
  };  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME |PCA_NOEXIT_ON_ARGS,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

        /* prec is in nr of decimal places, xtcprec is a multiplication factor: */
      xtcpr=1;
      for (i=0; i<prec; i++)
	xtcpr*=10;

  
  if(!opt2bSet("-f",NFILE,fnm)) {
      fprintf(stderr,"Concatenation/fix mode.\n");
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
    ftpin=fn2ftp(in_file);
    bVels= ((ftp==efTRR) ||(ftp==efTRJ) || (ftp==efGRO)) 
      && ((ftpin==efTRR) ||(ftpin==efTRJ) || (ftpin==efGRO));
    if (!bVels) {
      bHaveX=TRUE;
      bHaveV=FALSE;
    }

    if(!bVels)
	v=NULL;
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
	      if (!bSeparate) 
		  out=ffopen(out_file,"w");
	      break;
	  }
      }

      do {
	  /* set the new time */
	  t=t0+t1;
	  
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
	  cont_type[i+1]==TIME_EXPLICIT;
      }
      	else if(cont_type[i+1]==TIME_DITO)
	    begin=t+0.5*timestep;

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
    
  }
  else {
    top_file=ftp2fn(efTPS,NFILE,fnm);
    
    /* Check command line */
    in_file=opt2fn("-f",NFILE,fnm);
    if (ttrunc != -1) {
#ifndef _win_
      do_trunc(in_file,ttrunc);
#endif
    }
    else {
      if (bPFit)
	bFit = TRUE;
      bSetBox   = opt2parg_bSet("-box", asize(pa), pa);
      bChangeTime= opt2parg_bSet("-t0", asize(pa), pa);
      bExec     = opt2parg_bSet("-exec", asize(pa), pa);
      bTimeStep = opt2parg_bSet("-timestep", asize(pa), pa);
      bTDump    = opt2parg_bSet("-dump", asize(pa), pa);
      bPBC      = (strcmp(pbc_opt[0],"whole") == 0);
      bInBox    = (strcmp(pbc_opt[0],"inbox") == 0);
      bNoJump   = (strcmp(pbc_opt[0],"nojump") == 0);
      if (bFit && (strcmp(pbc_opt[0],"none") == 0))
	  bPBC = TRUE;
      
      bIndex=ftp2bSet(efNDX,NFILE,fnm);
      
      /* Determine output type */ 
      out_file=opt2fn("-o",NFILE,fnm);
      ftp=fn2ftp(out_file);
      fprintf(stderr,"Will write %s: %s\n",ftp2ext(ftp),ftp2desc(ftp));
      if (bVels) {
	/* check if velocities are possible in input and output files */
	ftpin=fn2ftp(in_file);
	bVels= ((ftp==efTRR) ||(ftp==efTRJ) || (ftp==efGRO)) 
	  && ((ftpin==efTRR) ||(ftpin==efTRJ) || (ftpin==efGRO));
      }
      if (!bVels) {
	bHaveX=TRUE;
	bHaveV=FALSE;
      }
      
      /* skipping */  
      if (skip_nr <= 0) {
	fprintf(stderr,"The number of frames to skip is <= 0. "
		"Writing out all frames.\n\n");
	skip_nr = 1;
      } 
      
      /* Determine whether to read a topology */
      bTPS = (ftp2bSet(efTPS,NFILE,fnm) || 
	      bPBC || bFit || (ftp == efGRO) || (ftp == efPDB));
      
      /* Determine if when can read index groups */
      bIndex = (bIndex || bTPS);
      
      if (bTPS) {
	read_tps_conf(top_file,top_title,&top,&xp,NULL,box,bFit);
	atoms=&top.atoms;
	/* top_title is only used for gro and pdb,
	 * the header in such a file is top_title t= ...
	 * to prevent a double t=, remove it from top_title
	 */
	if ((charpt=strstr(top_title," t= ")))
	  charpt[0]='\0';
      }
      
      if (bFit) {
	printf("Select group for root least squares fit\n");
	get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),
		  1,&ifit,&ind_fit,&gn_fit);
	
	if (ifit < 3) 
	  fatal_error(0,"Need at least 3 points to fit!\n");
      }
      
      if (bIndex) {
	printf("Select group for output\n");
	get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),
		  1,&nout,&index,&grpnm);
      }
      else {
	/* no index file, so read natoms from TRX */
	natoms=read_first_x(&status,in_file,&t,&x,box);
	close_trj(status);
	snew(index,natoms);
	for(i=0;i<natoms;i++)
	  index[i]=i;
	nout=natoms; 
      }
      
      if (natoms == 0)
	fatal_error(0,"No atoms found in file %s",in_file);
      
      /* if xp was not snew-ed before, do it now */
      if (!xp)
	snew(xp, natoms);
      
      if (bFit) {
	snew(w_rls,atoms->nr);
	for(i=0; (i<ifit); i++)
	  w_rls[ind_fit[i]]=atoms->atom[ind_fit[i]].m;
	
	/* Restore reference structure and set to origin, 
	   store original location (to put structure back) */
	rm_pbc(&(top.idef),atoms->nr,box,xp,xp);
	copy_rvec(xp[index[0]],x_shift);
	reset_x(ifit,ind_fit,nout,index,xp,w_rls);
	rvec_dec(x_shift,xp[index[0]]);
      }
      
      /* Make atoms struct for output in GRO or PDB files */
      if ((ftp == efGRO) || ((ftp == efG96) && bTPS) || (ftp == efPDB)) {
	/* get memory for stuff to go in pdb file */
	init_t_atoms(&useatoms,atoms->nr,FALSE);
	sfree(useatoms.resname);
	useatoms.resname=atoms->resname;
	for(i=0;(i<nout);i++) {
	  useatoms.atomname[i]=atoms->atomname[index[i]];
	  useatoms.atom[i].resnr=atoms->atom[index[i]].resnr;
	  useatoms.nres=max(useatoms.nres,useatoms.atom[i].resnr+1);
	}
	useatoms.nr=nout;
      }
      /* open trj file for reading */
      if (bVels)
	natoms=read_first_x_or_v(&status,in_file,&t,&x,&v,box);
      else
	natoms=read_first_x     (&status,in_file,&t,&x,   box);
      if (bChangeTime)
	tshift=tzero-t;
      else
	tzero=t;
      
      /* open output for writing */
      if ((bAppend) && (fexist(out_file))) {
	strcpy(filemode,"a");
	fprintf(stderr,"APPENDING to existing file %s\n",out_file);
      } else
	strcpy(filemode,"w");
      switch (ftp) {
      case efXTC:
	xdr = open_xtc(out_file,filemode);
	break;
      case efG87:
	out=ffopen(out_file,filemode);
	break;
      case efTRR:
      case efTRJ:
	out=NULL;
	trjout = open_tpx(out_file,filemode);
	break;
      case efGRO:
      case efG96:
      case efPDB:
	if (!bSeparate) 
	  out=ffopen(out_file,filemode);
	break;
      }
      
      bCopy=FALSE;
      if (bIndex)
	/* check if index is meaningful */
	for(i=0; i<nout; i++) {
	  if (index[i] >= natoms)
	    fatal_error(0,"Index[%d] %d is larger than the number of atoms in the trajectory file (%d)",i,index[i]+1,natoms);
	  bCopy = bCopy || (i != index[i]);
	}
      
      if (bCopy) {
	snew(xn,nout);
	xout=xn;
	snew(vn,nout);
	vout=vn;
      } else {
	xout=x;
	vout=v;
      }
      
      if (ftp == efG87)
	fprintf(out,"Generated by %s. #atoms=%d, %s BOX is stored in "
		"this file.\n",Program(),nout,bBox ? "a" : "NO");
      
      /* Start the big loop over frames */
      file_nr  =  0;  
      frame    =  0;
      outframe =  0;
      pt       = -666.0;
      bDTset   = FALSE;
      
      do {
	/* generate new box */
	if (bSetBox) {
	  clear_mat(box);
	  for (m=0; m<DIM; m++)
	    box[m][m] = newbox[m];
	}
	
	if (bTDump) {
	  /* determine timestep */
	  if (t0 == -1) {
	    t0=t;
	  } else {
	    if (!bDTset) {
	      dt=t-t0;
	      bDTset=TRUE;
	    }
	  }
	  bDumpFrame = (t >= tdump-(0.5*dt) ) && (t <= tdump+(0.5*dt) );
	}
	
	/* determine if an atom jumped across the box and reset it if so */
	if (bNoJump && (bTPS || frame!=0)) {
	  for(i=0; (i<natoms); i++)
	    for(d=0; (d<DIM); d++)
	      if ( x[i][d]-xp[i][d] > 0.5*box[d][d] )
		x[i][d] -= box[d][d];
	      else if ( x[i][d]-xp[i][d] < -0.5*box[d][d] )
		x[i][d] += box[d][d];
	}
	
	if (bPFit) {
	  /* Now modify the coords according to the flags,
	     for normal fit, this is only done for output frames */
	  rm_pbc(&(top.idef),natoms,box,x,x);
	  
	  reset_x(ifit,ind_fit,nout,index,x,w_rls);
	  do_fit(natoms,w_rls,xp,x);
	}
	
	/* store this set of coordinates for future use */
	if (bPFit || bNoJump) {
	  for(i=0; (i<natoms); i++) {
	    copy_rvec(x[i],xp[i]);
	    rvec_inc(x[i],x_shift);
	  }
	}
	
	if ( bCheckDouble && (t<=pt) && !bToldYouOnce ) {
	  fprintf(stderr,"\nRemoving duplicate frame(s) after t=%g "
		  "(t=%g, frame %d)\n",pt,t,frame);
	  bToldYouOnce=TRUE;
	}
	
	if ( ( ( !bTDump && (frame % skip_nr == 0) ) || bDumpFrame  ) &&
	     ( !bCheckDouble || ( bCheckDouble && (t > pt) ) ) ) {
	  
	  /* remember time from this frame */
	  pt = t; 
	  
	  /* calc new time */
	  if (bTimeStep) 
	    t=tzero+frame*timestep;
	  else
	    if (bChangeTime)
	      t=t+tshift;
	  
	  if (bTDump)
	    fprintf(stderr,"\nDumping frame at t= %g\n",t);
	  
	  /* check for writing at each delta_t */
	  bDoIt=(delta_t == 0);
	  if (!bDoIt)
	    bDoIt=bRmod(t-tzero, delta_t);
	  
	  if (bDoIt || bTDump) {
	    /* print sometimes */
	    if (bToldYouOnce) {
	      bToldYouOnce=FALSE;
	      fprintf(stderr,"\nContinue writing frames from t=%g, frame=%d\n",
		      t,outframe);
	    }
	    if ( ((outframe % SKIP) == 0) || (outframe < SKIP) )
	      fprintf(stderr," ->  frame %6d time %8.3f",outframe,t);
	    
	    if (!bPFit) {
	      /* Now modify the coords according to the flags,
		 for PFit we did this already! */
	      if (bPBC)
		rm_pbc(&(top.idef),natoms,box,x,x);
	      
	      if (bFit) {
		reset_x(ifit,ind_fit,nout,index,x,w_rls);
		do_fit(natoms,w_rls,xp,x);
		for(i=0; (i<natoms); i++)
		  rvec_inc(x[i],x_shift);
	      }
	    }
	    
	    if (bCenter) {
	      center_x(x,box,nout,index);
	    } 
	    
	    if (bCopy) {
	      for(i=0; (i<nout); i++) {
		copy_rvec(x[index[i]],xout[i]);
		if (bVels) copy_rvec(v[index[i]],vout[i]);
	      }
	    }
	    
	    /* this should make sure all atoms in output are really inside
	       a rectangular box. Was needed to make movies.
	       Peter Tieleman, Mon Jul 15 1996
	    */
	    if (bInBox)
	      put_atoms_in_box(nout,box,xout);
	    
	    if (opt2parg_bSet("-shift",asize(pa),pa))
	      for(i=0; (i<nout); i++)
		for (d=0; (d<DIM); d++)
		  xout[i][d] += outframe*shift[d];
	    
	    if (bVels) {
	      /* check if we have velocities and/or coordinates,
		 don't ask me why you can have a frame w/o coords !? */
	      bHaveV=FALSE;
	      bHaveX=FALSE;
	      for (i=0; (i<natoms); i++)
		for (d=0; (d<DIM); d++) {
		  bHaveV=bHaveV || v[i][d];
		  bHaveX=bHaveX || x[i][d];
		}
	    }
	    
	    switch(ftp) {
	    case efTRJ:
	    case efTRR:
	      fwrite_trn(trjout,frame,t,0,box,
			 nout,
			 bHaveX ? xout : NULL,
			 bHaveV ? vout : NULL,
			 NULL);
	      
	      break;
	    case efG87:
	      write_gms(out,nout,xout,bBox ? box : NULL );
	      break;
	    case efXTC:
	      write_xtc(xdr,nout,frame,t,box,xout,xtcpr);
	      break;
	    case efGRO:
	    case efG96:
	    case efPDB:
	      sprintf(title,"Generated by trjconv : %s t= %9.5f",top_title,t);
	      if (bSeparate) {
		sprintf(out_file2,"%d_%s",file_nr,out_file);
		out=ffopen(out_file2,"w");
	      }
	      switch(ftp) {
	      case efGRO:
		write_hconf_p(out,title,&useatoms,prec,xout,bHaveV?vout:NULL,box);
		break;
	      case efPDB:
		fprintf(out,"REMARK    GENERATED BY TRJCONV\n");
		sprintf(title,"%s t= %9.5f",top_title,t);
		write_pdbfile(out,title,&useatoms,xout,box,0,TRUE);
		break;
	      case efG96:
		if (bSeparate || bTDump)
		  write_g96_conf(out,title,&useatoms,xout,bHaveV?vout:NULL,box,
				 0,NULL);
		else {
		  if (outframe == 0)
		    fprintf(out,"TITLE\n%s\nEND\n",title);
		  fprintf(out,"TIMESTEP\n%9d%15.9f\nEND\nPOSITIONRED\n",
			  frame,t);
		  for(i=0; i<nout; i++)
		    fprintf(out,"%15.9f%15.9f%15.9f\n",
			    xout[i][XX],xout[i][YY],xout[i][ZZ]);
		  fprintf(out,"END\nBOX\n%15.9f%15.9f%15.9f\nEND\n",
			  box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
		}
		break;
	      }
	      if (bSeparate) {
		ffclose(out);
		out = NULL;
		file_nr++;
	      }
	      break;
	    default:
	      fatal_error(0,"DHE, ftp=%d\n",ftp);
	    }
	    
	    /* execute command */
	    if (bExec) {
	      char c[255];
	      sprintf(c,"%s  %d",exec_command,file_nr-1);
	      /*fprintf(stderr,"Executing '%s'\n",c);*/
	      system(c);
	    }
	    outframe++;
	  }
	}
	frame++;
	if (bVels) {
	  bHaveNextFrame=read_next_x_or_v(status,&t,natoms,x,v,box);
	} else
	  bHaveNextFrame=read_next_x     (status,&t,natoms,x,  box);
      } while (!bDumpFrame && bHaveNextFrame);
      
      if ( bTDump && !bDumpFrame )
	fprintf(stderr,"\nWARNING no frame found near specified time (%g), "
		"trajectory ended at %g\n",tdump,t);
      fprintf(stderr,"\n");
      
      close_trj(status);
      
      if (ftp == efXTC)
	close_xtc(xdr);
      else if (out != NULL)
	fclose(out);
    }
  }
  thanx(stdout);
  
  return 0;
}
