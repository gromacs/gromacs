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
#include "gstat.h"
#include "magic.h"
#include "binio.h"
#include "pbc.h"

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
  int          j,nre,step;
  real         t,lambda;
  
  if (t0 == -1)
    fatal_error(0,"You forgot to set the truncation time");
  
  /* Check whether this is a .trj file */
  check_trn(fn);
  
  in   = open_trn(fn,"r");
  fp   = fio_getfp(in);
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
	    "frame %d, time %g, bytes %d ??? (type YES if so)\n",
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
    "[BB]1.[bb] from one format to another[BR]",
    "[BB]2.[bb] select a subset of atoms[BR]",
    "[BB]3.[bb] remove periodicity from molecules[BR]",
    "[BB]4.[bb] keep multimeric molecules together[BR]",
    "[BB]5.[bb] center atoms in the box[BR]",
    "[BB]6.[bb] fit atoms to reference structure[BR]",
    "[BB]7.[bb] remove duplicate frames[BR]",
    "[BB]8.[bb] reduce the number of frames[BR]",
    "[BB]9.[bb] change the timestamps of the frames (e.g. t0 and delta-t)",
    "[PAR]",
    "Currently six formats are supported for input and output:",
    "[TT].xtc[tt], [TT].trr[tt], [TT].trj[tt], [TT].gro[tt], [TT].pdb[tt] and",
    "[TT].g87[tt].",
    "The file formats are detected from the file extension.",
    "For the [TT].pdb[tt] and [TT].gro[tt] files compression may be used",
    "using the regular UNIX compress command (this assumes the program ",
    "[TT]compress[tt] to be in your path which might not always be the ",
    "case). For [TT].gro[tt] and [TT].xtc[tt] files the output precision ",
    "can be given as a number of ",
    "decimal places. Note that velocities are only supported in ",
    "[TT].trr[tt], [TT].trj[tt] and [TT].gro[tt] files.[PAR]",
    "In the case of [TT].pdb[tt] or [TT].gro[tt] for output each frame is",
    "written to a separate file. The option [TT]-app[tt] can be used to",
    "append output to an existing trajectory file and,",
    "in case of .pdb or .gro, to write one file with all frames concatenated",
    "(you can use [TT]rasmol -nmrpdb[tt] to view such a .pdb file).", 
    "No checks are made to ensure integrity",
    "of the resulting combined trajectory file.[PAR]",
    "The program is supposed to be useful for making movies with graphics ",
    "programs like Grasp, Quanta or Molscript.[PAR]",
    "It is possible to select part of your trajectory and write it out",
    "to a new trajectory file in order to save disk space, e.g. for leaving",
    "out the water from a trajectory of a protein in water.",
    "[BB]ALWAYS[bb] put the original trajectory on tape!",
    "We recommend to use the portable [TT].xtc[tt] format for your analysis",
    "to save disk space and to have portable files.[PAR]",
    "There are two options for fitting the trajectory to a reference",
    "either for essential dynamics analysis or for whatever.",
    "The first option is just plain fitting to a reference structure",
    "in the run input file, the second option is a progressive fit",
    "in which the first timeframe is fitted to the reference structure ",
    "in the run input file to obtain and each subsequent timeframe is ",
    "fitted to the previously fitted structure. This way a continuous",
    "trajectory is generated, which might not be the case when using the",
    "regular fit method, e.g. when your protein undergoes large",
    "conformational transitions.[PAR]",
    "The option [TT]-removejump[tt] checks if atoms jump across",
    "the box and then puts them back. This has the effect that all molecules",
    "will remain whole (providing they were whole in the initial",
    "conformation), note that this ensures a continuous trajectory but",
    "molecules may (probably will) diffuse out of the box. Use",
    "[TT]-center[tt] to put the system in the center of the box.",
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
    "one specific time from your trajectory.[PAR]"
  };
  
  static bool  bPBC=FALSE,bNoJump=FALSE,bInBox=FALSE,bAppend=FALSE,bVels=TRUE;
  static bool  bCenter=FALSE,bCompress=FALSE;
  static bool  bFit=FALSE,bIFit=FALSE,bBox=TRUE;
  static bool  bCheckDouble=FALSE;
  static int   skip_nr=1,prec=3;
  static real  tzero=0.0,delta_t=0.0,timestep=0.0,ttrunc=-1,tdump=-1,toffset=0;
  static rvec  newbox = {0,0,0};
  static real  xshift=0.0;
  static char  *exec_command=NULL;
  
  t_pargs pa[] = {
    { "-inbox", FALSE,  etBOOL, &bInBox,
      "make sure all atoms are inside box" },
    { "-pbc", FALSE,  etBOOL, &bPBC,
      "make sure molecules are not broken into parts" },
    { "-removejump",FALSE,  etBOOL, &bNoJump,
      "make sure atoms don't jump across the box" },
    { "-center", FALSE,  etBOOL, &bCenter,
      "center atoms in box" },
    { "-xshift", FALSE, etREAL, &xshift,
      "all coordinates will be shifted by framenr*xshift" },
    { "-box", FALSE, etRVEC, &newbox,
      "size for new cubic box (default: read from input)" },
    { "-z", FALSE,  etBOOL, &bCompress,
      "compress output (for .pdb and .gro files)" },
    { "-fit", FALSE,  etBOOL, &bFit,
      "fit molecule to ref structure in .tpx file" },
    { "-pfit", FALSE,  etBOOL, &bIFit,
      "progressive fit, to the previous fitted structure" },
    { "-prec", FALSE,  etINT, &prec,
      "precision for .gro and .xtc writing" },
    { "-vel", FALSE, etBOOL, &bVels,
      "read and write velocities if possible" },
    { "-skip", FALSE,  etINT, &skip_nr,
      "only write out every nr-th frame" },
    { "-dt", FALSE,  etREAL, &delta_t,
      "only write out frame when (t MOD delta_t) == offset" },
    { "-offset", FALSE, etREAL, &toffset,
      "time offset for -dt option" },
    { "-t0", FALSE,  etREAL, &tzero,
      "starting time for trajectory"
      "(default: don't change)"},
#ifndef _win_
    { "-trunc", FALSE, etREAL, &ttrunc,
      "truncate input trj file after this amount of ps" },
#endif
    { "-dump", FALSE, etREAL, &tdump,
      "dump frame nearest specified time" },
    { "-g87box", FALSE,  etBOOL, &bBox,
      "write a box for .g87" },
    { "-exec", FALSE,  etSTR, &exec_command,
      "execute command for every output frame with the frame number "
      "as argument" },
    { "-timestep", FALSE,  etREAL, &timestep,
      "change time step between frames" },
    { "-app", FALSE,  etBOOL, &bAppend,
      "append output"},
    { "-checkdouble", FALSE, etBOOL, &bCheckDouble,
      "only write frames with time larger than previous frame" }
  };
      
  FILE         *out=NULL;
  int          trjout;
  int          status,ftp,ftpin,file_nr;
  rvec         *x,*xn,*xout,*v,*vn,*vout;
  rvec         *xp,shift;
  real         xtcpr, lambda,*w_rls;
  matrix       box;
  int          m,i,d,frame,outframe,natoms,nout,nre,step;
#define SKIP 10
  t_topology   *top=NULL;
  t_atoms      *atoms=NULL,useatoms;
  int          isize;
  atom_id      *index;
  char         *grpname;
  int          ifit,irms;
  atom_id      *ind_fit,*ind_rms;
  char         *gn_fit,*gn_rms;
  real         t,pt,tshift,t0=-1,dt=0.001;
  bool         bSelect,bDoIt,bIndex,bTDump,bSetTime,bTop=FALSE,bDTset=FALSE;
  bool         bExec,bTimeStep=FALSE,bDumpFrame=FALSE,bToldYouOnce=FALSE;
  bool         bHaveNextFrame,bHaveX,bHaveV,bSetBox;
  char         *grpnm;
  char         *top_file,*in_file,*out_file,out_file2[256];
  char         top_title[256],title[256],command[256],filemode[5];
  int          xdr;

  t_filenm fnm[] = {
    { efTRX, "-f",  NULL, ffREAD },
    { efTRX, "-o", "trajout", ffWRITE },
    { efTPS, NULL,  NULL, ffOPTRD },
    { efNDX, NULL,  NULL, ffOPTRD }
  };  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,FALSE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

  top_file=ftp2fn(efTPS,NFILE,fnm);

  /* Check command line */
  in_file=opt2fn("-f",NFILE,fnm);
  if (ttrunc != -1) {
#ifndef _win_
    do_trunc(in_file,ttrunc);
#endif
  }
  else {
    if (bIFit) {
      bFit=TRUE;
    }
    bSetBox   = opt2parg_bSet("-box", asize(pa), pa);
    bSetTime  = opt2parg_bSet("-t0", asize(pa), pa);
    bExec     = opt2parg_bSet("-exec", asize(pa), pa);
    bTimeStep = opt2parg_bSet("-timestep", asize(pa), pa);
    bTDump    = opt2parg_bSet("-dump", asize(pa), pa);
    bPBC = bPBC || bFit;
    if (bPBC && !fn_bTPX(top_file)) {
      fprintf(stderr,
	      "WARNING: can not remove periodicity without a run input file\n");
      bPBC=FALSE;
    }
    if (bNoJump && bPBC) {
      fprintf(stderr,
	      "WARNING: both -pbc and -removejump specified: ignoring -pbc\n");
      bPBC=FALSE;
    }
    /* prec is in nr of decimal places, xtcprec is a multiplication factor: */
    xtcpr=1;
    for (i=0; i<prec; i++)
      xtcpr*=10;

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
    }
    
    /* skipping */  
    if (skip_nr <= 0) {
      fprintf(stderr,"The number of frames to skip is <= 0. "
	      "Writing out all frames.\n\n");
      skip_nr = 1;
    } 
    
    /* Check whether we really should compress */
    bCompress = bCompress && ((ftp == efGRO) || (ftp == efPDB));
    
    /* Determine whether to read a topology */
    bTop = (ftp2bSet(efTPS,NFILE,fnm) || 
	    bPBC || bFit || (ftp == efGRO) || (ftp == efPDB));

    /* Determine if when can read index groups */
    bIndex = (bIndex || bTop);
     
    if (bTop)
      read_tps_conf(top_file,title,top,&atoms,&xp,NULL,box,bFit);

    if (bFit) {
      fprintf(stderr,"Select group for root least squares fit\n");
      get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),
		1,&ifit,&ind_fit,&gn_fit);

      if (ifit < 3) 
	fatal_error(0,"Need at least 3 points to fit!\n");
    }
    
    if (bIndex) {
      fprintf(stderr,"Select group for output\n");
      get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),
		1,&isize,&index,&grpnm);
    }
    else {
      /* no index file, so read natoms from TRX */
      natoms=read_first_x(&status,in_file,&t,&x,box);
      close_trj(status);
      snew(index,natoms);
      for(i=0;i<natoms;i++)
	index[i]=i;
      isize=natoms; 
    }
    
    /* if xp was not snew-ed before, do it now */
    if (!xp)
      snew(xp, natoms);
    
    if (bFit) {
      snew(w_rls,atoms->nr);
      for(i=0; (i<ifit); i++)
	w_rls[ind_fit[i]]=atoms->atom[ind_fit[i]].m;
      
      /* Restore reference structure and set to origin, 
         store original location (to put structure back) */
      if (bPBC)
	rm_pbc(&(top->idef),atoms->nr,box,xp,xp);
      copy_rvec(xp[index[0]],shift);
      reset_x(ifit,ind_fit,isize,index,xp,w_rls);
      rvec_dec(shift,xp[index[0]]);
    }
    
    /* Make atoms struct for output in GRO or PDB files */
    if ((ftp == efGRO) || (ftp == efPDB)) {
      /* get memory for stuff to go in pdb file */
      init_t_atoms(&useatoms,atoms->nr,FALSE);
      sfree(useatoms.resname);
      useatoms.resname=atoms->resname;
      for(i=0;(i<isize);i++) {
	useatoms.atomname[i]=atoms->atomname[index[i]];
	useatoms.atom[i].resnr=atoms->atom[index[i]].resnr;
	useatoms.nres=max(useatoms.nres,useatoms.atom[i].resnr+1);
      }
      useatoms.nr=isize;
    }
    /* open trj file for reading */
    if (bVels)
      natoms=read_first_x_or_v(&status,in_file,&t,&x,&v,box);
    else
      natoms=read_first_x     (&status,in_file,&t,&x,   box);
    if (bSetTime)
      tshift=tzero-t;
    else
      tzero=t;
    
    /* check if index is meaningful */
    if (isize > natoms)
      fatal_error(0,"Index has more atoms (%d) than trajectory file (%d)",
		  isize,natoms);
    else {
      bSelect=(natoms > isize);
      for(i=0; ((i<isize) && !bSelect); i++)
	bSelect=(i != index[i]);
      if (bSelect) {
	snew(xn,isize);
	xout=xn;
	snew(vn,isize);
	vout=vn;
	nout=isize;
      } else {
	xout=x;
	vout=v;
	nout=natoms;
      }
    }
    
    if (ftp == efG87)
      fprintf(out,"Generated by %s. #atoms=%d, %s BOX is stored in "
	      "this file.\n",Program(),nout,bBox ? "a" : "NO");
    
    /* Start the big loop over frames */
    file_nr=0;  
    frame=0;
    outframe=0;
    pt=-666.0;
    bDTset=FALSE;
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
      /* don't do this for first frame! */
      if (bNoJump && (frame!=0)) {
	for(i=0; (i<natoms); i++)
	  for(d=0; (d<DIM); d++)
	    if ( x[i][d]-xp[i][d] > 0.5*box[d][d] )
	      x[i][d]-=box[d][d];
	    else if ( x[i][d]-xp[i][d] < -0.5*box[d][d] )
	      x[i][d]+=box[d][d];
      }
      
      if (bIFit) {
	/* Now modify the coords according to the flags,
	   for normal fit, this is only done for output frames */
	if (bPBC)
	  rm_pbc(&(top->idef),natoms,box,x,x);
	
	reset_x(ifit,ind_fit,isize,index,x,w_rls);
	do_fit(natoms,w_rls,xp,x);
      }
      
      /* store this set of coordinates for future use */
      if (bIFit || bNoJump) {
	for(i=0; (i<natoms); i++) {
	  copy_rvec(x[i],xp[i]);
	  rvec_inc(x[i],shift);
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
	  if (bSetTime)
	    t=t+tshift;

	if (bTDump)
	  fprintf(stderr,"\nDumping frame at t= %g\n",t);

	/* check for writing at each delta_t */
	bDoIt=(delta_t == 0) || (t == 0);
	if (!bDoIt)
	  bDoIt=bRmod(t-toffset,delta_t);
	
	if (bDoIt || bTDump) {
	  /* print sometimes */
	  if (bToldYouOnce) {
	    bToldYouOnce=FALSE;
	    fprintf(stderr,"\nContinue writing frames from t=%g, frame=%d\n",
		    t,outframe);
	  }
	  if ( ((outframe % SKIP) == 0) || (outframe < SKIP) )
	    fprintf(stderr," ->  frame %6d time %8.3f",outframe,t);
	  
	  if (!bIFit) {
	    /* Now modify the coords according to the flags,
	       for IFit we did this already! */
	    if (bPBC) 
	      rm_pbc(&(top->idef),natoms,box,x,x);
	  
	    if (bFit) {
	      reset_x(ifit,ind_fit,isize,index,x,w_rls);
	      do_fit(natoms,w_rls,xp,x);
	      for(i=0; (i<natoms); i++)
		rvec_inc(x[i],shift);
	    }
	  }
	  
	  if (bCenter) {
	    center_x(x,box,isize,index);
	  } 
	  
	  if (bSelect) {
	    for(i=0; (i<isize); i++) {
	      copy_rvec(x[index[i]],xout[i]);
	      if (bVels) copy_rvec(v[index[i]],vout[i]);
	    }
	  }
	  
	  /* this should make sure all atoms in output are really inside
	     a rectangular box. Was needed to make movies.
	     Peter Tieleman, Mon Jul 15 1996
	     */
	  if (bInBox)
	    put_all_atoms_in_box(isize,box,xout);
	  
	  if (xshift != 0.0)
	    for(i=0; (i<isize); i++)
	      xout[i][XX]+=(outframe*xshift);
	  
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
	  case efGRO: {
	    FILE *fp;
	    
	    sprintf(title,"Generated by trjconv : %s t= %9.3f",top_title,t);
	    if (bAppend) {
	      sprintf(out_file2,"%s\0",out_file);
	      fp=ffopen(out_file2,"a");
	    } else {
	      if (bTDump)
		sprintf(out_file2,"%s\0",out_file);
	      else
		sprintf(out_file2,"%d_%s\0",file_nr,out_file);
	      fp=ffopen(out_file2,"w");
	    }
	    write_hconf_p(fp,title,&useatoms,prec,xout,bHaveV?vout:NULL,box);
	    ffclose(fp);
	    file_nr++;
	    break;
	  }
	  case efPDB: {
	    FILE      *fp;
	    int       nn;
	    
	    if (bAppend) {
	      sprintf(out_file2,"%s\0",out_file);
	      fp=ffopen(out_file2,"a");
	    }
	    else {
	      if (bTDump)
		sprintf(out_file2,"%s\0",out_file);
	      else
		sprintf(out_file2,"%d_%s\0",file_nr,out_file);
	      fp=ffopen(out_file2,"w");
	    }
	    fprintf(fp,"REMARK    GENERATED BY TRJCONV\n");
	    sprintf(title,"%s t= %9.3f",top_title,t);
	    write_pdbfile(fp,title,&useatoms,x,box,0,TRUE);
	    ffclose(fp);
	    
	    file_nr++;
	    
	    break;
	  }
	  default:
	    fatal_error(0,"DHE, ftp=%d\n",ftp);
	  }
	  /* trying to compress file using compress. This is by no means full- 
	   * proof, but feel free to add searches for environment variables, 
	   * alternative compressing programs (gzip) and more checks. 
	   */
	  if (bCompress) {
	    sprintf(command,"compress -f %s",out_file);
	    fprintf(stderr,"\rCompressing %s",out_file);
	    system(command);
	  }
	  
	  /* execute command */
	  if ( bExec ) {
	    char c[255];
	    sprintf(c,"%s %d",exec_command,file_nr - 1);
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
  
  thanx(stdout);
  
  return 0;
}
