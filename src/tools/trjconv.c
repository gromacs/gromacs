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


/*

void do_patch(char *fn)
{
  FILE         *stat;
  int          magic=GROMACS_MAGIC;
  t_statheader sh;
  long         pos;
  int          n,teller;
  real         r;

  is_trn(fn);
  stat=ffopen(fn,"r+");
  teller=0;
  while (!eof(stat)) {
    pos=ftell(stat);
    blockwrite(stat,magic);
    fseek(stat,pos,SEEK_SET);
    rd_header(stat,&sh);
    rd_hstatus(stat,&sh,&n,&r,&r,NULL,NULL,NULL,NULL,
	       &n,NULL,NULL,NULL,&n,NULL,NULL);
    teller++;
  }
  fclose(stat);
  fprintf(stderr,"\n");
}
*/

void check_trn(char *fn)
{
  if ((fn2ftp(fn) != efTRJ)  || (fn2ftp(fn) != efTRR))
    fatal_error(0,"%s is not a trj file, exiting\n",fn);
}

void do_trunc(char *fn, real t0)
{
  int          in;
  FILE         *fp;
  bool         bStop;
  t_trnheader  sh;
  long         fpos;
  char         yesno[256];
  int          j,natoms,nre,step;
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
  while (!bStop && fread_trnheader(in,&sh)) {
    fread_htrn(in,&sh,NULL,NULL,NULL,NULL);
    
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

char *ftpname(int ftp)
{
  static char name[10];
  
  switch (ftp) {
  case efTRJ: 
    strcpy(name,"trj");
    break;
  case efG87: 
    strcpy(name,"gromos 87");
    break;
  case efXTC: 
    strcpy(name,"xtc");
    break;
  case efGRO: 
    strcpy(name,"gro");
    break;
  case efPDB: 
    strcpy(name,"pdb");
    break;
  default:
    strcpy(name,"unknown");
  }
  
  return name;
}

bool check_ftp(int ftp, int nftp)
{
  char name[10], nname[10];
  
  if (ftp!=efTRJ) {
    strcpy( name,ftpname( ftp));
    strcpy(nname,ftpname(nftp));
    fprintf(stderr,"WARNING: multiple output filetypes set (%s and %s), "
	    "will use %s\n",name,nname,name);
    return FALSE;
  }
  return TRUE;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "trjconv can convert trajectory files in many ways:[BR]",
    "[BB]1.[bb] from one format to another[BR]",
    "[BB]2.[bb] select a subset of atoms[BR]",
    "[BB]3.[bb] remove periodicity from molecules[BR]",
    "[BB]4.[bb] center atoms in the box[BR]",
    "[BB]5.[bb] fit atoms to reference structure[BR]",
    "[BB]6.[bb] remove duplicate frames[BR]",
    "[BB]7.[bb] reduce the number of frames[BR]",
    "[BB]8.[bb] change the timestamps of the frames (e.g. t0 and delta-t)[BR]",
    "Currently five formats are supported for input and output:",
    "[TT].trj[tt], [TT].xtc[tt], [TT].g87[tt], [TT].pdb[tt], [TT].gro[tt].",
    "On input the file format is detected from the file extension,",
    "on output a [TT].trj[tt] file is generated unless one of the other file",
    "types is explicitly requested by a command line option.",
    "For the [TT].pdb[tt] and [TT].gro[tt] files compression may be used",
    "using the regular UNIX compress command (this assumes the program ",
    "[TT]compress[tt] to be in your path which might not always be the ",
    "case). For [TT].gro[tt] and [TT].xtc[tt] files the output precision ",
    "can be given as a number of ",
    "decimal places. Note that velocities are only supported in ",
    "[TT].trj[tt] and [TT].gro[tt] files.[PAR]",
    "In the case of [TT].pdb[tt] or [TT].gro[tt] for output each frame is",
    "written to a separate file. The option [TT]-app[tt] can be used to",
    "write one pdb file with all frames concatenated, or to append output",
    "to an existing trajectory file. No checks are made to ensure integrity",
    "of the resulting combined trajectory file.[PAR]",
    "The program is supposed to be useful",
    "for making movies with graphics programs like Grasp or Quanta.[PAR]",
    "It is possible to select part of your trajectory and write it out",
    "to a new trajectory file in order to save disk space, e.g. for leaving",
    "out the water from a trajectory of a protein in water.",
    "[BB]ALLWAYS[bb] put the original trajectory on tape!",
    "We recommend to use the portable [TT].xtc[tt] format for your analysis",
    "to save disk space and to have portable files.[PAR]",
    "There are two options for fitting the trajectory to a reference",
    "either for essential dynamics analysis or for whatever.",
    "The first option is just plain fitting to a reference structure",
    "in the [TT].tpx[tt] file, the second option is a progressive fit",
    "in which the first timeframe is fitted to the reference structure ",
    "in the [TT].tpx[tt] file to obtain and each subsequent timeframe is ",
    "fitted to the previously fitted structure. This way a continuous",
    "trajectory is generated, which might not be the case using the",
    "regular fit method, e.g. when your protein undergoes large",
    "conformational transitions.[PAR]",
    "With the option [TT]-dt[tt] it is possible to reduce the number of ",
    "frames in the output. This option relies on the accuracy of the times ",
    "in your input trajectory, so if these are inaccurate use the -timestep ",
    "option to modify the time (this can be done simultaneously).[PAR]",
    "Using [TT]-trunc[tt] trjconv can truncate [TT].trj[tt] in place, i.e. ",
    "without copying the file. This is useful when a run has crashed",
    "during disk I/O (one more disk full), or when two contiguous",
    "trajectories must be concatenated without have double frames.[PAR]",
    "Also the option [TT]-checkdouble[tt] may be used to remove all",
    "duplicate frames from such a concatenated trajectory, this is done",
    "by ignoring all frames with a time smaller than or equal to the previous",
    "frame.[PAR]",
    "The option [TT]-dump[tt] can be used to extract a frame at or near",
    "one specific time from your trajectory.[PAR]",
    "Finally, trjconv can patch the magic number of [TT].trj[tt] files",
    "in case a new version of [BB]GROMACS[bb] requires this."
  };
  
  static bool  bPBC=FALSE,bInBox=FALSE,bAppend=FALSE,bVels=TRUE;
  static bool  bCenter=FALSE,bCompress=FALSE;
  static bool  bFit=FALSE,bIFit=FALSE,bBox=TRUE;
  static bool  bTer=FALSE,/*bPatch=FALSE,*/bCheckDouble=FALSE;
  static int   skip_nr=1,prec=3;
  static real  tzero=0.0,delta_t=0.0,timestep=0.0,ttrunc=-1,tdump=-1,toffset=0;
  static real  newbox = -1, xshift=0.0;
  static char  *exec_command;
  
  t_pargs pa[] = {
    { "-inbox", FALSE,  etBOOL, &bInBox,
      "make sure all atoms are inside box." },
    { "-pbc", FALSE,  etBOOL, &bPBC,
      "make sure molecules are not broken into parts in output." },
    { "-center", FALSE,  etBOOL, &bCenter,
      "center atoms in box." },
    { "-xshift", FALSE, etREAL, &xshift,
      "all coordinates will be shifted by framenr*xshift" },
    { "-box", FALSE, etREAL, &newbox,
      "size for new cubic box (default: read box from input)" },
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
    /*{ "-patch", FALSE, etBOOL, &bPatch,
      "Set the new magic number in a trj file" },*/
    { "-skip", FALSE,  etINT, &skip_nr,
      "only write out every nr-th frame" },
    { "-dt", FALSE,  etREAL, &delta_t,
      "only write out frame when (t MOD delta_t) == offset." },
    { "-offset", FALSE, etREAL, &toffset,
      "time offset for -dt option" },
    { "-t0", FALSE,  etREAL, &tzero,
      "change times in trajectory by specified amount "
      "(default: don't change)"},
    { "-trunc", FALSE, etREAL, &ttrunc,
      "truncate input trj file after this amount of ps" },
    { "-dump", FALSE, etREAL, &tdump,
      "dump frame nearest specified time" },
    { "-g87box", FALSE,  etBOOL, &bBox,
      "write a box for .g87" },
    { "-exec", FALSE,  etSTR, &exec_command,
      "execute command every time frame with the frame number as argument" },
    { "-timestep", FALSE,  etREAL, &timestep,
      "change time step between frames (default: don't change)" },
    { "-app", FALSE,  etBOOL, &bAppend,
      "append output"},
    { "-ter", FALSE, etBOOL, &bTer,
      "Use 'TER' in pdb file as end of frame in stead of the "
      "default 'ENDMDL'" },
    { "-checkdouble", FALSE, etBOOL, &bCheckDouble,
      "only write frames with time larger than previous frame" }
  };
      
  FILE         *out=NULL;
  int          trjout;
  int          status,ftp,ftpout,file_nr;
  rvec         *x,*xn,*xout,*v;
  rvec         *xp,shift;
  real         xtcpr, lambda,*w_rls;
  matrix       box;
  int          m,i,d,frame,outframe,natoms,nout,nre,step;
  t_tpxheader  header;
  t_topology   top;
  t_atoms      useatoms;
  int          isize;
  atom_id      *index;
  char         *grpname;
  int          ifit,irms;
  atom_id      *ind_fit,*ind_rms;
  char         *gn_fit,*gn_rms;
  real         t,pt,t0=-1,dt=0.001;
  bool         bSelect,bDoIt,bIndex,bTDump,bSetTime,bTop=FALSE,bDTset=FALSE;
  bool         bExec,bTimeStep=FALSE,bDumpFrame=FALSE,bToldYouOnce=FALSE;
  bool         bHaveNextFrame,bHaveX,bHaveV;
  char         *grpnm;
  char         title[256],out_file[256],command[256],filemode[5];
  int          xdr;

  t_filenm fnm[] = {
    { efTRX, "-f",  NULL, ffREAD },
    { efTRJ, "-ot", "trajout", ffOPTWR },
    { efXTC, "-ox", "ctrajout", ffOPTWR },
    { efG87, "-og", "gtrajout", ffOPTWR },
    { efGRO, "-or", "confout", ffOPTWR },
    { efPDB, "-op", NULL, ffOPTWR },
    { efTPX, NULL,  NULL, ffOPTRD },
    { efNDX, NULL,  NULL, ffOPTRD }
  };  
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,FALSE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);

  /* Check command line */
  if (ttrunc != -1) {
    do_trunc(ftp2fn(efTRX,NFILE,fnm),ttrunc);
  }
  /*else if (bPatch) {
    do_patch(ftp2fn(efTRX,NFILE,fnm));
  }*/
  else {
    if (bIFit) {
      bFit=TRUE;
    }
    bSetTime  = opt2parg_bSet("-t0", asize(pa), pa);
    bExec     = opt2parg_bSet("-exec", asize(pa), pa);
    bTimeStep = opt2parg_bSet("-timestep", asize(pa), pa);
    bTDump    = opt2parg_bSet("-dump", asize(pa), pa);
    pdb_use_ter(bTer);
    xtcpr=1;
    for (i=0; i<prec; i++)
      xtcpr*=10;

    bIndex=ftp2bSet(efNDX,NFILE,fnm);
    
    /* Determine output type */  
    ftp=efTRJ;
    if (ftp2bSet(efXTC,NFILE,fnm))
      if (check_ftp(ftp,efXTC))
	ftp=efXTC;
    if (ftp2bSet(efPDB,NFILE,fnm))
      if (check_ftp(ftp,efPDB))
	ftp=efPDB;
    if (ftp2bSet(efGRO,NFILE,fnm))
      if (check_ftp(ftp,efGRO))
	ftp=efGRO;
    if (ftp2bSet(efG87,NFILE,fnm))
      if (check_ftp(ftp,efG87))
	ftp=efG87;
    if (bVels) {
      /* check if velocities are possible in input and output files */
      ftpout=fn2ftp(ftp2fn(efTRX,NFILE,fnm));
      bVels=((ftp==efTRJ) || (ftp==efGRO)) 
	&& ((ftpout==efTRJ) || (ftpout==efGRO));
    } else {
      bHaveX=TRUE;
      bHaveV=FALSE;
    }
    
    if (bAppend) {
      strcpy(filemode,"a");
      if (fexist(ftp2fn_null(ftp,NFILE,fnm))) {
	fprintf(stderr,"APPENDING to existing file %s\n",
		ftp2fn(ftp,NFILE,fnm));
      }
    } else {
      strcpy(filemode,"w");
    }
    switch (ftp) {
    case efXTC:
      xdr = open_xtc(ftp2fn(ftp,NFILE,fnm),filemode);
      break;
    case efG87:
      out=ftp2FILE(ftp,NFILE,fnm,filemode);
      break;
    case efTRJ:
      out=NULL;
      trjout = open_tpx(ftp2fn(ftp,NFILE,fnm),filemode);
      break;
    }
    
    /* skipping */  
    if (skip_nr <= 0) {
      fprintf(stderr,"The number of files to skip is <= 0. "
	      "Writing out all frames.\n\n");
      skip_nr = 1;
    } 
    
    /* Check whether we really should compress */
    bCompress = bCompress && ((ftp == efGRO) || (ftp == efPDB));
    
    /* Determine whether to read a topology */
    bTop = (ftp2bSet(efTPX,NFILE,fnm) || 
	    bPBC || bFit || (ftp == efGRO) || (ftp == efPDB));

    /* Determine if when can read index groups */
    bIndex = (bIndex || bTop);
     
    if (bTop) {  
      read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header);
      snew(xp,header.natoms);
      read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,NULL,box,
	       &natoms,xp,NULL,NULL,&top);
    }
    if (bFit) {
      fprintf(stderr,"Select group for root least squares fit\n");
      get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
		1,&ifit,&ind_fit,&gn_fit);

      if (ifit < 3) 
	fatal_error(0,"Need >= 3 points to fit!\n");
    }  

    if (bIndex) {
      fprintf(stderr,"Select group for output\n");
      get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
		1,&isize,&index,&grpnm);
    }
    else {
      /* no index file, so read natoms from TRX */
      natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
      close_trj(status);
      snew(index,natoms);
      for(i=0;i<natoms;i++)
	index[i]=i;
      isize=natoms; 
    }

    if (bFit) {
      snew(w_rls,header.natoms);
      for(i=0; (i<ifit); i++)
	w_rls[ind_fit[i]]=top.atoms.atom[ind_fit[i]].m;
      
      /* Restore reference structure and set to origin, 
         store original location (to put structure back) */
      rm_pbc(&(top.idef),top.atoms.nr,box,xp,xp);
      clear_rvec(shift);
      for (i=0; (i<header.natoms); i++)
	rvec_inc(shift,xp[i]);
      svmul(1./header.natoms,shift,shift);
      reset_x(ifit,ind_fit,isize,index,xp,w_rls);
    }
    
    /* Make atoms struct for output in GRO or PDB files */
    if ((ftp == efGRO) || (ftp == efPDB)) {
      /* get memory for stuff to go in pdb file */
      snew(useatoms.atom,top.atoms.nr);
      snew(useatoms.atomname,top.atoms.nr);
      
      useatoms.nres=0;
      useatoms.resname=top.atoms.resname;
      for(i=0;(i<isize);i++) {
	useatoms.atomname[i]=top.atoms.atomname[index[i]];
	useatoms.atom[i].resnr=top.atoms.atom[index[i]].resnr;
	useatoms.nres=max(useatoms.nres,useatoms.atom[i].resnr+1);
      }
      useatoms.nr=isize;
    }
    /* open trj file for reading */
    if (bVels)
      natoms=read_first_x_v(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,&v,box);
    else
      natoms=read_first_x  (&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,   box);
    if (!bSetTime)
      tzero=t;
    
    /* check if index is meaningful */
    if (isize > natoms)
      fatal_error(0,"Index has more atoms (%d) than trajectory file (%d)",
		  isize,natoms);
    else {
      bSelect=(natoms != isize);
      for(i=0; ((i<isize) && !bSelect); i++)
	bSelect=(i != index[i]);
      if (bSelect) {
	snew(xn,isize);
	xout=xn;
	nout=isize;
      } else {
	xout=x;
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
      if (newbox >= 0) {
	clear_mat(box);
	for (m=0; m<DIM; m++)
	  box[m][m] = newbox;
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
      
      if (bIFit) {
	/* Now modify the coords according to the flags,
	   for normal fit, this is only done for output frames */
	rm_pbc(&(top.idef),natoms,box,x,x);
	
	reset_x(ifit,ind_fit,isize,index,x,w_rls);
	do_fit(natoms,w_rls,xp,x);
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
	    t=t+tzero;

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
	  if ((outframe % 10) == 0)
	    fprintf(stderr,"->  frame %6d time %8.3f",outframe,t);
	 
	  if (!bIFit) {
	    /* Now modify the coords according to the flags,
	       for IFit we did this already! */
	    if (bPBC || bFit) 
	      rm_pbc(&(top.idef),natoms,box,x,x);
	  
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
	    for(i=0; (i<isize); i++)
	      copy_rvec(x[index[i]],xn[i]);
	  }
	  
	  /* this should make sure all atoms in output are really inside
	     a rectangular box. Was needed to make movies.
	     Peter Tieleman, Mon Jul 15 1996
	     */
	  if (bInBox) {
	    for(i=0; (i<isize); i++)
	      for(m=0; m < DIM; m++)
		{
		  if (xout[i][m] < 0) xout[i][m] += box[m][m];
		  if (xout[i][m] > box[m][m]) xout[i][m] -= box[m][m];
		}
	  } 
	  
	  if (xshift != 0.0)
	    for(i=0; (i<isize); i++)
	      xout[i][XX]+=(outframe*xshift);
	  
	  if (bVels) {
	  /* check if we have velocities and/or coordinates,
	     don't ask me why you can have a frame w/o coords !? */
	    bHaveV=FALSE;
	    for (i=0; (i<natoms); i++)
	      for (d=0; (d<DIM); d++)
		bHaveV=bHaveV || v[i][d];
	    bHaveX=FALSE;
	    for (i=0; (i<natoms); i++)
	      for (d=0; (d<DIM); d++)
		bHaveX=bHaveX || x[i][d];
	  }
	  
	  switch(ftp) {
	  case efTRJ:
	  case efTRR:
	    fwrite_trn(trjout,frame,t,0,box,
		       nout,
		       bHaveX ? xout : NULL,
		       bHaveV ? v    : NULL,
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
	    
	    sprintf(title,"Generated by trjconv : %s t= %9.3f",*(top.name),t);
	    if (bAppend) {
	      sprintf(out_file,"%s\0",ftp2fn(efGRO,NFILE,fnm));
	      fp=ffopen(out_file,"a");
	    } else {
	      if (bTDump)
		sprintf(out_file,"%s\0",ftp2fn(efGRO,NFILE,fnm));
	      else
		sprintf(out_file,"%d_%s\0",file_nr,ftp2fn(efGRO,NFILE,fnm));
	      fp=ffopen(out_file,"w");
	    }
	    write_hconf_p(fp,title,&useatoms,prec,xout,bHaveV ? v : NULL,box);
	    ffclose(fp);
	    file_nr++;
	    break;
	  }
	  case efPDB: {
	    FILE      *fp;
	    t_pdbatom *pdba;
	    int       nn;
	    char      cc;
	    
	    pdba=atoms2pdba(&useatoms,xout);
	    
	    if (bAppend) {
	      cc='A'+(file_nr % 26); /* counts A..Z,A..Z,etc. */
	      for(nn=0; (nn<useatoms.nr); nn++)
		pdba[nn].chain=cc;
	      
	      sprintf(out_file,"%s\0",ftp2fn(efPDB,NFILE,fnm));
	      fp=ffopen(out_file,"a");
	    }
	    else {
	      if (bTDump)
		sprintf(out_file,"%s\0",ftp2fn(efPDB,NFILE,fnm));
	      else
		sprintf(out_file,"%d_%s\0",file_nr,ftp2fn(efPDB,NFILE,fnm));
	      fp=ffopen(out_file,"w");
	    }
	    fprintf(fp,"REMARK    GENERATED BY TRJCONV : %s t= %9.3f\n",
		    *(top.name),t);
	    print_pdbatoms(fp,useatoms.nr,pdba,box);
	    ffclose(fp);
	    sfree(pdba);
	    
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
