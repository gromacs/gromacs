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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "sysstuff.h"
#include "string2.h"
#include "smalloc.h"
#include "pbc.h"
#include "statutil.h"
#include "gmxfio.h"
#include "trnio.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "gmxfio.h"
#include "xtcio.h"
#include "pdbio.h"
#include "confio.h"
#include "wgms.h"
#include <math.h>

/* defines for frame counter output */
static int __frame=NOTSET;
#define SKIP 10
#define INITCOUNT __frame=-1

/* frames for read_first/next_x */
static t_trxframe *xframe=NULL;
static int nxframe=0;

int nframes_read(void)
{
  return __frame;
}
  
static void printcount_(char *l,real t)
{
  fprintf(stderr,"\r%-14s %6d time %8.3f   ",l,__frame,convert_time(t));
}

static void printcount2(real t,real t0,real tpf,real tppf)
{
  __frame++;
  if (__frame % SKIP == 0 || __frame < SKIP)
    printcount_(check_times2(t,t0,tpf,tppf) < 0 ? "Skipping frame": "Reading frame",t);
}

static void printcount(real t)
{
  printcount2(t,t,t,t);
}

static void printlast(real t)
{
  printcount_("Last frame",t);
  fprintf(stderr,"\n");
}

static void printincomp(t_trxframe *fr)
{
  if (fr->not_ok & HEADER_NOT_OK)
    fprintf(stderr,"WARNING: Incomplete header: nr %d time %g\n",
	    __frame+1,fr->time);
  else if (fr->not_ok)
    fprintf(stderr,"WARNING: Incomplete frame: nr %d time %g\n",
	    __frame+1,fr->time);
}

int prec2ndec(real prec)
{
  if (prec <= 0)
    fatal_error(0,"DEATH HORROR prec (%g) <= 0 in prec2ndec",prec);
  
  return (int)(log(prec)/log(10.0)+0.5);
}

/* Globals for gromos-87 input */
typedef enum { effXYZ, effXYZBox, effG87, effG87Box, effNR } eFileFormat;
static eFileFormat eFF;
static int         NATOMS;
static double      DT,BOX[3];
static bool        bReadBox;

void clear_trxframe(t_trxframe *fr,bool bFirst)
{
  fr->not_ok  = 0;
  fr->bTitle  = FALSE;
  fr->bStep   = FALSE;
  fr->bTime   = FALSE;
  fr->bLambda = FALSE;
  fr->bAtoms  = FALSE;
  fr->bPrec   = FALSE;
  fr->bX      = FALSE;
  fr->bV      = FALSE;
  fr->bF      = FALSE;
  fr->bBox    = FALSE;
  if (bFirst) {
    fr->flags  = 0;
    fr->natoms = -1;
    fr->t0     = 0;
    fr->tpf    = 0;
    fr->tppf   = 0;
    fr->title  = NULL;
    fr->step   = 0;
    fr->time   = 0;
    fr->lambda = 0;
    fr->atoms  = NULL;
    fr->prec   = 0;
    fr->x      = NULL;
    fr->v      = NULL;
    fr->f      = NULL;
    clear_mat(fr->box);
  }
}

int write_trxframe_indexed(int fnum,t_trxframe *fr,int nind,atom_id *ind)
{
  char title[STRLEN];
  rvec *xout=NULL,*vout=NULL,*fout=NULL;
  int  i;
  real prec;

  if (fr->bPrec)
    prec = fr->prec;
  else
    prec = 1000.0;
  
  switch (fio_getftp(fnum)) {
  case efTRJ:
  case efTRR:
    break;
  default:
    if (!fr->bX)
      fatal_error(0,"Need coordinates to write a %s trajectory",
		  ftp2ext(fio_getftp(fnum)));
    break;
  }

  switch (fio_getftp(fnum)) {
  case efTRJ:
  case efTRR:
    if (fr->bV) {
      snew(vout,nind);
      for(i=0; i<nind; i++) 
	copy_rvec(fr->v[ind[i]],vout[i]);
    }
    if (fr->bF) {
      snew(fout,nind);
      for(i=0; i<nind; i++) 
	copy_rvec(fr->f[ind[i]],fout[i]);
    }
  case efXTC:
  case efG87:
    if (fr->bX) {
      snew(xout,nind);
      for(i=0; i<nind; i++) 
	copy_rvec(fr->x[ind[i]],xout[i]);
    }
    break;
  default:
    break;
  }

  switch (fio_getftp(fnum)) {
  case efXTC: 
    write_xtc(fnum,nind,fr->step,fr->time,fr->box,xout,prec);
    break;
  case efTRJ:
  case efTRR:  
    fwrite_trn(fnum,nframes_read(),
	       fr->time,fr->step,fr->box,nind,xout,vout,fout);
    break;
  case efGRO:
  case efPDB:
  case efBRK:
  case efENT:
    if (!fr->bAtoms)
      fatal_error(0,"Can not write a %s file without atom names",
		  ftp2ext(fio_getftp(fnum)));
    sprintf(title,"frame t= %.3f",fr->time);
    if (fio_getftp(fnum) == efGRO)
      write_hconf_indexed_p(fio_getfp(fnum),title,fr->atoms,nind,ind,
			    prec2ndec(prec),
			    fr->x,fr->bV ? fr->v : NULL,fr->box);
    else
      write_pdbfile_indexed(fio_getfp(fnum),title,fr->atoms,
			    fr->x,fr->box,0,fr->step,nind,ind);
    break;
  case efG87:
    write_gms(fio_getfp(fnum),nind,xout,fr->box);
    break;
  case efG96:
    write_g96_conf(fio_getfp(fnum),fr,nind,ind); 
    break;
  default:
    fatal_error(0,"Sorry, write_trxframe_indexed can not write %s",
		ftp2ext(fio_getftp(fnum)));
    break;
  }

  switch (fio_getftp(fnum)) {
  case efTRN:
  case efTRJ:
  case efTRR:
    if (vout) sfree(vout);
    if (fout) sfree(fout);
  case efXTC:
  case efG87:
    sfree(xout);
    break;
  default:
    break;
  }
  
  return 0;
}

int write_trxframe(int fnum,t_trxframe *fr)
{
  char title[STRLEN];
  real prec;

  if (fr->bPrec)
    prec = fr->prec;
  else
    prec = 1000.0;

  switch (fio_getftp(fnum)) {
  case efTRJ:
  case efTRR:
    break;
  default:
    if (!fr->bX)
      fatal_error(0,"Need coordinates to write a %s trajectory",
		  ftp2ext(fio_getftp(fnum)));
    break;
  }

  switch (fio_getftp(fnum)) {
  case efXTC:
    write_xtc(fnum,fr->natoms,fr->step,fr->time,fr->box,fr->x,prec);
    break;
  case efTRJ:
  case efTRR:  
    fwrite_trn(fnum,fr->step,fr->time,fr->lambda,fr->box,fr->natoms,
	       fr->bX ? fr->x:NULL,fr->bV ? fr->v:NULL ,fr->bF ? fr->f:NULL);
    break;
  case efGRO:
  case efPDB:
  case efBRK:
  case efENT:
    if (!fr->bAtoms)
      fatal_error(0,"Can not write a %s file without atom names",
		  ftp2ext(fio_getftp(fnum)));
    sprintf(title,"frame t= %.3f",fr->time);
    if (fio_getftp(fnum) == efGRO)
      write_hconf_p(fio_getfp(fnum),title,fr->atoms,
		    prec2ndec(prec),fr->x,fr->bV ? fr->v : NULL,fr->box);
    else
      write_pdbfile(fio_getfp(fnum),title,
		    fr->atoms,fr->x,fr->box,0,fr->step);
    break;
  case efG87:
    write_gms(fio_getfp(fnum),fr->natoms,fr->x,fr->box);
    break;
  case efG96:
    write_g96_conf(fio_getfp(fnum),fr,-1,NULL); 
    break;
  default:
    fatal_error(0,"Sorry, write_trxframe can not write %s",
		ftp2ext(fio_getftp(fnum)));
    break;
  }

  return 0;
}

int write_trx(int fnum,int nind,atom_id *ind,t_atoms *atoms,
	      int step,real time,matrix box,rvec x[],rvec *v)
{
  t_trxframe fr;
  
  clear_trxframe(&fr,TRUE);
  fr.bStep = TRUE;
  fr.step = step;
  fr.bTime = TRUE;
  fr.time = time;
  fr.bAtoms = atoms!=NULL;
  fr.atoms = atoms;
  fr.bX = TRUE;
  fr.x = x;
  fr.bV = v!=NULL;
  fr.v = v;
  fr.bBox = TRUE;
  copy_mat(box,fr.box);
  
  return write_trxframe_indexed(fnum,&fr,nind,ind);
}

void close_trx(int status)
{
  fio_close(status);
}

int open_trx(char *outfile,char *filemode)
{
  if (filemode[0]!='w' && filemode[0]!='a')
    fatal_error(0,"Sorry, write_trx can only write");

  return fio_open(outfile,filemode);
}

static bool gmx_next_frame(int status,t_trxframe *fr)
{
  t_trnheader sh;
  bool bOK,bRet;
  
  bRet = FALSE;

  if (fread_trnheader(status,&sh,&bOK)) {
    fr->natoms=sh.natoms;
    fr->bStep=TRUE;
    fr->step=sh.step;
    fr->bTime=TRUE;
    fr->time=sh.t;
    fr->bLambda = TRUE;
    fr->lambda = sh.lambda;
    fr->bBox = sh.box_size>0;
    if (fr->flags & (TRX_READ_X | TRX_NEED_X)) {
      if (fr->x==NULL)
	snew(fr->x,sh.natoms);
      fr->bX = sh.x_size>0;
    }
    if (fr->flags & (TRX_READ_V | TRX_NEED_V)) {
      if (fr->v==NULL)
	snew(fr->v,sh.natoms);
      fr->bV = sh.v_size>0;
    }
    if (fr->flags & (TRX_READ_F | TRX_NEED_F)) {
      if (fr->f==NULL)
	snew(fr->f,sh.natoms);
      fr->bF = sh.f_size;
    }
    if (fread_htrn(status,&sh,fr->box,fr->x,fr->v,fr->f))
      bRet = TRUE;
    else
      fr->not_ok = DATA_NOT_OK;
  } else
    if (!bOK)
      fr->not_ok = HEADER_NOT_OK;

  return bRet;    
}

static void choose_ff(FILE *status)
{
  int i,m,c;

  printf("\n\n");
  printf("   Select File Format\n");
  printf("---------------------------\n");
  printf("1. XYZ File\n");
  printf("2. XYZ File with Box\n");
  printf("3. Gromos-87 Ascii Trajectory\n");
  printf("4. Gromos-87 Ascii Trajectory with Box\n");

  do {
    printf("\nChoice: ");
    fflush(stdout);
    scanf("%d",&i);
    i--;
  } while ((i < 0) || (i >= effNR));
  printf("\n");
  
  eFF = (eFileFormat) i;

  for(m=0; (m<DIM); m++) BOX[m]=0;
  
  bReadBox = (eFF == effG87Box) || (eFF == effXYZBox);
    
  switch (eFF) {
  case effXYZ:
  case effXYZBox:
    fscanf(status,"%d%lf%lf%lf%lf",&NATOMS,&BOX[XX],&BOX[YY],&BOX[ZZ],&DT);
    break;
  case effG87:
  case effG87Box:
    printf("GROMOS! OH DEAR...\n\n");
    printf("Number of atoms ? ");
    fflush(stdout);
    scanf("%d",&NATOMS);

    printf("Time between timeframes ? ");
    fflush(stdout);
    scanf("%lf",&DT);

    if (eFF == effG87) {
      printf("Box X Y Z ? ");
      fflush(stdout);
      scanf("%lf%lf%lf",&BOX[XX],&BOX[YY],&BOX[ZZ]);
    }
    do {
      c=fgetc(status);
      printf("%c",c);
    } while (c != '\n');
    printf("\n");
    fflush(stdout);
    break;
  default:
    printf("Hellow World\n");
  }
}

static bool do_read_xyz(FILE *status,int natoms,rvec x[],matrix box)
{
  int    i,m;
  double x0;

  for(i=0; (i<natoms); i++) {
    for(m=0; (m<DIM); m++) {
      if (fscanf(status,"%lf",&x0) != 1) {
	if (i || m)
	  fprintf(stderr,"error reading statusfile: x[%d][%d]\n",i,m);
	/* else eof! */
	return FALSE;
      }
      x[i][m]=x0;
    }
  }
  if (bReadBox) {
    for(m=0; (m<DIM); m++) {
      if (fscanf(status,"%lf",&x0) != 1) 
	return FALSE;
      box[m][m]=x0;
    }
  }
  return TRUE;
}

static bool xyz_next_x(FILE *status, real *t, int natoms, rvec x[], matrix box)
     /* Reads until a new x can be found (return TRUE)
      * or eof (return FALSE)
      */
{
  extern real tbegin,tend;
  real pt;
  
  pt=*t;
  while ((tbegin >= 0) && (*t < tbegin)) {
    if (!do_read_xyz(status,natoms,x,box))
      return FALSE;
    printcount(*t);
    *t+=DT;
    pt=*t;
  }
  if (((tend >= 0) && (*t < tend)) || (tend < 0.0)) {
    if (!do_read_xyz(status,natoms,x,box)) {
      printlast(*t);
      return FALSE;
    }
    printcount(*t);
    pt=*t;
    *t+=DT;
    return TRUE;
  }
  printlast(pt);
  return FALSE;
}

static int xyz_first_x(FILE *status, real *t, rvec **x, matrix box)
/* Reads status, mallocs x, and returns x and box
 * Returns natoms when succesful, FALSE otherwise
 */
{
  int    m;
  
  INITCOUNT;

  clear_mat(box);
  choose_ff(status);

  for(m=0; (m<DIM); m++)
    box[m][m]=BOX[m];

  snew(*x,NATOMS);
  *t=DT;
  if (!xyz_next_x(status,t,NATOMS,*x,box)) 
    return 0;
  *t=0.0;
  
  return NATOMS;
}

static bool pdb_next_x(FILE *status,t_trxframe *fr)
{
  t_atoms   atoms;
  int       model_nr, na;
  char      title[STRLEN],*time;
  double    dbl;
  
  atoms.nr = fr->natoms;
  atoms.atom=NULL;
  atoms.pdbinfo=NULL;
  /* the other pointers in atoms should not be accessed if these are NULL */
  model_nr=NOTSET;
  na=read_pdbfile(status, title, &model_nr, &atoms, fr->x, fr->box, TRUE);
  if (nframes_read()==0)
    fprintf(stderr," '%s', %d atoms\n",title, fr->natoms);
  fr->bPrec = TRUE;
  fr->prec = 10000;
  fr->bX = TRUE;
  fr->bBox = fr->box[XX][XX] == 0;
  
  if (model_nr!=NOTSET) {
    fr->bStep = TRUE;
    fr->step = model_nr;
  }
  time=strstr(title," t= ");
  if (time) {
    fr->bTime = TRUE;
    sscanf(time+4,"%lf",&dbl);
    fr->time=(real)dbl;
  } else {
    fr->bTime = FALSE;
    /* this is a bit dirty, but it will work: if no time is read from 
       comment line in pdb file, set time to current frame number */
    if (fr->bStep)
      fr->time=(real)fr->step;
    else
      fr->time=(real)nframes_read();
  }
  if (na==0) {
    return FALSE;
  } else { 
    if (na != fr->natoms)
      fatal_error(0,"Number of atoms in pdb frame %d is %d instead of %d",
		  nframes_read(),na,fr->natoms);
    return TRUE;
  }
}

static int pdb_first_x(FILE *status, t_trxframe *fr)
{
  INITCOUNT;
  
  fprintf(stderr,"Reading frames from pdb file");
  frewind(status);
  get_pdb_coordnum(status, &fr->natoms);
  if (fr->natoms==0)
    fatal_error(0,"\nNo coordinates in pdb file\n");
  frewind(status);
  snew(fr->x,fr->natoms);
  pdb_next_x(status, fr);

  return fr->natoms;
}

bool read_next_frame(int status,t_trxframe *fr)
{
  real pt;
  int  ct;
  bool bOK,bRet,bMissingData=FALSE,bSkip=FALSE;

  bRet = FALSE;
  pt=fr->time; 

  do {
    clear_trxframe(fr,FALSE);
    fr->tppf = fr->tpf;
    fr->tpf  = fr->time;
    
    switch (fio_getftp(status)) {
    case efTRJ:
    case efTRR:
    case efMTX:
      bRet = gmx_next_frame(status,fr);
      break;
    case efG96:
      read_g96_conf(fio_getfp(status),NULL,fr);
      bRet = (fr->natoms > 0);
      break;
    case efG87:
      bRet = xyz_next_x(fio_getfp(status),&fr->time,fr->natoms,fr->x,fr->box);
      fr->bTime = bRet;
      fr->bX    = bRet;
      fr->bBox  = bRet;
      break;
    case efXTC:
      bRet = read_next_xtc(status,fr->natoms,&fr->step,&fr->time,fr->box,
			   fr->x,&fr->prec,&bOK);
      fr->bPrec = bRet;
      fr->bStep = bRet;
      fr->bTime = bRet;
      fr->bX    = bRet;
      fr->bBox  = bRet;
      if (!bOK) {
	/* Actually the header could also be not ok,
	   but from bOK from read_next_xtc this can't be distinguished */
	fr->not_ok = DATA_NOT_OK;
      }
      break;
    case efPDB:
      bRet = pdb_next_x(fio_getfp(status),fr);
      break;
    case efGRO:
      bRet = gro_next_x_or_v(fio_getfp(status),fr);
      break;
    default:
      fatal_error(0,"DEATH HORROR in read_next_frame ftp=%s,status=%d",
		  ftp2ext(fio_getftp(status)),status);
    }
    
    if (bRet) {
      bMissingData = ((fr->flags & TRX_NEED_X && !fr->bX) ||
		      (fr->flags & TRX_NEED_V && !fr->bV) ||
		      (fr->flags & TRX_NEED_F && !fr->bF));
      bSkip = FALSE;
      if (!bMissingData) {
	ct=check_times2(fr->time,fr->t0,fr->tpf,fr->tppf);
	if (ct == 0 || (fr->flags & TRX_DONT_SKIP && ct<0)) {
	  printcount2(fr->time,fr->t0,fr->tpf,fr->tppf);
	} else if (ct > 0)
	  bRet = FALSE;
	else {
	  printcount2(fr->time,fr->t0,fr->tpf,fr->tppf);
	  bSkip = TRUE;
	}
      }
    }
    
  } while (bRet && (bMissingData || bSkip));
  
  if (!bRet) {
    printlast(pt);
    if (fr->not_ok)
      printincomp(fr);
  }
  
  return bRet;
}

int read_first_frame(int *status,char *fn,t_trxframe *fr,int flags)
{
  int  fp;
  bool bFirst,bOK;

  clear_trxframe(fr,TRUE);
  fr->flags = flags;

  bFirst = TRUE;
  INITCOUNT;
  
  fp = *status =fio_open(fn,"r");
  switch (fio_getftp(fp)) {
  case efTRJ:
  case efTRR:
  case efMTX:
    break;
  case efG96:
    /* Can not rewind a compressed file, so open it twice */
    read_g96_conf(fio_getfp(fp),fn,fr);
    fio_close(fp);
    clear_trxframe(fr,FALSE);
    if (flags & (TRX_READ_X | TRX_NEED_X))
      snew(fr->x,fr->natoms);
    if (flags & (TRX_READ_V | TRX_NEED_V))
      snew(fr->v,fr->natoms);
    fp = *status =fio_open(fn,"r");
    break;
  case efG87:
    fr->natoms=xyz_first_x(fio_getfp(fp),&fr->time,&fr->x,fr->box);
    if (fr->natoms) {
      fr->bTime = TRUE;
      fr->bX    = TRUE;
      fr->bBox  = TRUE;
      printcount(fr->time);
    }
    bFirst = FALSE;
    break;
  case efXTC:
    if (read_first_xtc(fp,&fr->natoms,&fr->step,&fr->time,fr->box,&fr->x,
		       &fr->prec,&bOK) == 0) {
      if (bOK) {
	fatal_error(0,"No XTC!\n");
      } else {
	fr->not_ok = DATA_NOT_OK;
      }
    }
    if (fr->not_ok) {
      fr->natoms = 0;
      printincomp(fr);
    } else {
      fr->bPrec = TRUE;
      fr->bStep = TRUE;
      fr->bTime = TRUE;
      fr->bX    = TRUE;
      fr->bBox  = TRUE;
      printcount(fr->time);
    }
    bFirst = FALSE;
    break;
  case efPDB:
    pdb_first_x(fio_getfp(fp),fr);
    if (fr->natoms)
      printcount(fr->time);
    bFirst = FALSE;
    break;
  case efGRO:
    if (gro_first_x_or_v(fio_getfp(fp),fr))
      printcount(fr->time);
    bFirst = FALSE;
    break;
  default:
    fatal_error(0,"Not supported in read_first_frame: %s",fn);
    break;
  }
  
  if (bFirst || 
      (!(fr->flags & TRX_DONT_SKIP) && (check_times(fr->time) < 0)))
    /* Read a frame when no frame was read or the first was skipped */
    if (!read_next_frame(*status,fr))
      return FALSE;
  fr->t0 = fr->time;
  
  return (fr->natoms > 0);
}

/***** C O O R D I N A T E   S T U F F *****/

int read_first_x(int *status,char *fn,
		 real *t,rvec **x,matrix box)
{
  t_trxframe fr;

  read_first_frame(status,fn,&fr,TRX_NEED_X);
  if (*status >= nxframe) {
    nxframe = *status+1;
    srenew(xframe,nxframe);
  }
  xframe[*status] = fr;
  *t = xframe[*status].time;
  *x = xframe[*status].x;
  copy_mat(xframe[*status].box,box);
  
  return xframe[*status].natoms;
}

bool read_next_x(int status,real *t, int natoms, rvec x[], matrix box)
{
  bool bRet;
  
  xframe[status].x = x;
  bRet = read_next_frame(status,&xframe[status]);
  *t = xframe[status].time;
  copy_mat(xframe[status].box,box);
  
  return bRet;
}

void close_trj(int status)
{
  fio_close(status);
}

void rewind_trj(int status)
{
  INITCOUNT;
  
  fio_rewind(status);
}

/***** V E L O C I T Y   S T U F F *****/

static void clear_v(t_trxframe *fr)
{
  int i;

  if (!fr->bV)
    for(i=0; i<fr->natoms; i++)
      clear_rvec(fr->v[i]);
}

int read_first_v(int *status,char *fn,real *t,rvec **v,matrix box)
{
  t_trxframe fr;

  read_first_frame(status,fn,&fr,TRX_NEED_V);
  *t = fr.time;
  clear_v(&fr);
  *v = fr.v;
  copy_mat(fr.box,box);
  
  return fr.natoms;
}

bool read_next_v(int status,real *t,int natoms,rvec v[],matrix box)
{
  t_trxframe fr;
  bool bRet;

  clear_trxframe(&fr,TRUE);
  fr.flags = TRX_NEED_V;
  fr.natoms = natoms;
  fr.time = *t;
  fr.v = v;
  bRet = read_next_frame(status,&fr);
  *t = fr.time;
  clear_v(&fr);
  copy_mat(fr.box,box);

  return bRet;
}
