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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_trxio_c = "$Id$";

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

/* defines for frame counter output */
static int frame=-666;
#define SKIP 10
#define INITCOUNT frame=-1

static void printcount(char *l,real t)
{
  fprintf(stderr,"\r%s frame %6d time %8.3f   ",l,frame,t);
}

#define CHECKCOUNT(l,t) if ( ((frame % SKIP)==0) || (frame < SKIP)) printcount(l,t)

static void printskip(real t)
{
  frame++;
  CHECKCOUNT("Skipping",t);
}

static void printread(real t)
{
  frame++;
  CHECKCOUNT("Reading",t);
}

static void printlast(real t)
{
  printcount("Last",t);
  fprintf(stderr,"\n");
}

static void printincomp(real t, real pt)
{
  fprintf(stderr,"WARNING: Incomplete frame: nr %d time %g\n",frame+1,t);
}

static void printincomph(real t)
{
  fprintf(stderr,"WARNING: Incomplete frame header: nr %d time %g\n",
	  frame+1,t);
}

/* Globals for gromos-87 input */
typedef enum { effXYZ, effXYZBox, effG87, effG87Box, effNR } eFileFormat;
static eFileFormat eFF;
static int         NATOMS;
static double      DT,BOX[3];
static bool        bReadBox;

typedef struct {
  int  fnum;
  FILE *fp;
  int  ftp;
} t_trx_out;

static int nfile=0;
static t_trx_out *trx_out=NULL;

int write_trx(int fnum,int nind,atom_id *ind,t_atoms *atoms,
	      int step,real time,matrix box,rvec x[],rvec *v)
{
  char title[STRLEN];
  rvec *xout,*vout;
  int i;

  switch (trx_out[fnum].ftp) {
  case efTRN: 
    if (v) {
      snew(vout,nind);
      for(i=0; i<nind; i++) 
	copy_rvec(v[ind[i]],vout[i]);
    }
  case efXTC:
  case efG87:
    snew(xout,nind);
    for(i=0; i<nind; i++) 
      copy_rvec(x[ind[i]],xout[i]);
    break;
  default:
    break;
  }

  switch (trx_out[fnum].ftp) {
  case efXTC: 
    write_xtc(trx_out[fnum].fnum,nind,step,time,box,xout,1000);
    break;
  case efTRN:
    fwrite_trn(trx_out[fnum].fnum,frame,time,step,box,nind,xout,vout,NULL);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    sprintf(title,"frame t= %.3f",time);
    write_pdbfile_indexed(trx_out[fnum].fp,title,atoms,x,box,0,TRUE,nind,ind);
    break;
  case efGRO:
    sprintf(title,"frame t= %.3f",time);
    write_hconf_indexed(trx_out[fnum].fp,title,atoms,nind,ind,x,v,box);
    break;
  case efG87:
    write_gms(trx_out[fnum].fp,nind,xout,box);
    break;
  default:
    fatal_error(0,"Sorry, write_trx_x can not write %s",
		ftp2ext(trx_out[fnum].ftp));
    break;
  }

  switch (trx_out[fnum].ftp) {
  case efTRN:
    if (v) sfree(vout);
  case efXTC:
  case efG87:
    sfree(xout);
    break;
  default:
    break;
  }
  
  return 0;
}

int close_trx(int fnum)
{
  switch (trx_out[fnum].ftp) {
  case efXTC: 
    close_xtc(trx_out[fnum].fnum);
    break;
  case efTRN:
    close_trn(trx_out[fnum].fnum);
    break;
  case efPDB:
  case efBRK:
  case efENT:
  case efGRO:
  case efG87:  
    fclose(trx_out[fnum].fp);
    break;
  default:
    fatal_error(0,"Sorry, write_trx_x can not close file %d",fnum);
    break;
  }

  return 0;
}

int open_trx(char *outfile,char *filemode)
{
  if (filemode[0] != 'w')
    fatal_error(0,"Sorry, write_trx_x can only write");
  
  srenew(trx_out,nfile+1);
  
  trx_out[nfile].ftp = fn2ftp(outfile);
  
  switch (trx_out[nfile].ftp) {
  case efXTC: 
    trx_out[nfile].fnum = open_xtc(outfile,filemode);
    break;
  case efTRN:
    trx_out[nfile].fnum = open_trn(outfile,filemode);
    break;
  case efPDB:
  case efBRK:
  case efENT:
  case efGRO:
  case efG87:  
    trx_out[nfile].fp = ffopen(outfile,filemode);
    break;
  default:
    fatal_error(0,"Sorry write_trx_x can not write %s",outfile);
    break;
  }
  
  nfile++;

  return nfile-1;
}

static bool gmx_next_x(int status,real *t,int natoms,rvec x[],matrix box)
{
  t_trnheader sh;
  real pt;
  int  ct;
  bool bB,bX,bOK;
  
  pt=*t;
  while (fread_trnheader(status,&sh,&bOK)) {
    bX = sh.x_size;
    bB = sh.box_size;
    pt = *t;
    *t = sh.t;
    if (!fread_htrn(status,&sh,
		    bB ? box : NULL,
		    bX ? x : NULL,
		    NULL,
		    NULL)) {
      printlast(pt);
      printincomp(*t,pt);
      return FALSE;
    }
    if ((ct=check_times(*t))==0) {
      printread(*t);
      if (bB)
	init_pbc(box,FALSE);  
      if (bX)
	return TRUE;
    } else if (ct > 0) {
      printlast(pt);
      return FALSE;
    } else {
      printskip(*t);
    }
  }
  printlast(pt);
  if (!bOK) printincomph(sh.t);

  return FALSE;    
}

static bool gmx_next_x_or_v(int status,real *t,int natoms,
			    rvec x[],rvec v[],matrix box)
{
  t_trnheader sh;
  real pt;
  int  i,d,ct;
  bool bB,bX,bV,bOK;
    
  pt=*t;
  while (fread_trnheader(status,&sh,&bOK)) {
    bX=sh.x_size;
    bV=sh.v_size;
    bB=sh.box_size;
    pt=*t;
    *t=sh.t;
    if (!fread_htrn(status,&sh,
		    bB ? box : NULL,
		    bX ? x : NULL,
		    bV ? v : NULL,
		    NULL)) {
      printlast(pt);
      printincomp(*t,pt);
      return FALSE;
    }
    if ((ct=check_times(*t))==0) {
      printread(*t);
      if (bB)
	init_pbc(box,FALSE);  
      if (bX || bV ) {
	if (!bX)
	  for (i=0; (i<sh.natoms); i++)
	    for (d=0; (d<DIM); d++)
	      x[i][d]=0;
	if (!bV)
	  for (i=0; (i<sh.natoms); i++)
	    for (d=0; (d<DIM); d++)
	      v[i][d]=0;
	return TRUE;
      }
    } else if (ct > 0) {
      printlast(pt);
      return FALSE;
    } else {
      printskip(*t);
    }
  }
  printlast(pt);
  if (!bOK)
    printincomph(sh.t);
      
  return FALSE;    
}
  
static bool gmx_next_x_v(int status,real *t,int natoms,
			 rvec x[],rvec v[],matrix box)
{
  t_trnheader sh;
  real pt;
  int  ct;
  bool bB,bX,bV,bOK;
    
  pt=*t;
  while (fread_trnheader(status,&sh,&bOK)) {
    bX=sh.x_size;
    bV=sh.v_size;
    bB=sh.box_size;
    pt=*t;
    *t=sh.t;
    if (!fread_htrn(status,&sh,
		    bB ? box : NULL,
		    bX ? x : NULL,
		    bV ? v : NULL,
		    NULL)) {
      printlast(pt);
      printincomp(*t,pt);
      return FALSE;
    }
    if ((ct=check_times(*t))==0) {
      printread(*t);
      if (bB)
	init_pbc(box,FALSE);  
      if (bX && bV )
	return TRUE;
    } else if (ct > 0) {
      printlast(pt);
      return FALSE;
    } else {
      printskip(*t);
    }
  }
  printlast(pt);
  if (!bOK) printincomph(sh.t);

  return FALSE;    
}

static int gmx_first_x(int status, real *t, rvec **x, matrix box)
{
  t_trnheader sh;
  bool bOK;
  
  INITCOUNT;

  if (fread_trnheader(status,&sh,&bOK)) {
    snew(*x,sh.natoms);
    rewind_trj(status);
    if (!gmx_next_x(status,t,sh.natoms,*x,box)) {
      fprintf(stderr,"No coordinates in trajectory\n");
      exit(1);
    }
  }
  if (!bOK) printincomph(sh.t);

  return sh.natoms;
}

static int gmx_first_x_v(int status, real *t, rvec **x,rvec **v,matrix box)
{
  t_trnheader sh;
  bool bOK;

  INITCOUNT;
  
  if (fread_trnheader(status,&sh,&bOK)) {
    snew(*x,sh.natoms);
    snew(*v,sh.natoms);
    rewind_trj(status);
    if (!gmx_next_x_v(status,t,sh.natoms,*x,*v,box)) {
      fprintf(stderr,"No coordinates and velocities in trajectory\n");
      exit(1);
    }
  }
  if (!bOK) printincomph(sh.t);

  return sh.natoms;
}

static int gmx_first_x_or_v(int status, real *t,rvec **x,rvec **v,matrix box)
{
  t_trnheader sh;
  bool bOK;
  
  INITCOUNT;
  
  if (fread_trnheader(status,&sh,&bOK)) {
    snew(*x,sh.natoms);
    snew(*v,sh.natoms);
    rewind_trj(status);
    if (!gmx_next_x_or_v(status,t,sh.natoms,*x,*v,box)) {
      fprintf(stderr,"No coordinates and velocities in trajectory\n");
      exit(1);
    }
  }
  if (!bOK) printincomph(sh.t);
  
  return sh.natoms;
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
    printskip(*t);
    *t+=DT;
    pt=*t;
  }
  if (((tend >= 0) && (*t < tend)) || (tend < 0.0)) {
    if (!do_read_xyz(status,natoms,x,box)) {
      printlast(*t);
      return FALSE;
    }
    printread(*t);
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

  init_pbc(box,FALSE);
  
  return NATOMS;
}

static bool pdb_next_x(FILE *status,real *t,int natoms,rvec x[],matrix box)
{
  t_atoms   atoms;
  int       i,na;
  char      title[STRLEN],*time;
  double    dbl;

  init_t_atoms(&atoms,natoms,FALSE);
  na=read_pdbfile(status, title, &atoms, x, box, TRUE);
  if (box[XX][XX] == 0)
    box[XX][XX]=box[YY][YY]=box[ZZ][ZZ]=10;
  time=strstr(title," t= ");
  if (time) {
    sscanf(time+4,"%lf",&dbl);
    *t=(real)dbl;
  } else
    *t=0.0;
  /* free_t_atoms(&atoms); */
  if (na==0) {
    return FALSE;
  } else { 
    if (na != natoms)
      fatal_error(0,"Number of atoms in pdb frame %d is %d instead of %d",
		  frame,na,natoms);
    return TRUE;
  }
}

static int pdb_first_x(FILE *status, real *t, rvec **x, matrix box)
{
  int   i, natoms;
  char  title[STRLEN];
  
  INITCOUNT;
  
  *t=0.0;
  fprintf(stderr,"Reading frames from pdb file.\n");
  frewind(status);
  get_pdb_coordnum(status, &natoms);
  fprintf(stderr,"Nr of atoms: %d.\n", natoms);
  if (natoms==0)
    fatal_error(0,"No coordinates in pdb file\n");
  frewind(status);
  snew(*x,natoms);
  pdb_next_x(status, t, natoms, *x, box);

  return natoms;
}

/***** C O O R D I N A T E   S T U F F *****/

int read_first_x(int *status,char *fn,
		 real *t,rvec **x,matrix box)
{
  int  fp;
  int  natoms,step;
  real prec;
  bool bOK;

  INITCOUNT;
  
  fp = *status =fio_open(fn,"r");

  switch (fio_getftp(fp)) {
  case efTRJ:
  case efTRR:
  case efMTX:
    natoms=gmx_first_x(fp,t,x,box);
    break;
  case efG87:
    natoms=xyz_first_x(fio_getfp(fp),t,x,box);
    break;
  case efXTC:
    if (read_first_xtc(fp,&natoms,&step,t,box,x,&prec,&bOK) == 0)
      fatal_error(0,"No XTC!\n");
    if (check_times(*t) < 0)
      if (!read_next_x(*status,t,natoms,*x,box))
	return 0;
    printread(*t);
    break;
  case efPDB:
    natoms=pdb_first_x(fio_getfp(fp),t,x,box);
    if (natoms)
      printread(*t);
    break;
  case efGRO:
    natoms=gro_first_x(fio_getfp(fp),t,x,box);
    break;
  default:
    fatal_error(0,"Not supported in read_first_x: %s",fn);
  }
  
  return natoms;
}

bool read_next_x(int status,real *t, int natoms, rvec x[], matrix box)
{
  int step,ct;
  real prec,pt;
  bool bOK;
  
  switch (fio_getftp(status)) {
  case efTRJ:
  case efTRR:
  case efMTX:
    return gmx_next_x(status,t,natoms,x,box);
  case efG87:
    return xyz_next_x(fio_getfp(status),t,natoms,x,box);
  case efXTC:
    pt=*t;
    while (read_next_xtc(status,&natoms,&step,t,box,x,&prec,&bOK)) {
      if ((ct=check_times(*t)) == 0) {
	printread(*t);
	init_pbc(box,FALSE);  
	return TRUE;
      }
      else if (ct > 0) {
	printlast(pt);
	return FALSE;
      } else {
	printskip(*t);
      }
    }
    printlast(pt);
    if (!bOK)
      printincomp(*t,pt);
    return FALSE;
  case efPDB:
    pt=*t;
    if (pdb_next_x(fio_getfp(status),t,natoms,x,box)) {
      printread(*t);
      return TRUE;
    } else {
      printlast(pt);
      return FALSE;
    }
  case efGRO:
    pt=*t;
    if (gro_next_x(fio_getfp(status),t,natoms,x,box)) {
      printread(*t);
      return TRUE;
    } else {
      printlast(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x ftp=%s,status=%d",
		ftp2ext(fio_getftp(status)),status);
  }
  return FALSE;
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

int read_first_v(int *status,char *fn,real *t,rvec **v,matrix box)
{
  t_trnheader sh;
  int  fp;
  bool bOK;

  INITCOUNT;
  
  fp = *status = fio_open(fn,"r");
  
  switch (fio_getftp(fp)) {
  case efTRJ:
  case efTRR:
    if (fread_trnheader(fp,&sh,&bOK)) {
      snew(*v,sh.natoms);
      *t = sh.t;
      if (!fread_htrn(fp,&sh,NULL,NULL,*v,NULL)) {
	printincomp(*t,-1);
	return FALSE;
      }
      printread(*t);
      return sh.natoms;
    }
    if (!bOK) printincomph(sh.t);
    return -1;
  case efGRO:
    return gro_first_v(fio_getfp(fp),t,v,box);
  default:
    fatal_error(0,"Not supported in read_first_v: %s",fn);
  }
    
  return 0;
}

bool read_next_v(int status,real *t,int natoms,rvec v[],matrix box)
{
  t_trnheader sh;
  real pt;
  bool bV,bOK;
  
  pt=*t;
  switch (fio_getftp(status)) {
  case efTRJ:
  case efTRR:
    while (fread_trnheader(status,&sh,&bOK)) {
      bV=sh.v_size;
      *t = sh.t;
      if (!fread_htrn(status,&sh,NULL,NULL,bV ? v : NULL,NULL)) {
	printlast(pt);
	printincomp(*t,pt);
	return FALSE;
      } 
      if ((check_times(*t)==0) && (bV)) {
	printread(*t);
	return TRUE;
      }
      if (check_times(*t) > 0) {
	printlast(pt);
	return FALSE;
      }
    }
    printlast(pt);
    if (!bOK) printincomph(sh.t);
    break;
  case efGRO: 
    if (gro_next_v(fio_getfp(status),t,natoms,v,box)) {
      printread(*t);
      return TRUE;
    } else {
      printlast(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"DEATH HORROR in read_next_v: ftp=%s,status=%d",
		ftp2ext(fio_getftp(status)),status);
  }
  
  return FALSE;
}


/* coordinates and velocities */

int read_first_x_v(int *status,char *fn,
		   real *t,rvec **x,rvec **v,matrix box)
{
  int  fp,natoms;
  
  INITCOUNT;

  fp = *status = fio_open(fn,"r");

  switch (fio_getftp(fp)) {
  case efTRJ:
  case efTRR:
    natoms=gmx_first_x_v(fp,t,x,v,box);
    break;
  case efGRO:
    natoms=gro_first_x_v(fio_getfp(fp),t,x,v,box);
    break;
  default:
    fatal_error(0,"Not supported in read_first_x_v: %s",fn);
  }

  return natoms;
}

bool read_next_x_v(int status,real *t, int natoms, 
		   rvec x[],rvec v[],matrix box)
{
  real pt;
  
  switch (fio_getftp(status)) {
  case efTRJ:
  case efTRR:
    return gmx_next_x_v(status,t,natoms,x,v,box);
  case efGRO: 
    pt=*t;
    if (gro_next_x_v(fio_getfp(status),t,natoms,x,v,box)) {
      printread(*t);
      return TRUE;
    } else {
      printlast(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x_v ftp=%s,status=%d",
		ftp2ext(fio_getftp(status)),status);
  }
  return FALSE;
}

int read_first_x_or_v(int *status,char *fn,
		   real *t,rvec **x,rvec **v,matrix box)
{
  int  fp,natoms;
  
  INITCOUNT;

  fp = *status = fio_open(fn,"r");

  switch (fio_getftp(fp)) {
  case efTRJ:
  case efTRR:
    natoms=gmx_first_x_or_v(fp,t,x,v,box);
    break;
  case efGRO:
    natoms=gro_first_x_or_v(fio_getfp(fp),t,x,v,box);
    break;
  default:
    fatal_error(0,"Not supported in read_first_x_or_v: %s",fn);
  }

  return natoms;
}

bool read_next_x_or_v(int status,real *t, int natoms, 
		      rvec x[],rvec v[],matrix box)
{
  real pt;

  switch (fio_getftp(status)) {
  case efTRJ:
  case efTRR:
    return gmx_next_x_or_v(status,t,natoms,x,v,box);
  case efGRO:
    pt=*t;
    if (gro_next_x_or_v(fio_getfp(status),t,natoms,x,v,box)) {
      printread(*t);
      return TRUE;
    } else {
      printlast(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x_or_v ftp=%s,status=%d",
		ftp2ext(fio_getftp(status)),status);
  }
  return FALSE;
}
