/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_trxio_c = "$Id$";

#include <ctype.h>
#include "sysstuff.h"
#include "string2.h"
#include "smalloc.h"
#include "pbc.h"
#include "statusio.h"
#include "statutil.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "xtcio.h"
#include "pdbio.h"

/* defines for frame counter output */
static int frame=-666;
#define SKIP 10
#define INITCOUNT frame=-1;
#define PRINTCOUNT(l,t) fprintf(stderr,"\r%s frame %6d time %8.3f   ",l,frame,t)
#define CHECKCOUNT(l,t) if ( ((frame % SKIP)==0) || (frame < SKIP)) PRINTCOUNT(l,t)
#define PRINTSKIP(t) {frame++; CHECKCOUNT("Skipping",t);}
#define PRINTREAD(t) {frame++; CHECKCOUNT("Reading",t);}
#define PRINTLAST(t) { PRINTCOUNT("Last",t); fprintf(stderr,"\n"); }

/* Globals for gromos-87 input */
typedef enum { effXYZ, effXYZBox, effG87, effG87Box, effNR } eFileFormat;
static eFileFormat eFF;
static int         NATOMS;
static double      DT,BOX[3];
static bool        bReadBox;

static bool gmx_next_x(FILE *status,real *t,int natoms,rvec x[],matrix box)
{
  t_statheader sh;
  real lambda,pt;
  int  step,nre,ct;
  bool bB,bX;
  
  while (!eof(status)) {
    rd_header(status,&sh);
    bX=sh.x_size;
    bB=sh.box_size;
    pt=*t;
    rd_hstatus(status,&sh,&step,t,&lambda,NULL,
	       bB ? box : NULL,NULL,NULL,&sh.natoms,
	       bX ? x : NULL,NULL,NULL,
	       &nre,NULL,NULL);
    if (ct=check_times(*t)==0) {
      PRINTREAD(*t);
      if (bB)
	init_pbc(box,FALSE);  
      if (bX)
	return TRUE;
    } else if (ct > 0) {
      PRINTLAST(pt);
      return FALSE;
    } else {
      PRINTSKIP(*t);
    }
  }
  PRINTLAST(pt);
  return FALSE;    
}

static bool gmx_next_x_or_v(FILE *status,real *t,int natoms,
			    rvec x[],rvec v[],matrix box)
{
  t_statheader sh;
  real lambda,pt;
  int  i,d,step,nre,ct;
  bool bB,bX,bV;
    
  while (!eof(status)) {
    rd_header(status,&sh);
    bX=sh.x_size;
    bV=sh.v_size;
    bB=sh.box_size;
    pt=*t;
    rd_hstatus(status,&sh,&step,t,&lambda,NULL,
	       bB ? box : NULL,
	       NULL,NULL,&sh.natoms,
	       bX ? x : NULL,
	       bV ? v : NULL,
	       NULL,&nre,NULL,NULL);
	       
    if (ct=check_times(*t)==0) {
      PRINTREAD(*t);
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
      PRINTLAST(pt);
      return FALSE;
    } else {
      PRINTSKIP(*t);
    }
  }
  PRINTLAST(pt);
  return FALSE;    
}

static bool gmx_next_x_v(FILE *status,real *t,int natoms,
			 rvec x[],rvec v[],matrix box)
{
  t_statheader sh;
  real lambda,pt;
  int  step,nre,ct;
  bool bB,bX,bV;
    
  while (!eof(status)) {
    rd_header(status,&sh);
    bX=sh.x_size;
    bV=sh.v_size;
    bB=sh.box_size;
    pt=*t;
    rd_hstatus(status,&sh,&step,t,&lambda,NULL,
	       bB ? box : NULL,
	       NULL,NULL,&sh.natoms,
	       bX ? x : NULL,
	       bV ? v : NULL,
	       NULL,&nre,NULL,NULL);
	       
    if (ct=check_times(*t)==0) {
      PRINTREAD(*t);
      if (bB)
	init_pbc(box,FALSE);  
      if (bX && bV )
	return TRUE;
    } else if (ct > 0) {
      PRINTLAST(pt);
      return FALSE;
    } else {
      PRINTSKIP(*t);
    }
  }
  PRINTLAST(pt);
  return FALSE;    
}

static int gmx_first_x(FILE *status, real *t, rvec **x, matrix box)
{
  t_statheader sh;
  
  INITCOUNT;

  fprintf(stderr,"\nReading statusfile, version: %s\n",rd_header(status,&sh));

  snew(*x,sh.natoms);
  frewind(status);
  if (!gmx_next_x(status,t,sh.natoms,*x,box)) {
    fprintf(stderr,"No coordinates in trajectory\n");
    exit(1);
  }

  return sh.natoms;
}

static int gmx_first_x_v(FILE *status, real *t, rvec **x,rvec **v,matrix box)
{
  t_statheader sh;
  
  INITCOUNT;

  fprintf(stderr,"\nReading statusfile, version: %s\n",rd_header(status,&sh));
  snew(*x,sh.natoms);
  snew(*v,sh.natoms);
  frewind(status);
  if (!gmx_next_x_v(status,t,sh.natoms,*x,*v,box)) {
    fprintf(stderr,"No coordinates and velocities in trajectory\n");
    exit(1);
  }

  return sh.natoms;
}

static int gmx_first_x_or_v(FILE *status, real *t, 
			    rvec **x,rvec **v,matrix box)
{
  t_statheader sh;
  
  INITCOUNT;

  fprintf(stderr,"\nReading statusfile, version: %s\n",rd_header(status,&sh));
  snew(*x,sh.natoms);
  snew(*v,sh.natoms);
  frewind(status);
  if (!gmx_next_x_or_v(status,t,sh.natoms,*x,*v,box)) {
    fprintf(stderr,"No coordinates and velocities in trajectory\n");
    exit(1);
  }
  
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
      };
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
  
  while ((tbegin >= 0) && (*t < tbegin)) {
    if (!do_read_xyz(status,natoms,x,box))
      return FALSE;
    PRINTSKIP(*t);
    *t+=DT;
  }
  if (((tend >= 0) && (*t < tend)) || (tend < 0.0)) {
    if (!do_read_xyz(status,natoms,x,box)) {
      PRINTLAST(*t);
      return FALSE;
    }
    PRINTREAD(*t);
    pt=*t;
    *t+=DT;
    return TRUE;
  }
  PRINTLAST(pt);
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
  t_pdbatom *pdb;
  int i;

  *t=0.0;  
  natoms = read_pdbatoms(status, &pdb, box, FALSE);
  fprintf(stderr,"\rRead frame %6d",frame++);
  
  for(i=0; (i<natoms); i++)
    copy_rvec(pdb[i].x, x[i]);
    
  sfree(pdb);

  if (natoms==0) {
    return FALSE;
  }
  else { 
    return TRUE;
  }
}

static int pdb_first_x(FILE *status, real *t, rvec **x, matrix box)
{
  t_pdbatom *pdb;
  int i, j, natoms;
  
  INITCOUNT;
  
  *t=0.0;
  fprintf(stderr,"Reading frames from pdb file.\n");
  frewind(status);
  natoms = read_pdbatoms(status, &pdb, box, FALSE);
  fprintf(stderr,"No of atoms: %d.\n", natoms);

  if (box[XX][XX] == 0)
    box[XX][XX]=box[YY][YY]=box[ZZ][ZZ]=10;
    
  if (natoms==0) {
    fprintf(stderr,"No coordinates in pdb file\n");
    exit(1);
  }
  else {
    snew(*x,natoms);
    for(i=0; (i<natoms); i++)
      copy_rvec(pdb[i].x, (*x)[i] );
  }
  sfree(pdb);

  return natoms;
}

/***** C O O R D I N A T E   S T U F F *****/

typedef struct {
  int  ftp;
  char *fn;
  FILE *fp;
  XDR  xd;
} t_trjf;

static t_trjf *trjf=NULL;
static int    ntrjf=0;

int read_first_x(int *status,char *fn,
		 real *t,rvec **x,matrix box)
{
  char buf[256];
  int  natoms,step,ftp;
  real prec;

  INITCOUNT;
  
  ftp=fn2ftp(fn);
  *status=ntrjf;
  srenew(trjf,ntrjf+1);
  trjf[ntrjf].ftp=ftp;
  trjf[ntrjf].fn=strdup(fn);
  ntrjf++;

  switch (ftp) {
  case efTRJ:
    trjf[ntrjf-1].fp=ffopen(fn,"r");
    natoms=gmx_first_x(trjf[ntrjf-1].fp,t,x,box);
    break;
  case efMTX:
    trjf[ntrjf-1].fp=ffopen(fn,"r");
    natoms=gmx_first_x(trjf[ntrjf-1].fp,t,x,box);
    break;
  case efG87:
    trjf[ntrjf-1].fp=ffopen(fn,"r");
    natoms=xyz_first_x(trjf[ntrjf-1].fp,t,x,box);
    break;
  case efXTC:
    if (read_first_xtc(&(trjf[ntrjf-1].xd),fn,
		       &natoms,&step,t,box,x,&prec) == 0)
      fatal_error(0,"No XTC!\n");
    if (check_times(*t) < 0)
      if (!read_next_x(*status,t,natoms,*x,box))
	return 0;
    PRINTREAD(*t);
    break;
  case efPDB:
    trjf[ntrjf-1].fp=ffopen(fn,"r");
    natoms=pdb_first_x(trjf[ntrjf-1].fp,t,x,box);
    break;
  case efGRO:
    trjf[ntrjf-1].fp=ffopen(fn,"r");
    natoms=gro_first_x(trjf[ntrjf-1].fp,t,x,box);
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
  
  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_x)",status);
    
  switch (trjf[status].ftp) {
  case efTRJ:
    return gmx_next_x(trjf[status].fp,t,natoms,x,box);
  case efMTX:
    return gmx_next_x(trjf[status].fp,t,natoms,x,box);
  case efG87:
    return xyz_next_x(trjf[status].fp,t,natoms,x,box);
  case efXTC:
    pt=*t;
    while (read_next_xtc(&(trjf[status].xd),&natoms,&step,t,box,x,&prec)) {
      if ((ct=check_times(*t)) == 0) {
	PRINTREAD(*t);
	init_pbc(box,FALSE);  
	return TRUE;
      }
      else if (ct > 0) {
	PRINTLAST(pt);
	return FALSE;
      } else {
	PRINTSKIP(*t);
      }
    }
    PRINTLAST(pt);
    return FALSE;
  case efPDB:
    return pdb_next_x(trjf[status].fp,t,natoms,x,box);
  case efGRO:
    pt=*t;
    if (gro_next_x(trjf[status].fp,t,natoms,x,box)) {
      PRINTREAD(*t);
      return TRUE;
    } else {
      PRINTLAST(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x ftp=%d,status=%d",
		trjf[status].ftp,status);
  }
  return FALSE;
}

void close_trj(int status)
{
  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_x)",status);
    
  switch (trjf[status].ftp) {
  case efTRJ:
  case efMTX:
  case efG87:
  case efPDB:
  case efGRO:
    fclose(trjf[status].fp);
    break;
  case efXTC:
    xdrclose(&trjf[status].xd);
    break;
  default:
    fatal_error(0,"\nDeath HORROR in close_trj ftp=%d",trjf[status].ftp);
  }
  trjf[status].ftp=-1;
}

void rewind_trj(int status)
{
  char buf[256];
  
  INITCOUNT;
  
  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_x)",status);
  
  switch (trjf[status].ftp) {
  case efTRJ:
    frewind(trjf[status].fp);
    break;
  case efG87:
    frewind(trjf[status].fp);
    fgets2(buf,254,trjf[status].fp);
    break;
  case efXTC:
    xdrclose(&trjf[status].xd);
    xdropen(&trjf[status].xd,trjf[status].fn,"r");
    break;
  case efPDB:
    frewind(trjf[status].fp);
    break;
  default:
    fatal_error(0,"\nDeath HORROR in rewind_trj ftp=%d",trjf[status].ftp);
  }
}

/***** V E L O C I T Y   S T U F F *****/

int read_first_v(int *status,char *fn,real *t,rvec **v,matrix box)
{
  t_statheader sh;
  real lambda;
  int  step,nre,natoms;

  INITCOUNT;
  
  *status=ntrjf;
  srenew(trjf,ntrjf+1);
  trjf[ntrjf].ftp=fn2ftp(fn);
  trjf[ntrjf].fn=strdup(fn);
  trjf[ntrjf].fp=ffopen(fn,"r");
  ntrjf++;
  
  switch (trjf[ntrjf-1].ftp) {
  case efTRJ:
    fprintf(stderr,"Reading trj file, version: %s\n",
	    rd_header(trjf[ntrjf-1].fp,&sh));
    snew(*v,sh.natoms);
    rd_hstatus(trjf[ntrjf-1].fp,
	       &sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
	       &sh.natoms,NULL,*v,NULL,
	       &nre,NULL,NULL);
    return sh.natoms;
  case efGRO:
    return gro_first_v(trjf[ntrjf-1].fp,t,v,box);
  default:
    fatal_error(0,"Not supported in read_first_v: %s",trjf[ntrjf-1].fn);
  }
    
  return 0;
}

bool read_next_v(int status,real *t,int natoms,rvec v[],matrix box)
{
  t_statheader sh;
  real lambda,pt;
  int  step,nre;
  bool bV;

  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_*)",status);
    
  switch (trjf[status].ftp) {
  case efTRJ:
    while (!eof(trjf[status].fp)) {
      rd_header(trjf[status].fp,&sh);
      bV=sh.v_size;
      rd_hstatus(trjf[status].fp,
		 &sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
		 &sh.natoms,NULL,bV ? v : NULL,NULL,
		 &nre,NULL,NULL);
      if ((check_times(*t)==0) && (bV))
	return TRUE;
      if (check_times(*t) > 0)
	return FALSE;
    }
    PRINTREAD(*t);
    break;
  case efGRO: 
    pt=*t;
    if (gro_next_v(trjf[status].fp,t,natoms,v,box)) {
      PRINTREAD(*t);
      return TRUE;
    } else {
      PRINTLAST(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"DEATH HORROR in read_next_v: ftp=%d,status=%d",
		trjf[status].ftp,status);
  }
  
  return FALSE;
}


/* coordinates and velocities */

int read_first_x_v(int *status,char *fn,
		   real *t,rvec **x,rvec **v,matrix box)
{
  char buf[256];
  int  natoms,step,ftp;
  real prec;
  
  INITCOUNT;

  *status=ntrjf;
  srenew(trjf,ntrjf+1);
  trjf[ntrjf].ftp=fn2ftp(fn);
  trjf[ntrjf].fn=strdup(fn);
  trjf[ntrjf].fp=ffopen(fn,"r");
  ntrjf++;

  switch (trjf[ntrjf-1].ftp) {
  case efTRJ:
    natoms=gmx_first_x_v(trjf[ntrjf-1].fp,t,x,v,box);
    break;
  case efGRO:
    natoms=gro_first_x_v(trjf[ntrjf-1].fp,t,x,v,box);
    break;
  default:
    fatal_error(0,"Not supported in read_first_x_v: %s",trjf[ntrjf-1].fn);
  }

  return natoms;
}

bool read_next_x_v(int status,real *t, int natoms, 
		   rvec x[],rvec v[],matrix box)
{
  int step,ct;
  real prec,pt;
  
  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_*)",status);
    
  switch (trjf[status].ftp) {
  case efTRJ:
    return gmx_next_x_v(trjf[status].fp,t,natoms,x,v,box);
  case efGRO: 
    pt=*t;
    if (gro_next_x_v(trjf[status].fp,t,natoms,x,v,box)) {
      PRINTREAD(*t);
      return TRUE;
    } else {
      PRINTLAST(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x_v ftp=%d,status=%d",
		trjf[status].ftp,status);
  }
  return FALSE;
}

int read_first_x_or_v(int *status,char *fn,
		   real *t,rvec **x,rvec **v,matrix box)
{
  char buf[256];
  int  natoms,step,ftp;
  real prec;
  
  INITCOUNT;

  *status=ntrjf;
  srenew(trjf,ntrjf+1);
  trjf[ntrjf].ftp=fn2ftp(fn);
  trjf[ntrjf].fn=strdup(fn);
  trjf[ntrjf].fp=ffopen(fn,"r");
  ntrjf++;

  switch (trjf[ntrjf-1].ftp) {
  case efTRJ:
    natoms=gmx_first_x_or_v(trjf[ntrjf-1].fp,t,x,v,box);
    break;
  case efGRO:
    natoms=gro_first_x_or_v(trjf[ntrjf-1].fp,t,x,v,box);
    break;
  default:
    fatal_error(0,"Not supported in read_first_x_or_v: %s",trjf[ntrjf-1].fn);
  }

  return natoms;
}

bool read_next_x_or_v(int status,real *t, int natoms, 
		      rvec x[],rvec v[],matrix box)
{
  int step,ct;
  real prec,pt;
  
  if (status >= ntrjf)
    fatal_error(0,"File %d not opened yet... (call read_first_*)",status);
    
  switch (trjf[status].ftp) {
  case efTRJ:
    return gmx_next_x_or_v(trjf[status].fp,t,natoms,x,v,box);
  case efGRO:
    pt=*t;
    if (gro_next_x_or_v(trjf[status].fp,t,natoms,x,v,box)) {
      PRINTREAD(*t);
      return TRUE;
    } else {
      PRINTLAST(pt);
      return FALSE;
    }
  default:
    fatal_error(0,"\nDEATH HORROR ERROR in read_next_x_or_v ftp=%d,status=%d",
		trjf[status].ftp,status);
  }
  return FALSE;
}
