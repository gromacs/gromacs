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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_xtcio_c = "$Id$";

#include <string.h>
#include "typedefs.h"
#include "xdrf.h"
#include "gmxfio.h"
#include "xtcio.h"
#include "smalloc.h"
#include "vec.h"
#include "futil.h"
#include "fatal.h"

#define XTC_MAGIC 1995

int open_xtc(char *fn,char *mode)
{
  return fio_open(fn,mode);
}

void close_xtc(int fp)
{
  fio_close(fp);
}

static void check_xtc_magic(int magic)
{
  if (magic != XTC_MAGIC) 
    fatal_error(0,"Magic Number Error in XTC file (read %d, should be %d)",
		magic,XTC_MAGIC);
}

int xtc_check(char *str,bool bResult,char *file,int line)
{
  if (!bResult) {
    fprintf(stderr,"XTC error: read/write of %s failed, "
	    "source file %s, line %d\n",str,file,line);
    return 0;
  }
  return 1;
}

void xtc_check_fat_err(char *str,bool bResult,char *file,int line)
{
  if (!bResult) {
    fatal_error(0,"XTC error: read/write of %s failed, "
		"source file %s, line %d\n",str,file,line);
  }
}

static int xtc_header(XDR *xd,int *magic,int *natoms,int *step,real *time)
{
  int result;

  if (xdr_int(xd,magic) == 0)
    return 0;
  result=XTC_CHECK("natoms", xdr_int(xd,natoms));  /* number of atoms */
  if (result)
    result=XTC_CHECK("step",   xdr_int(xd,step));    /* frame number    */
  if (result)
    result=XTC_CHECK("time",   xdr_real(xd,time));   /* time            */
    
  return result;
}

static int xtc_coord(XDR *xd,
		     int *natoms,
		     matrix box,rvec *x,real *prec)
{
  int i,j,result;
  
  /* box */
  result=1;
  for(i=0; ((i<DIM) && result); i++)
    for(j=0; ((j<DIM) && result); j++)
      result=XTC_CHECK("box",xdr_real(xd,&(box[i][j])));
  
  if (result)
    /* coordinates     */
    result=XTC_CHECK("x",xdr3drcoord(xd,x[0],natoms,prec)); 
  
  return result;
}

static int xtc_io(XDR *xd,int *magic,
		  int *natoms,int *step,real *time,
		  matrix box,rvec *x,real *prec)
{
  if (!xtc_header(xd,magic,natoms,step,time))
    return 0;
  return xtc_coord(xd,natoms,box,x,prec);
}

int write_xtc(int fp,
	      int natoms,int step,real time,
	      matrix box,rvec *x,real prec)
{
  int magic_number = XTC_MAGIC;
  XDR *xd;
  
  xd = fio_getxdr(fp);
  /* write magic number and xtc identidier */
  if (!xtc_header(xd,&magic_number,&natoms,&step,&time))
    return 0;
    
  /* write data */
  return xtc_coord(xd,&natoms,box,x,&prec);
}

int read_first_xtc(int fp,int *natoms,int *step,real *time,
		   matrix box,rvec **x,real *prec)
{
  int magic;
  XDR *xd;
  
  xd = fio_getxdr(fp);
  
  /* read header and malloc x */
  if ( !xtc_header(xd,&magic,natoms,step,time))
    return 0;
    
  /* Check magic number */
  check_xtc_magic(magic);
  
  snew(*x,*natoms);
  
  return xtc_coord(xd,natoms,box,*x,prec);
}

int read_next_xtc(int fp,
		  int *natoms,int *step,real *time,
		  matrix box,rvec *x,real *prec)
{
  int magic;
  XDR *xd;
  
  xd = fio_getxdr(fp);
  
  /* read header */
  if ( !xtc_header(xd,&magic,natoms,step,time))
    return 0;
    
  /* Check magic number */
  check_xtc_magic(magic);
  
  return xtc_coord(xd,natoms,box,x,prec);
}

