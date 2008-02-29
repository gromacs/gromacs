/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "typedefs.h"
#include "xdrf.h"
#include "gmxfio.h"
#include "xtcio.h"
#include "smalloc.h"
#include "vec.h"
#include "futil.h"
#include "gmx_fatal.h"

#define XTC_MAGIC 1995


static int xdr_r2f(XDR *xdrs,real *r,bool bRead)
{
#ifdef GMX_DOUBLE
    float f;
    int   ret;
    
    if (!bRead)
      f = *r;
    ret = xdr_float(xdrs,&f);
    if (bRead)
      *r = f;
    
    return ret;
#else
    return xdr_float(xdrs,(float *)r);
#endif
}


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
    gmx_fatal(FARGS,"Magic Number Error in XTC file (read %d, should be %d)",
		magic,XTC_MAGIC);
}

int xtc_check(char *str,bool bResult,char *file,int line)
{
  if (!bResult) {
    if (debug)
      fprintf(debug,"\nXTC error: read/write of %s failed, "
	      "source file %s, line %d\n",str,file,line);
    return 0;
  }
  return 1;
}

void xtc_check_fat_err(char *str,bool bResult,char *file,int line)
{
  if (!bResult) {
    gmx_fatal(FARGS,"XTC read/write of %s failed, "
		"source file %s, line %d\n",str,file,line);
  }
}

static int xtc_header(XDR *xd,int *magic,int *natoms,int *step,real *time,
		      bool bRead,bool *bOK)
{
  int result;

  if (xdr_int(xd,magic) == 0)
    return 0;
  result=XTC_CHECK("natoms", xdr_int(xd,natoms));  /* number of atoms */
  if (result)
    result=XTC_CHECK("step",   xdr_int(xd,step));    /* frame number    */
  if (result)
    result=XTC_CHECK("time",   xdr_r2f(xd,time,bRead));   /* time */
  *bOK=(result!=0);

  return result;
}

static int xtc_coord(XDR *xd,int *natoms,matrix box,rvec *x,real *prec, bool bRead)
{
  int i,j,result;
#ifdef GMX_DOUBLE
  float *ftmp;
  float fprec;
#endif
    
  /* box */
  result=1;
  for(i=0; ((i<DIM) && result); i++)
    for(j=0; ((j<DIM) && result); j++)
      result=XTC_CHECK("box",xdr_r2f(xd,&(box[i][j]),bRead));

  if (!result)
      return result;
  
#ifdef GMX_DOUBLE
  /* allocate temp. single-precision array */
  snew(ftmp,(*natoms)*DIM);
  
  /* Copy data to temp. array if writing */
  if(!bRead)
  {
      for(i=0; (i<*natoms); i++)
      {
          ftmp[DIM*i+XX]=x[i][XX];      
          ftmp[DIM*i+YY]=x[i][YY];      
          ftmp[DIM*i+ZZ]=x[i][ZZ];      
      }
      fprec = *prec;
  }
  result=XTC_CHECK("x",xdr3dfcoord(xd,ftmp,natoms,&fprec));
  
  /* Copy from temp. array if reading */
  if(bRead)
  {
      for(i=0; (i<*natoms); i++)
      {
          x[i][XX] = ftmp[DIM*i+XX];      
          x[i][YY] = ftmp[DIM*i+YY];      
          x[i][ZZ] = ftmp[DIM*i+ZZ];      
      }
      *prec = fprec;
  }  
  sfree(ftmp);
#else
    result=XTC_CHECK("x",xdr3dfcoord(xd,x[0],natoms,prec)); 
#endif 
    
  return result;
}



int write_xtc(int fp,
	      int natoms,int step,real time,
	      matrix box,rvec *x,real prec)
{
  int magic_number = XTC_MAGIC;
  XDR *xd;
  bool bDum;

  xd = fio_getxdr(fp);
  /* write magic number and xtc identidier */
  if (!xtc_header(xd,&magic_number,&natoms,&step,&time,FALSE,&bDum))
    return 0;
    
  /* write data */
  return xtc_coord(xd,&natoms,box,x,&prec,FALSE);
}

int read_first_xtc(int fp,int *natoms,int *step,real *time,
		   matrix box,rvec **x,real *prec,bool *bOK)
{
  int magic;
  XDR *xd;
  
  *bOK=TRUE;
  xd = fio_getxdr(fp);
  
  /* read header and malloc x */
  if ( !xtc_header(xd,&magic,natoms,step,time,TRUE,bOK))
    return 0;
    
  /* Check magic number */
  check_xtc_magic(magic);
  
  snew(*x,*natoms);

  *bOK=xtc_coord(xd,natoms,box,*x,prec,TRUE);
  
  return *bOK;
}

int read_next_xtc(int fp,
		  int natoms,int *step,real *time,
		  matrix box,rvec *x,real *prec,bool *bOK)
{
  int magic;
  int n;
  XDR *xd;

  *bOK=TRUE;
  xd = fio_getxdr(fp);
  
  /* read header */
  if (!xtc_header(xd,&magic,&n,step,time,TRUE,bOK))
    return 0;
  if (n>natoms)
    gmx_fatal(FARGS, "Frame contains more atoms (%d) than expected (%d)", 
		n, natoms);
    
  /* Check magic number */
  check_xtc_magic(magic);

  *bOK=xtc_coord(xd,&natoms,box,x,prec,TRUE);

  return *bOK;
}


