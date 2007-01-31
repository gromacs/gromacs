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

#include <string.h>
/*#include "typedefs.h"*/
#include "smalloc.h"
/*#include "vec.h"
  #include "futil.h"*/
#include "xdrio/xdrfile.h"
#include "xdrio/xtcio.h"

#define XTC_MAGIC 1995

#ifdef HAVE_FSEEKO
#define gmx_fseek(A,B,C) fseeko(A,B,C)
#define gmx_ftell(A) ftello(A)
#define gmx_off_t off_t
#else
#define gmx_fseek(A,B,C) fseek(A,B,C)
#define gmx_ftell(A) ftell(A)
#define gmx_off_t int
#endif

#define xdrfile_int xdrfile_read_int

typedef struct {
  int code;
  char *msg;
} xtc_error_msg_t;

char *xtc_error_msg(int xtc_error)
{
  xtc_error_msg_t xem[] = {
    { -1, "Something is rotten in the state of Denmark" },
    { 0, "Everything is just fine. Just fine." }
  };
#define NXEM sizeof(xem)/sizeof(xem[0])
  int  i;
  char *ptr;
  
  for(i=0; (i<NXEM); i++) {
    if (xtc_error == xem[i].code)
      break;
  }
  if (i == NXEM)
    i == 0;
  ptr = strdup(xem[i]);
  
  return ptr;
}

static void check_xtc_magic(int magic)
{
  if (magic != XTC_MAGIC) 
    gmx_fatal(FARGS,"Magic Number Error in XTC file (read %d, should be %d)",
		magic,XTC_MAGIC);
}

int xtc_check(char *str,int bResult,char *file,int line)
{
  if (!bResult) {
    if (debug)
      fprintf(debug,"\nXTC error: read/write of %s failed, "
	      "source file %s, line %d\n",str,file,line);
    return 0;
  }
  return 1;
}

void xtc_check_fat_err(char *str,int bResult,char *file,int line)
{
  if (!bResult) {
    gmx_fatal(FARGS,"XTC read/write of %s failed, "
		"source file %s, line %d\n",str,file,line);
  }
}

static int xtc_header(XDRFILE *xd,int *magic,int *natoms,int *step,float *time,
		      int *bOK)
{
  int result;

  if (xdrfile_int(xd,magic) == 0)
    return 0;
  result=XTC_CHECK("natoms", xdrfile_int(xd,natoms));  /* number of atoms */
  if (result)
    result=XTC_CHECK("step",   xdrfile_int(xd,step));    /* frame number    */
  if (result)
    result=XTC_CHECK("time",   xdrfile_real(xd,time));   /* time            */
  *bOK=(result!=0);

  return result;
}

static int xtc_coord(XDRFILE *xd,int *natoms,matrix box,rvec *x,real *prec, int bRead)
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
      result=XTC_CHECK("box",xdr_real(xd,&(box[i][j])));
  
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



int write_xtc(XDRFILE *xdr,
	      int natoms,int step,real time,
	      matrix box,rvec *x,real prec)
{
  int magic_number = XTC_MAGIC;
  int bDum;

  /* write magic number and xtc identidier */
  if (!xtc_header(xdr,&magic_number,&natoms,&step,&time,&bDum))
    return 0;
    
  /* write data */
  return xtc_coord(xdr,&natoms,box,x,&prec,FALSE);
}

int read_first_xtc(XDRFILE *xdr,int *natoms,int *step,real *time,
		   matrix box,rvec **x,real *prec,int *bOK)
{
  int magic;
  
  *bOK=TRUE;
  
  /* read header and malloc x */
  if ( !xtc_header(xdr,&magic,natoms,step,time,bOK))
    return 0;
    
  /* Check magic number */
  check_xtc_magic(magic);
  
  snew(*x,*natoms);

  *bOK=xtc_coord(xdr,natoms,box,*x,prec,TRUE);
  
  return *bOK;
}

int read_next_xtc(XDRFILE *xdr,
		  int natoms,int *step,real *time,
		  matrix box,rvec *x,real *prec,int *bOK)
{
  int magic;
  int n;

  *bOK=TRUE;
  
  /* read header */
  if (!xtc_header(xdr,&magic,&n,step,time,bOK))
    return 0;
  if (n>natoms)
    gmx_fatal(FARGS, "Frame contains more atoms (%d) than expected (%d)", 
		n, natoms);
    
  /* Check magic number */
  check_xtc_magic(magic);

  *bOK=xtc_coord(xdr,&natoms,box,x,prec,TRUE);

  return *bOK;
}


/******************************************************************

  XTC files have a relatively simple structure.
  They have a header of 16 bytes and the rest are
  the compressed coordinates of the files. Due to the
  compression 00 is not present in the coordinates.
  The first 4 bytes of the header are the magic number
  1995 (0x000007CB). If we find this number we are guaranteed
  to be in the header, due to the presence of so many zeros.
  The second 4 bytes are the number of atoms in the frame, and is
  assumed to be constant. The third 4 bytes are the frame number.
  The last 4 bytes are a floating point representation of the time.

********************************************************************/

static const int header_size = 16;


/* Check if we are at the header start.
   At the same time it will also read 1 int */
static int xtc_at_header_start(XDRFILE *xdr, int natoms, int * timestep, float * time){
  int i_inp[3];
  float f_inp[10];
  int i;
  gmx_off_t off;


  if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
    return -1;
  }
  /* read magic natoms and timestep */
  for(i = 0;i<3;i++){
    if(!xdr_int(xdridptr[fp+1], &(i_inp[i]))){
      gmx_fseek(xdrfiles[fp+1],off+sizeof(int),SEEK_SET);
      return -1;
    }    
  }
  /* quick return */
  if(i_inp[0] != XTC_MAGIC){
    if(gmx_fseek(xdrfiles[fp+1],off+sizeof(int),SEEK_SET)){
      return -1;
    }
    return 0;
  }
  /* read time and box */
  for(i = 0;i<10;i++){
    if(!xdr_float(xdridptr[fp+1], &(f_inp[i]))){
      gmx_fseek(xdrfiles[fp+1],off+sizeof(int),SEEK_SET);
      return -1;
    }    
  }
  /* Make a rigourous check to see if we are in the beggining of a header
     Hopefully there are no ambiguous cases */
  /* This check makes use of the fact that the box matrix has 3 zeroes on the upper
     right triangle and that the first element must be nonzero unless the entire matrix is zero
  */
  if(i_inp[1] == natoms && 
     ((f_inp[1] != 0 && f_inp[6] == 0) ||
      (f_inp[1] == 0 && f_inp[5] == 0 && f_inp[9] == 0))){
    if(gmx_fseek(xdrfiles[fp+1],off+sizeof(int),SEEK_SET)){
      return -1;
    }
    *time = f_inp[0];
    *timestep = i_inp[2];
    return 1;
  }
  if(gmx_fseek(xdrfiles[fp+1],off+sizeof(int),SEEK_SET)){
    return -1;
  }
  return 0;         
}

static int 
xtc_get_next_frame_number(XDRFILE *xdr,int natoms)
{
    gmx_off_t off;
    int step;  
    float time;
    int ret;

    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }

    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdridptr[fp+1],&step);
    while(1){
      ret = xtc_at_header_start(fp,natoms,&step,&time);
      if(ret == 1){
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
	return step;
      }else if(ret == -1){
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
      }
    }
    return -1;
}


static float 
xtc_get_next_frame_time(XDRFILE *xdr,int natoms, int * bOK)
{
    gmx_off_t off;
    float time;
    int step;
    int ret;
    *bOK = 0;
    
    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdridptr[fp+1],&step);
    while(1){
      ret = xtc_at_header_start(fp,natoms,&step,&time);
      if(ret == 1){
	*bOK = 1;
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  *bOK = 0;
	  return -1;
	}
	return time;
      }else if(ret == -1){
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
	return -1;
      }
    }
    return -1;
}


static float 
xtc_get_current_frame_time(XDRFILE *xdr,int natoms, int * bOK)
{
    gmx_off_t off;
    int step;  
    float time;
    int ret;
    *bOK = 0;

    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }

    
    while(1){
      ret = xtc_at_header_start(fp,natoms,&step,&time);
      if(ret == 1){
	*bOK = 1;
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  *bOK = 0;
	  return -1;
	}
	return time;
      }else if(ret == -1){
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
	return -1;
      }else if(ret == 0){
	/*Go back.*/
	if(gmx_fseek(xdrfiles[fp+1],-8,SEEK_CUR)){
	  return -1;
	}
      }
    }
    return -1;
}

/* Currently not used, just for completeness */
static int 
xtc_get_current_frame_number(XDRFILE *xdr,int natoms, int * bOK)
{
    gmx_off_t off;
    int ret;  
    int step;
    float time;
    *bOK = 0;
    
    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }


    while(1){
      ret = xtc_at_header_start(fp,natoms,&step,&time);
      if(ret == 1){
	*bOK = 1;
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
	return step;
      }else if(ret == -1){
	if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
	  return -1;
	}
	return -1;
      }
    }
    return -1;
}


static gmx_off_t xtc_get_next_frame_start(XDRFILE *xdr, int natoms)
{
    int inp;
    gmx_off_t res;
    int ret;
    int step;
    float time;
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdridptr[fp+1],&step);
    while(1)
    {
      ret = xtc_at_header_start(fp,natoms,&step,&time);
      if(ret == 1){
	if((res = gmx_ftell(xdrfiles[fp+1])) >= 0){
	  return res - sizeof(int);
	}else{
	  return res;
	}
      }else if(ret == -1){
	return -1;
      }
    }
    return -1;
}



static float 
xtc_estimate_dt(XDRFILE *xdr, int natoms, int * bOK)
{
  float  res;
  float  tinit;
  gmx_off_t off;
  
  if((off   = gmx_ftell(xdrfiles[fp+1])) < 0){
    return -1;
  }
  
    tinit = xtc_get_current_frame_time(fp,natoms,bOK);
    
    *bOK = 1;
    
    if(!bOK)
    {
        return -1;
    }
    
    res = xtc_get_next_frame_time(fp,natoms,bOK);
    
    if(!bOK)
    {
        return -1;
    }
    
    res -= tinit;
    if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
      *bOK = 0;
      return -1;
    }
    return res;
}


int 
xtc_seek_frame(int frame, XDRFILE *xdr, int natoms)
{
    gmx_off_t low = 0;
    gmx_off_t high,pos;

    
    /* round to 4 bytes */
    int fr;
    gmx_off_t  offset;
    if(gmx_fseek(xdrfiles[fp+1],0,SEEK_END)){
      return -1;
    }

    if((high = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }
    
    /* round to 4 bytes  */
    high /= sizeof(int);
    high *= sizeof(int);
    offset = ((high/2)/sizeof(int))*sizeof(int);
    
    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return -1;
    }
    
    while(1)
    {
        fr = xtc_get_next_frame_number(fp,natoms);
        if(fr < 0)
        {
            return -1;
        }
        if(fr != frame && abs(low-high) > header_size)
        {
            if(fr < frame)
            {
                low = offset;      
            }
            else
            {
                high = offset;      
            }
            /* round to 4 bytes */
            offset = (((high+low)/2)/4)*4;
            
            if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
	      return -1;
	    }
        }
        else
        {
            break;
        }
    }
    if(offset <= header_size)
    {
        offset = low;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return -1;
    }

    if((pos = xtc_get_next_frame_start(fp,natoms))< 0){
    /* we probably hit an end of file */
      return -1;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],pos,SEEK_SET)){
      return -1;
    }
    
    return 0;
}

     

int 
xtc_seek_time(real time, XDRFILE *xdr, int natoms)
{
    float t;
    float dt;
    int bOK;
    gmx_off_t low = 0;
    gmx_off_t high,offset,pos;
    int res;
    int dt_sign = 0;
  
    if(gmx_fseek(xdrfiles[fp+1],0,SEEK_END)){
      return -1;
    }
    
    if((high = gmx_ftell(xdrfiles[fp+1])) < 0){
      return -1;
    }
    /* round to int  */
    high /= sizeof(int);
    high *= sizeof(int);
    offset = ((high/2)/sizeof(int))*sizeof(int);

    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return -1;	
   }

    dt = xtc_estimate_dt(fp,natoms,&bOK);
      
      
    
    if(!bOK)
    {
        return -1;
    }
    
    while(1)
    {
	dt = xtc_estimate_dt(fp,natoms,&bOK);
        if(!bOK)
        {
            return -1;
        }else{
	  if(dt > 0){
	    if(dt_sign == -1){
	      /* Found a place in the trajectory that has positive time step while
		 other has negative time step */
	      return -2;
	    }
	    dt_sign = 1;
	  }else if(dt < 0){
	    if(dt_sign == 1){
	      /* Found a place in the trajectory that has positive time step while
		 other has negative time step */
	      return -2;
	    }
	    dt_sign = -1;
	  }	  
	}
        t = xtc_get_next_frame_time(fp,natoms,&bOK);
        if(!bOK)
        {
            return -1;
        }

	/* If we are before the target time and the time step is positive or 0, or we have
	 after the target time and the time step is negative, or the difference between 
	the current time and the target time is bigger than dt and above all the distance between high
	and low is bigger than 1 frame, then do another step of binary search. Otherwise stop and check
	if we reached the solution */
        if((((t < time && dt_sign >= 0) || (t > time && dt_sign == -1)) || ((t-time) >= dt && dt_sign >= 0) || ((time-t) >= -dt && dt_sign < 0))  && (abs(low-high) > header_size))
        {
	  if(dt >= 0 && dt_sign != -1)
            {
                if(t < time)
                {
                    low = offset;      
                }
                else
                {
                    high = offset;      
                }
	    }
	  else if(dt <= 0 && dt_sign == -1)
            {
	      if(t >= time)
                {
		  low = offset;      
                }
	      else
                {
		  high = offset;      
                }
	    }else{
	      /* We should never reach here */
	      return -1;
	    }
            /* round to 4 bytes and subtract header*/
            offset = (((high+low)/2)/sizeof(int))*sizeof(int);
            if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
	      return -1;
	    }
        }
        else
        {
            if(abs(low-high) <= header_size)
            {
                break;
            }
            /* reestimate dt */
            if(xtc_estimate_dt(fp,natoms,&bOK) != dt)
            {
                if(bOK)
                {
                    dt = xtc_estimate_dt(fp,natoms,&bOK);
                }
            }
            if(t >= time && t-time < dt)
            {
                break;
            }
        }        
    }
    
    if(offset <= header_size)
    {
        offset = low;
    }
    
    gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET);

    if((pos = xtc_get_next_frame_start(fp,natoms)) < 0){
      return -1;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],pos,SEEK_SET)){
      return -1;
    }
    return 0;
}

float 
xtc_get_last_frame_time(XDRFILE *xdr, int natoms, int * bOK)
{
    float  time;
    gmx_off_t  off;
    int res;
    *bOK = 1;
    off = gmx_ftell(xdrfiles[fp+1]);  
    if(off < 0){
      *bOK = 0;
      return -1;
    }
    
    if( (res = gmx_fseek(xdrfiles[fp+1],-4,SEEK_END)) != 0){
      *bOK = 0;
      return -1;
    }

    time = xtc_get_current_frame_time(fp, natoms, bOK);
    if(!(*bOK)){
      return -1;
    }
    
    if( (res = gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)) != 0){
      *bOK = 0;
      return -1;
    } 
    return time;
}


int
xtc_get_last_frame_number(XDRFILE *xdr, int natoms, int * bOK)
{
    int    frame;
    gmx_off_t  off;
    int res;
    *bOK = 1;
    
    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      *bOK = 0;
      return -1;
    }

    
    if(gmx_fseek(xdrfiles[fp+1],-4,SEEK_END)){
      *bOK = 0;
      return -1;
    }

    frame = xtc_get_current_frame_number(fp, natoms, bOK);
    if(!bOK){
      return -1;
    }


    if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
      *bOK = 0;
      return -1;
    }    

    return frame;
}
  
