/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003-2007.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * IN contrast to the rest of Gromacs, XDRFILE is distributed under the 
 * BSD license, so you can use it any way you wish, including closed source:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a 
 * copy of this software and associated documentation files (the "Software"), 
 * to deal in the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 */


/* The external C module interface is documented in xdrfile.h, and the FORTRAN
 * interface in xdrfile_fortran.txt.
 *
 * Macros that will affect the generated code when defined:
 *
 * F77_FUNC         Macro to define name mangling for the Fortran77 interface.
 *                  By default we keep the names in lower case and append an
 *                  underscore. Alternative name mangling can be defined like: 
 *
 *                  #define F77_FUNC(name,NAME) _ ## NAME
 *                  (Use upper case, prepend an underscore)
 *
 * HAVE_RPC_XDR_H   Use system RPC/XDR libraries instead of our built-in 
 *                  implementation at the end of this file. This is useful if 
 *                  you have hardware with non-IEEE floating-point, or use 
 *		            EBCDIC instead of ASCII strings.
 *                  Obviously it wont hurt to enable in other cases too.
 *                  On some systems you might have to link with -lnsl
 *
 * Please note: The library is only threadsafe when called from C, not Fortran.
 */


/* Get HAVE_RPC_XDR_H, F77_FUNC from config.h if available */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>


/* get fixed-width types if we are using ANSI C99 */
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#elif (defined HAVE_INTTYPES_H)
#  include <inttypes.h>
#endif


#ifdef HAVE_RPC_XDR_H
#  include <rpc/rpc.h>
#  include <rpc/xdr.h>
#endif

#include "xdrfile.h"

/* This program tests reading and writing to XDR files */

static void _die(char *msg,int line,char *file)
{
  fprintf(stderr,"Fatal error: %s\n",msg);
  fprintf(stderr,"Death occurred at %s, line %d\n",file,line);
  exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

int main(int argc,char *argv[])
{
#define BUFLEN 37
  XDRFILE *xfp;
  int     i,j,k,len,ncoord=BUFLEN/3;
  char    ptr[BUFLEN],*buf = "abcdefghijklmnopqrstuvwxyz";
  unsigned char uptr[BUFLEN];
  short   sptr[BUFLEN],sptr2[BUFLEN];
  unsigned short usptr[BUFLEN],usptr2[BUFLEN];
  int     iptr[BUFLEN],iptr2[BUFLEN];
  unsigned int uiptr[BUFLEN],uiptr2[BUFLEN];
  float   fptr[BUFLEN],fptr2[BUFLEN];
  double  dptr[BUFLEN],dptr2[BUFLEN];
  char    optr[BUFLEN],optr2[BUFLEN];
#define NPREC 4
  float   fprec[] = { 0, 10, 243, 1003 };
  double  dprec[] = { -1, 13, 245, 998 };
  
  printf("Going to test the xdrfile library routines\n");

  /* Can not write a string that's on the stack since all data is
     treated as variables.
  */
  len = strlen(buf)+1;
  if (len >= BUFLEN)
    die("Increase BUFLEN");
  strcpy(ptr,buf);
  strcpy((char *)uptr,buf);
  /* Initiate float arrays */
  for(i=0; (i<BUFLEN); i++) {
    fptr[i] = cos(i*13.0/M_PI);
    dptr[i] = sin(i*13.0/M_PI);
  }
  /* Initiate opaque array */
  memcpy(optr,dptr,BUFLEN);
  
  /*************************************/
  /*           WRITING BIT             */   
  /*************************************/
  printf("Writing xdrfile\n");

  if ((xfp = xdrfile_open("test.xdr","w")) == NULL)
    die("Can not open file for writing");
    
  if (xdrfile_write_char(ptr,len,xfp) != len)
    die("Writing char string");
  if (xdrfile_write_uchar(uptr,len,xfp) != len)
    die("Writing uchar string");
  if (xdrfile_write_short(sptr,BUFLEN,xfp) != BUFLEN) 
    die("Writing short array");
  if (xdrfile_write_ushort(usptr,BUFLEN,xfp) != BUFLEN) 
    die("Writing ushort array");
  if (xdrfile_write_int(iptr,BUFLEN,xfp) != BUFLEN) 
    die("Writing int array");
  if (xdrfile_write_uint(uiptr,BUFLEN,xfp) != BUFLEN) 
    die("Writing uint array");
  if (xdrfile_write_float(fptr,BUFLEN,xfp) != BUFLEN)
    die("Writing float array");
  if (xdrfile_write_double(dptr,BUFLEN,xfp) != BUFLEN)
    die("Writing double array");
  if (xdrfile_write_string(buf,xfp) != len)
    die("Writing string");
  if (xdrfile_write_opaque(optr,BUFLEN,xfp) != BUFLEN)
    die("Writing opaque");
  for(k=0; (k<NPREC); k++) {
    if (xdrfile_compress_coord_float(fptr,ncoord,fprec[k],xfp) != ncoord)
      die("Writing compress_coord_float");
    if (xdrfile_compress_coord_double(dptr,ncoord,dprec[k],xfp) != ncoord)
      die("Writing compress_coord_double");
  }
  if (xdrfile_close(xfp) != 0)
    die("Can not close xdr file");
  
  /*************************************/
  /*          READING BIT              */
  /*************************************/
  printf("Reading xdrfile\n");
  if ((xfp = xdrfile_open("test.xdr","r")) == NULL)
    die("Can not open file for reading");
  
  if ((xdrfile_read_char(ptr,len,xfp)) != len)
    die("Not the right number of chars read from string");
  if (strcmp(ptr,buf) != 0)
    printf("did not read the expected chars");
  if (xdrfile_read_uchar(uptr,len,xfp) != len)
    die("Not the right number of uchars read from string");
  if (strcmp((char *)uptr,buf) != 0)
    printf("did not read the expected uchars");
  if (xdrfile_read_short(sptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading short array");
  for(i=0; (i<BUFLEN); i++) 
    if (sptr2[i] != sptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %10d, read: %10d\n",i,sptr[i],sptr2[i]);
      die("Comparing short array");
    }
  if (xdrfile_read_ushort(usptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading ushort array");
  for(i=0; (i<BUFLEN); i++) 
    if (usptr2[i] != usptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %10d, read: %10d\n",i,usptr[i],usptr2[i]);
      die("Comparing ushort array");
    }
  if (xdrfile_read_int(iptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading int array");
  for(i=0; (i<BUFLEN); i++) 
    if (iptr2[i] != iptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %10d, read: %10d\n",i,iptr[i],iptr2[i]);
      die("Comparing int array");
    }
  if (xdrfile_read_uint(uiptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading uint array");
  for(i=0; (i<BUFLEN); i++) 
    if (uiptr2[i] != uiptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %10d, read: %10d\n",i,uiptr[i],uiptr2[i]);
      die("Comparing uint array");
    }
  if (xdrfile_read_float(fptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading float array");
  for(i=0; (i<BUFLEN); i++) 
    if (fptr2[i] != fptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %12g, read: %12g\n",i,fptr[i],fptr2[i]);
      die("Comparing float array");
    }
  if (xdrfile_read_double(dptr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading double array");
  for(i=0; (i<BUFLEN); i++) 
    if (dptr2[i] != dptr[i]) {
      fprintf(stderr,"i: %5d, wrote: %12g, read: %12g\n",i,dptr[i],dptr2[i]);
      die("Comparing double array");
    }
  if (xdrfile_read_string(ptr,BUFLEN,xfp) != len)
    die("Reading string");
  if (strcmp(ptr,buf) != 0)
    die("Comparing strings");
  if (xdrfile_read_opaque(optr2,BUFLEN,xfp) != BUFLEN) 
    die("Reading opaque array");
  for(i=0; (i<BUFLEN); i++) 
    if (optr2[i] != optr[i]) {
      fprintf(stderr,"i: %5d, wrote: %2d, read: %2d\n",i,optr[i],optr2[i]);
      die("Comparing opaque array");
    }
  for(k=0; (k<NPREC); k++) {
    float  ff,fx;
    double dd,dx;
    int    nc=ncoord;
    if (xdrfile_decompress_coord_float(fptr2,&nc,&ff,xfp) != ncoord)
      die("Reading compress_coord_float");
    if (ff != fprec[k])
      die("Float precision");
    if (ff <= 0)
      ff = 1000;
    for(i=0; (i<ncoord); i++) 
      for(j=0; (j<3); j++) {
	fx = rint(fptr[3*i+j]/ff) * ff;
	if (fx != fptr2[3*i+j]) {
	  printf("prec: %10g, i: %3d, j: %d, fx: %10g, fptr2: %12g, fptr: %12g\n",
		 ff,i,j,fx,fptr2[3*i+j],fptr[3*i+j]);
	  die("Reading decompressed float coordinates");
	}
      }
    if (xdrfile_decompress_coord_double(dptr2,&nc,&dd,xfp) != ncoord)
      die("Reading compress_coord_double");
    if (dd != dprec[i])
      die("Double precision");
  }
    
    
  if (xdrfile_close(xfp) != 0)
    die("Can not close xdr file");
  
  printf("No errors\n");
  
  return 0;
}
