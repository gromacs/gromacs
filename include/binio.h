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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _binio_h
#define _binio_h

static char *SRCID_binio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) binio.h 1.5 11/23/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "sysstuff.h"
#include "fatal.h"

/*
 *    Module to binary write and read.
 *
 *                                        @                   
 *    @@@   @@                            @
 *     @     @
 *     @  @  @   @@@   @@ @@@  @@ @@    @@@    @@ @@     @@ @@
 *     @  @  @  @   @   @@   @  @@  @     @     @@  @   @  @@
 *     @ @ @ @      @   @       @   @     @     @   @  @    @
 *      @@ @@   @@@@@   @       @   @     @     @   @  @    @
 *      @   @  @    @   @       @   @     @     @   @  @    @
 *      @   @  @   @@   @       @   @     @     @   @   @  @@
 *      @   @   @@@ @@ @@@@    @@@ @@@  @@@@@  @@@ @@@   @@ @
 *                                                          @
 *                                                         @
 *                                                      @@@
 *
 *    Use this module only to write and read simple types or array(s)
 *    of simple types. STRUCTURES ARE DEFINITELY NOT ALLOWED.
 */

#define nblockwrite(fp,nelem,data) \
  _blockwrite(fp,nelem,sizeof(*data),(data),#data,__FILE__,__LINE__)
#define blockwrite(fp,data) \
  _blockwrite(fp,1,sizeof(data),&(data),#data,__FILE__,__LINE__)
#define cblockwrite(fp,ptr,nchars) \
  _blockwrite(fp,1,(nchars),(ptr),#ptr,__FILE__,__LINE__)
#define nblockread(fp,nelem,data) \
  _blockread(fp,nelem,sizeof(*data),(data),#data,__FILE__,__LINE__)
#define blockread(fp,data) \
  _blockread(fp,1,sizeof(data),&(data),#data,__FILE__,__LINE__)
#define cblockread(fp,ptr,nchars) \
  _blockread(fp,1,(nchars),(ptr),#ptr,__FILE__,__LINE__)

#define patch(fp,fpos,write) \
  do \
    { \
      int result,fhere; \
      \
      fhere=ftell(fp); \
      if ((result=fseek(fp,fpos,SEEK_SET))!=0) \
        fatal_error(errno,"could not seek to position %d from file %s, " \
                    "line %d, result=%d",(fpos),__FILE__,__LINE__,result); \
      write; \
      if ((result=fseek(fp,fhere,SEEK_SET))!=0) \
        fatal_error(errno,"could not seek back to %d from file %s, line %d," \
                    " result=%d",fhere,__FILE__,__LINE__,result); \
    } \
  while (0)

extern void _blockwrite(FILE *fp,int nelem,int size,void *data,
                        char *what,char *file,int line);

extern void _blockread(FILE *fp,int nelem,int size,void *data,
                       char *what,char *file,int line);

#endif	/* _binio_h */
