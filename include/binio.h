/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _binio_h
#define _binio_h

static char *SRCID_binio_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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


extern void _blockwrite(FILE *fp,int nelem,int size,void *data,
                        char *what,char *file,int line);

extern void _blockread(FILE *fp,int nelem,int size,void *data,
                       char *what,char *file,int line);

#endif	/* _binio_h */
