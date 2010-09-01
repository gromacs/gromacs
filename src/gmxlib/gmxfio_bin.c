/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include <stdio.h>
#include <errno.h>
#ifdef HAVE_IO_H
#include <io.h>
#endif

#include "gmx_fatal.h"
#include "macros.h"
#include "smalloc.h"
#include "futil.h"
#include "filenm.h"
#include "string2.h"
#include "gmxfio.h"
#include "md5.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

#include "gmxfio_int.h"

/* This is the part that reads dummy and ascii files.  */


static gmx_bool do_binread(t_fileio *fio, void *item, int nitem, int eio, 
                       const char *desc, const char *srcfile, int line);
static gmx_bool do_binwrite(t_fileio *fio, const void *item, int nitem, int eio, 
                        const char *desc, const char *srcfile, int line);


const t_iotype bin_iotype={do_binread, do_binwrite};


static gmx_bool do_binwrite(t_fileio *fio, const void *item, int nitem, int eio, 
                        const char *desc, const char *srcfile, int line)
{
    size_t size = 0, wsize;
    int ssize;

    gmx_fio_check_nitem(fio, eio, nitem, srcfile, line);
    switch (eio)
    {
    case eioREAL:
        size = sizeof(real);
        break;
    case eioFLOAT:
        size = sizeof(float);
        break;
    case eioDOUBLE:
        size = sizeof(double);
        break;
    case eioINT:
        size = sizeof(int);
        break;
    case eioGMX_LARGE_INT:
        size = sizeof(gmx_large_int_t);
        break;
    case eioUCHAR:
        size = sizeof(unsigned char);
        break;
    case eioNUCHAR:
        size = sizeof(unsigned char);
        break;
    case eioUSHORT:
        size = sizeof(unsigned short);
        break;
    case eioRVEC:
        size = sizeof(rvec);
        break;
    case eioNRVEC:
        size = sizeof(rvec);
        break;
    case eioIVEC:
        size = sizeof(ivec);
        break;
    case eioSTRING:
        size = ssize = strlen((char *) item) + 1;
        do_binwrite(fio, &ssize, 1, eioINT, desc, srcfile, line);
        break;
    default:
        gmx_fio_fe(fio, eio, desc, srcfile, line);
    }

    wsize = fwrite(item, size, nitem, fio->fp);

    if ((wsize != nitem) && fio->bDebug)
    {
        fprintf(stderr,
                "Error writing %s %s to file %s (source %s, line %d)\n",
                eioNames[eio], desc, fio->fn, srcfile, line);
        fprintf(stderr, "written size %u bytes, source size %u bytes\n",
                (unsigned int) wsize, (unsigned int) size);
    }
    return (wsize == nitem);
}

static gmx_bool do_binread(t_fileio *fio, void *item, int nitem, int eio, 
                       const char *desc, const char *srcfile, int line)
{
    size_t size = 0, rsize;
    int ssize;

    gmx_fio_check_nitem(fio, eio, nitem, srcfile, line);
    switch (eio)
    {
    case eioREAL:
        if (fio->bDouble)
            size = sizeof(double);
        else
            size = sizeof(float);
        break;
    case eioFLOAT:
        size = sizeof(float);
        break;
    case eioDOUBLE:
        size = sizeof(double);
        break;
    case eioINT:
        size = sizeof(int);
        break;
    case eioGMX_LARGE_INT:
        size = sizeof(gmx_large_int_t);
        break;
    case eioUCHAR:
        size = sizeof(unsigned char);
        break;
    case eioNUCHAR:
        size = sizeof(unsigned char);
        break;
    case eioUSHORT:
        size = sizeof(unsigned short);
        break;
    case eioRVEC:
    case eioNRVEC:
        if (fio->bDouble)
            size = sizeof(double) * DIM;
        else
            size = sizeof(float) * DIM;
        break;
    case eioIVEC:
        size = sizeof(ivec);
        break;
    case eioSTRING:
        do_binread(fio, &ssize, 1, eioINT, desc, srcfile, line);
        size = ssize;
        break;
    default:
        gmx_fio_fe(fio, eio, desc, srcfile, line);
    }
    if (item)
        rsize = fread(item, size, nitem, fio->fp);
    else
    {
        /* Skip over it if we have a NULL pointer here */
        gmx_fseek(fio->fp, (gmx_off_t)(size*nitem), SEEK_CUR);
        rsize = nitem;
    }
    if ((rsize != nitem) && (fio->bDebug))
        fprintf(stderr,
                "Error reading %s %s from file %s (source %s, line %d)\n",
                eioNames[eio], desc, fio->fn, srcfile, line);

    return (rsize == nitem);
}


