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

/* This is the part that reads xdr files.  */


/* file type functions */
static gmx_bool do_xdrread(t_fileio *fio, void *item, int nitem, int eio, 
                       const char *desc, const char *srcfile, int line);
static gmx_bool do_xdrwrite(t_fileio *fio, const void *item, int nitem, int eio, 
                        const char *desc, const char *srcfile, int line);


const t_iotype xdr_iotype={do_xdrread, do_xdrwrite};


#ifdef USE_XDR

static gmx_bool do_xdr(t_fileio *fio, void *item, int nitem, int eio, 
                   const char *desc, const char *srcfile, int line)
{
    unsigned char ucdum, *ucptr;
    bool_t res = 0;
    float fvec[DIM];
    double dvec[DIM];
    int j, m, *iptr, idum;
    gmx_large_int_t sdum;
    real *ptr;
    unsigned short us;
    double d = 0;
    float f = 0;

    gmx_fio_check_nitem(fio, eio, nitem, srcfile, line);
    switch (eio)
    {
    case eioREAL:
        if (fio->bDouble)
        {
            if (item && !fio->bRead)
                d = *((real *) item);
            res = xdr_double(fio->xdr, &d);
            if (item)
                *((real *) item) = d;
        }
        else
        {
            if (item && !fio->bRead)
                f = *((real *) item);
            res = xdr_float(fio->xdr, &f);
            if (item)
                *((real *) item) = f;
        }
        break;
    case eioFLOAT:
        if (item && !fio->bRead)
            f = *((float *) item);
        res = xdr_float(fio->xdr, &f);
        if (item)
            *((float *) item) = f;
        break;
    case eioDOUBLE:
        if (item && !fio->bRead)
            d = *((double *) item);
        res = xdr_double(fio->xdr, &d);
        if (item)
            *((double *) item) = d;
        break;
    case eioINT:
        if (item && !fio->bRead)
            idum = *(int *) item;
        res = xdr_int(fio->xdr, &idum);
        if (item)
            *(int *) item = idum;
        break;
    case eioGMX_LARGE_INT:
        /* do_xdr will not generate a warning when a 64bit gmx_large_int_t
         * value that is out of 32bit range is read into a 32bit gmx_large_int_t.
         */
        if (item && !fio->bRead)
            sdum = *(gmx_large_int_t *) item;
        res = xdr_gmx_large_int(fio->xdr, &sdum, NULL);
        if (item)
            *(gmx_large_int_t *) item = sdum;
        break;
    case eioUCHAR:
        if (item && !fio->bRead)
            ucdum = *(unsigned char *) item;
        res = xdr_u_char(fio->xdr, &ucdum);
        if (item)
            *(unsigned char *) item = ucdum;
        break;
    case eioNUCHAR:
        ucptr = (unsigned char *) item;
        res = 1;
        for (j = 0; (j < nitem) && res; j++)
        {
            res = xdr_u_char(fio->xdr, &(ucptr[j]));
        }
        break;
    case eioUSHORT:
        if (item && !fio->bRead)
            us = *(unsigned short *) item;
        res = xdr_u_short(fio->xdr, (unsigned short *) &us);
        if (item)
            *(unsigned short *) item = us;
        break;
    case eioRVEC:
        if (fio->bDouble)
        {
            if (item && !fio->bRead)
                for (m = 0; (m < DIM); m++)
                    dvec[m] = ((real *) item)[m];
            res = xdr_vector(fio->xdr, (char *) dvec, DIM,
                             (unsigned int) sizeof(double),
                             (xdrproc_t) xdr_double);
            if (item)
                for (m = 0; (m < DIM); m++)
                    ((real *) item)[m] = dvec[m];
        }
        else
        {
            if (item && !fio->bRead)
                for (m = 0; (m < DIM); m++)
                    fvec[m] = ((real *) item)[m];
            res = xdr_vector(fio->xdr, (char *) fvec, DIM,
                             (unsigned int) sizeof(float),
                             (xdrproc_t) xdr_float);
            if (item)
                for (m = 0; (m < DIM); m++)
                    ((real *) item)[m] = fvec[m];
        }
        break;
    case eioNRVEC:
        ptr = NULL;
        res = 1;
        for (j = 0; (j < nitem) && res; j++)
        {
            if (item)
                ptr = ((rvec *) item)[j];
            res = do_xdr(fio, ptr, 1, eioRVEC, desc, srcfile, line);
        }
        break;
    case eioIVEC:
        iptr = (int *) item;
        res = 1;
        for (m = 0; (m < DIM) && res; m++)
        {
            if (item && !fio->bRead)
                idum = iptr[m];
            res = xdr_int(fio->xdr, &idum);
            if (item)
                iptr[m] = idum;
        }
        break;
    case eioSTRING:
    {
        char *cptr;
        int slen;

        if (item)
        {
            if (!fio->bRead)
                slen = strlen((char *) item) + 1;
            else
                slen = 0;
        }
        else
            slen = 0;

        if (xdr_int(fio->xdr, &slen) <= 0)
            gmx_fatal(FARGS, "wrong string length %d for string %s"
                      " (source %s, line %d)",slen,desc,srcfile,line);
        if (!item && fio->bRead)
            snew(cptr,slen);
        else
            cptr=(char *)item;
        if (cptr)
            res = xdr_string(fio->xdr,&cptr,slen);
        else
            res = 1;
        if (!item && fio->bRead)
            sfree(cptr);
        break;
    }
    default:
        gmx_fio_fe(fio, eio, desc, srcfile, line);
    }
    if ((res == 0) && (fio->bDebug))
        fprintf(stderr,"Error in xdr I/O %s %s to file %s (source %s, line %d)\n",
                eioNames[eio],desc,fio->fn,srcfile,line);

    return (res != 0);
}


static gmx_bool do_xdrread(t_fileio *fio, void *item, int nitem, int eio, 
                       const char *desc, const char *srcfile, int line)
{
    return do_xdr(fio, item, nitem, eio, desc, srcfile, line);
}


static gmx_bool do_xdrwrite(t_fileio *fio, const void *item, int nitem, int eio, 
                        const char *desc, const char *srcfile, int line)
{
    void *it=(void*)item; /* ugh.. */
    return do_xdr(fio, it, nitem, eio, desc, srcfile, line);
}

#endif


