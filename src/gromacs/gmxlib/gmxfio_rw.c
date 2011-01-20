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


/*******************************************************************
 *
 * READ/WRITE FUNCTIONS 
 *
*******************************************************************/

gmx_bool gmx_fio_reade_real(t_fileio *fio, real *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_float(t_fileio *fio, float *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_reade_double(t_fileio *fio, double *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_int(t_fileio *fio, int *item,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_gmx_large_int(t_fileio *fio, gmx_large_int_t *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioGMX_LARGE_INT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_uchar(t_fileio *fio, unsigned char *item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_reade_ushort(t_fileio *fio, unsigned short *item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_rvec(t_fileio *fio, rvec *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_ivec(t_fileio *fio, ivec *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_string(t_fileio *fio, char *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


/* Write */

gmx_bool gmx_fio_writee_real(t_fileio *fio, real item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_float(t_fileio *fio, float item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_double(t_fileio *fio, double item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_writee_int(t_fileio *fio, int item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_gmx_large_int(t_fileio *fio, gmx_large_int_t item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioGMX_LARGE_INT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_uchar(t_fileio *fio, unsigned char item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_ushort(t_fileio *fio, unsigned short item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, &item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_rvec(t_fileio *fio, rvec *item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_ivec(t_fileio *fio, ivec *item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_string(t_fileio *fio, const char *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}



/* Read/write functions */

gmx_bool gmx_fio_doe_real(t_fileio *fio, real *item,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioREAL, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;

}

gmx_bool gmx_fio_doe_float(t_fileio *fio, float *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioFLOAT, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_double(t_fileio *fio, double *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_doe_gmx_bool(t_fileio *fio, gmx_bool *item,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    int itmp;
    
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret=fio->iotp->nread(fio, &itmp, 1, eioINT, desc, srcfile, line);
        *item = itmp;
    }
    else
    {
        itmp = *item;
        ret=fio->iotp->nwrite(fio, &itmp, 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int(t_fileio *fio, int *item,
                     const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioINT, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_gmx_large_int(t_fileio *fio, gmx_large_int_t *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioGMX_LARGE_INT, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioGMX_LARGE_INT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_uchar(t_fileio *fio, unsigned char *item,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioUCHAR, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ushort(t_fileio *fio, unsigned short *item,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioUSHORT, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_rvec(t_fileio *fio, rvec *item,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioRVEC, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ivec(t_fileio *fio, ivec *item,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioIVEC, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_string(t_fileio *fio, char *item,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
        ret=fio->iotp->nread(fio, item, 1, eioSTRING, desc, srcfile, line);
    else
        ret=fio->iotp->nwrite(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}







/* Array reading & writing */

gmx_bool gmx_fio_nreade_real(t_fileio *fio, real *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_float(t_fileio *fio, float *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret= ret && fio->iotp->nread(fio, &(item[i]), 1, eioFLOAT, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nreade_double(t_fileio *fio, double *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret= ret && fio->iotp->nread(fio, &(item[i]), 1, eioDOUBLE, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_int(t_fileio *fio, int *item, int n,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_gmx_large_int(t_fileio *fio, gmx_large_int_t *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioGMX_LARGE_INT, desc, 
                              srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nreade_uchar(t_fileio *fio, unsigned char *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, n, eioNUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_ushort(t_fileio *fio, unsigned short *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioUSHORT, desc, 
                                    srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_rvec(t_fileio *fio, rvec *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nread(fio, item, n, eioNRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_ivec(t_fileio *fio, ivec *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, item[i], 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_string(t_fileio *fio, char *item[], int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nread(fio, item[i], 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}




/* Array writing */

gmx_bool gmx_fio_nwritee_real(t_fileio *fio, const real *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioREAL, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_float(t_fileio *fio, const float *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioFLOAT, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_double(t_fileio *fio, const double *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioDOUBLE, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_int(t_fileio *fio, const int *item, int n,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT, desc, srcfile, 
                                     line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_gmx_large_int(t_fileio *fio, 
                               const gmx_large_int_t *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioGMX_LARGE_INT, 
                                     desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_uchar(t_fileio *fio, const unsigned char *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, item, n, eioNUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_ushort(t_fileio *fio, const unsigned short *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioUSHORT, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nwritee_rvec(t_fileio *fio, const rvec *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret=fio->iotp->nwrite(fio, item, n, eioNRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_ivec(t_fileio *fio, const ivec *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioIVEC, desc, 
                                     srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nwritee_string(t_fileio *fio, const char *item[], int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
        ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioSTRING, desc, srcfile, 
                               line);
    gmx_fio_unlock(fio);
    return ret;
}



/* array read/write functions */

gmx_bool gmx_fio_ndoe_real(t_fileio *fio, real *item, int n,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioREAL, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioREAL, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_float(t_fileio *fio, float *item, int n,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioFLOAT, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioFLOAT, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_double(t_fileio *fio, double *item, int n,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioDOUBLE, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioDOUBLE, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_gmx_bool(t_fileio *fio, gmx_bool *item, int n,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i,itmp;
    
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &itmp, 1, eioINT, desc, 
                                  srcfile, line);
            item[i] = itmp;
        }
        else
        {
            itmp = item[i];
            ret=ret && fio->iotp->nwrite(fio, &itmp, 1, eioINT, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_ndoe_int(t_fileio *fio, int *item, int n,
                      const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_gmx_large_int(t_fileio *fio, gmx_large_int_t *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioGMX_LARGE_INT, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioGMX_LARGE_INT, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_uchar(t_fileio *fio, unsigned char *item, int n,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret=ret && fio->iotp->nread(fio, item, n, eioNUCHAR, desc, 
                                    srcfile, line);
    }
    else
    {
        ret=ret && fio->iotp->nwrite(fio, item, n, eioNUCHAR, desc, 
                                     srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_ushort(t_fileio *fio, unsigned short *item, int n,
                        const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioUSHORT, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioUSHORT, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_rvec(t_fileio *fio, rvec *item, int n,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret=ret && fio->iotp->nread(fio, item, n, eioNRVEC, desc, srcfile, line);
    }
    else
    {
        ret=ret && fio->iotp->nwrite(fio, item, n, eioNRVEC, desc, srcfile, 
                               line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_ivec(t_fileio *fio, ivec *item, int n,
                       const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioIVEC, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioIVEC, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_string(t_fileio *fio, char *item[], int n,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret=TRUE;
    int i;
    gmx_fio_lock(fio);
    for(i=0;i<n;i++)
    {
        if (fio->bRead)
        {
            ret=ret && fio->iotp->nread(fio, &(item[i]), 1, eioSTRING, desc, 
                                  srcfile, line);
        }
        else
        {
            ret=ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioSTRING, desc, 
                                   srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}





