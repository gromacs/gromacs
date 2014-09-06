/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "config.h"

#include <errno.h>
#include <stdio.h>

#ifdef HAVE_IO_H
#include <io.h>
#endif

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_int.h"
#include "gromacs/fileio/md5.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


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
    ret = fio->iotp->nread(fio, item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_float(t_fileio *fio, float *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_reade_double(t_fileio *fio, double *item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_int(t_fileio *fio, int *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_int64(t_fileio *fio, gmx_int64_t *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioINT64, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_uchar(t_fileio *fio, unsigned char *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_reade_ushort(t_fileio *fio, unsigned short *item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_rvec(t_fileio *fio, rvec *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_ivec(t_fileio *fio, ivec *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_reade_string(t_fileio *fio, char *item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


/* Write */

gmx_bool gmx_fio_writee_real(t_fileio *fio, real item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_float(t_fileio *fio, float item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_double(t_fileio *fio, double item,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_writee_int(t_fileio *fio, int item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_int64(t_fileio *fio, gmx_int64_t item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioINT64, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_uchar(t_fileio *fio, unsigned char item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_ushort(t_fileio *fio, unsigned short item,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, &item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_rvec(t_fileio *fio, rvec *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_ivec(t_fileio *fio, ivec *item,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_writee_string(t_fileio *fio, const char *item,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, item, 1, eioSTRING, desc, srcfile, line);
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
    {
        ret = fio->iotp->nread(fio, item, 1, eioREAL, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioREAL, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;

}

gmx_bool gmx_fio_doe_float(t_fileio *fio, float *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioFLOAT, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioFLOAT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_double(t_fileio *fio, double *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_doe_gmx_bool(t_fileio *fio, gmx_bool *item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    int      itmp;

    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret   = fio->iotp->nread(fio, &itmp, 1, eioINT, desc, srcfile, line);
        *item = itmp;
    }
    else
    {
        itmp = *item;
        ret  = fio->iotp->nwrite(fio, &itmp, 1, eioINT, desc, srcfile, line);
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
    {
        ret = fio->iotp->nread(fio, item, 1, eioINT, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int64(t_fileio *fio, gmx_int64_t *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioINT64, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioINT64, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_uchar(t_fileio *fio, unsigned char *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioUCHAR, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioUCHAR, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ushort(t_fileio *fio, unsigned short *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioUSHORT, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioUSHORT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_rvec(t_fileio *fio, rvec *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioRVEC, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioRVEC, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ivec(t_fileio *fio, ivec *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioIVEC, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioIVEC, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_string(t_fileio *fio, char *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = fio->iotp->nread(fio, item, 1, eioSTRING, desc, srcfile, line);
    }
    else
    {
        ret = fio->iotp->nwrite(fio, item, 1, eioSTRING, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}







/* Array reading & writing */

gmx_bool gmx_fio_nreade_real(t_fileio *fio, real *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioREAL, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_float(t_fileio *fio, float *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioFLOAT, desc,
                                      srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nreade_double(t_fileio *fio, double *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioDOUBLE, desc,
                                      srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_int(t_fileio *fio, int *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_int64(t_fileio *fio, gmx_int64_t *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT64, desc,
                                      srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nreade_uchar(t_fileio *fio, unsigned char *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, n, eioNUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_ushort(t_fileio *fio, unsigned short *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioUSHORT, desc,
                                      srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_rvec(t_fileio *fio, rvec *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nread(fio, item, n, eioNRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_ivec(t_fileio *fio, ivec *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, item[i], 1, eioIVEC, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nreade_string(t_fileio *fio, char *item[], int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nread(fio, item[i], 1, eioSTRING, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}




/* Array writing */

gmx_bool gmx_fio_nwritee_real(t_fileio *fio, const real *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioREAL, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_float(t_fileio *fio, const float *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioFLOAT, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_double(t_fileio *fio, const double *item, int n,
                                const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioDOUBLE, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_int(t_fileio *fio, const int *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT, desc, srcfile,
                                       line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_int64(t_fileio *fio,
                               const gmx_int64_t *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT64,
                                       desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_uchar(t_fileio *fio, const unsigned char *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, item, n, eioNUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_ushort(t_fileio *fio, const unsigned short *item, int n,
                                const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioUSHORT, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nwritee_rvec(t_fileio *fio, const rvec *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = fio->iotp->nwrite(fio, item, n, eioNRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_nwritee_ivec(t_fileio *fio, const ivec *item, int n,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioIVEC, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_nwritee_string(t_fileio *fio, const char *item[], int n,
                                const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioSTRING, desc, srcfile,
                                       line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



/* array read/write functions */

gmx_bool gmx_fio_ndoe_real(t_fileio *fio, real *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioREAL, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioREAL, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_float(t_fileio *fio, float *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioFLOAT, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioFLOAT, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_double(t_fileio *fio, double *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioDOUBLE, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioDOUBLE, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_gmx_bool(t_fileio *fio, gmx_bool *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i, itmp;

    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &itmp, 1, eioINT, desc,
                                          srcfile, line);
            item[i] = itmp;
        }
        else
        {
            itmp = item[i];
            ret  = ret && fio->iotp->nwrite(fio, &itmp, 1, eioINT, desc,
                                            srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_ndoe_int(t_fileio *fio, int *item, int n,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_int64(t_fileio *fio, gmx_int64_t *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioINT64, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioINT64, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_uchar(t_fileio *fio, unsigned char *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = ret && fio->iotp->nread(fio, item, n, eioNUCHAR, desc,
                                      srcfile, line);
    }
    else
    {
        ret = ret && fio->iotp->nwrite(fio, item, n, eioNUCHAR, desc,
                                       srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_ushort(t_fileio *fio, unsigned short *item, int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioUSHORT, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioUSHORT, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_rvec(t_fileio *fio, rvec *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        ret = ret && fio->iotp->nread(fio, item, n, eioNRVEC, desc, srcfile, line);
    }
    else
    {
        ret = ret && fio->iotp->nwrite(fio, item, n, eioNRVEC, desc, srcfile,
                                       line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_ivec(t_fileio *fio, ivec *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioIVEC, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioIVEC, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_string(t_fileio *fio, char *item[], int n,
                             const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            ret = ret && fio->iotp->nread(fio, &(item[i]), 1, eioSTRING, desc,
                                          srcfile, line);
        }
        else
        {
            ret = ret && fio->iotp->nwrite(fio, &(item[i]), 1, eioSTRING, desc,
                                           srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}
