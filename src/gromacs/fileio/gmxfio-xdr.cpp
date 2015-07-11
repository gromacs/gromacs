/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gmxfio-xdr.h"

#include <cstdio>
#include <cstring>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "gmxfio-impl.h"

/* Enumerated for data types in files */
enum {
    eioREAL, eioFLOAT, eioDOUBLE, eioINT, eioINT64,
    eioUCHAR, eioNUCHAR, eioUSHORT,
    eioRVEC, eioNRVEC, eioIVEC, eioSTRING, eioNR
};

static const char *eioNames[eioNR] =
{
    "REAL", "INT", "GMX_STE_T", "UCHAR", "NUCHAR", "USHORT", "RVEC", "NRVEC",
    "IVEC", "STRING"
};

void gmx_fio_setprecision(t_fileio *fio, gmx_bool bDouble)
{
    gmx_fio_lock(fio);
    fio->bDouble = bDouble;
    gmx_fio_unlock(fio);
}

XDR *gmx_fio_getxdr(t_fileio *fio)
{
    XDR *ret = NULL;
    gmx_fio_lock(fio);
    GMX_RELEASE_ASSERT( fio->xdr != NULL, "Implementation error: NULL XDR pointers");
    ret = fio->xdr;
    gmx_fio_unlock(fio);
    return ret;
}

/* check the number of items given against the type */
static void gmx_fio_check_nitem(int eio, int nitem, const char *file, int line)
{
    if ((nitem != 1) && !((eio == eioNRVEC) || (eio == eioNUCHAR)))
    {
        gmx_fatal(FARGS,
                  "nitem (%d) may differ from 1 only for %s or %s, not   for %s"
                  "(%s, %d)", nitem, eioNames[eioNUCHAR], eioNames[eioNRVEC],
                  eioNames[eio], file, line);
    }
}

/* output a data type error. */
static void gmx_fio_fe(t_fileio *fio, int eio, const char *desc,
                       const char *srcfile, int line)
{
    gmx_fatal(FARGS, "Trying to %s %s type %d (%s), src %s, line %d",
              fio->bRead ? "read" : "write", desc, eio,
              ((eio >= 0) && (eio < eioNR)) ? eioNames[eio] : "unknown",
              srcfile, line);
}

/* This is the part that reads xdr files.  */

static gmx_bool do_xdr(t_fileio *fio, void *item, int nitem, int eio,
                       const char *desc, const char *srcfile, int line)
{
    unsigned char   ucdum, *ucptr;
    bool_t          res = 0;
    float           fvec[DIM];
    double          dvec[DIM];
    int             j, m, *iptr, idum;
    gmx_int64_t     sdum;
    real           *ptr;
    unsigned short  us;
    double          d = 0;
    float           f = 0;

    GMX_RELEASE_ASSERT( fio->xdr != NULL, "Implementation error: NULL XDR pointers");
    gmx_fio_check_nitem(eio, nitem, srcfile, line);
    switch (eio)
    {
        case eioREAL:
            if (fio->bDouble)
            {
                if (item && !fio->bRead)
                {
                    d = *((real *) item);
                }
                res = xdr_double(fio->xdr, &d);
                if (item)
                {
                    *((real *) item) = d;
                }
            }
            else
            {
                if (item && !fio->bRead)
                {
                    f = *((real *) item);
                }
                res = xdr_float(fio->xdr, &f);
                if (item)
                {
                    *((real *) item) = f;
                }
            }
            break;
        case eioFLOAT:
            if (item && !fio->bRead)
            {
                f = *((float *) item);
            }
            res = xdr_float(fio->xdr, &f);
            if (item)
            {
                *((float *) item) = f;
            }
            break;
        case eioDOUBLE:
            if (item && !fio->bRead)
            {
                d = *((double *) item);
            }
            res = xdr_double(fio->xdr, &d);
            if (item)
            {
                *((double *) item) = d;
            }
            break;
        case eioINT:
            if (item && !fio->bRead)
            {
                idum = *(int *) item;
            }
            res = xdr_int(fio->xdr, &idum);
            if (item)
            {
                *(int *) item = idum;
            }
            break;
        case eioINT64:
            if (item && !fio->bRead)
            {
                sdum = *(gmx_int64_t *) item;
            }
            res = xdr_int64(fio->xdr, &sdum);
            if (item)
            {
                *(gmx_int64_t *) item = sdum;
            }
            break;
        case eioUCHAR:
            if (item && !fio->bRead)
            {
                ucdum = *(unsigned char *) item;
            }
            res = xdr_u_char(fio->xdr, &ucdum);
            if (item)
            {
                *(unsigned char *) item = ucdum;
            }
            break;
        case eioNUCHAR:
            ucptr = (unsigned char *) item;
            res   = 1;
            for (j = 0; (j < nitem) && res; j++)
            {
                res = xdr_u_char(fio->xdr, &(ucptr[j]));
            }
            break;
        case eioUSHORT:
            if (item && !fio->bRead)
            {
                us = *(unsigned short *) item;
            }
            res = xdr_u_short(fio->xdr, (unsigned short *) &us);
            if (item)
            {
                *(unsigned short *) item = us;
            }
            break;
        case eioRVEC:
            if (fio->bDouble)
            {
                if (item && !fio->bRead)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        dvec[m] = ((real *) item)[m];
                    }
                }
                res = xdr_vector(fio->xdr, (char *) dvec, DIM,
                                 static_cast<unsigned int>(sizeof(double)),
                                 (xdrproc_t) xdr_double);
                if (item)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        ((real *) item)[m] = dvec[m];
                    }
                }
            }
            else
            {
                if (item && !fio->bRead)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        fvec[m] = ((real *) item)[m];
                    }
                }
                res = xdr_vector(fio->xdr, (char *) fvec, DIM,
                                 static_cast<unsigned int>(sizeof(float)),
                                 (xdrproc_t) xdr_float);
                if (item)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        ((real *) item)[m] = fvec[m];
                    }
                }
            }
            break;
        case eioNRVEC:
            ptr = NULL;
            res = 1;
            for (j = 0; (j < nitem) && res; j++)
            {
                if (item)
                {
                    ptr = ((rvec *) item)[j];
                }
                res = do_xdr(fio, ptr, 1, eioRVEC, desc, srcfile, line);
            }
            break;
        case eioIVEC:
            iptr = (int *) item;
            res  = 1;
            for (m = 0; (m < DIM) && res; m++)
            {
                if (item && !fio->bRead)
                {
                    idum = iptr[m];
                }
                res = xdr_int(fio->xdr, &idum);
                if (item)
                {
                    iptr[m] = idum;
                }
            }
            break;
        case eioSTRING:
        {
            char *cptr;
            int   slen;

            if (item)
            {
                if (!fio->bRead)
                {
                    slen = strlen((char *) item) + 1;
                }
                else
                {
                    slen = 0;
                }
            }
            else
            {
                slen = 0;
            }

            if (xdr_int(fio->xdr, &slen) <= 0)
            {
                gmx_fatal(FARGS, "wrong string length %d for string %s"
                          " (source %s, line %d)", slen, desc, srcfile, line);
            }
            if (!item && fio->bRead)
            {
                snew(cptr, slen);
            }
            else
            {
                cptr = (char *)item;
            }
            if (cptr)
            {
                res = xdr_string(fio->xdr, &cptr, slen);
            }
            else
            {
                res = 1;
            }
            if (!item && fio->bRead)
            {
                sfree(cptr);
            }
            break;
        }
        default:
            gmx_fio_fe(fio, eio, desc, srcfile, line);
    }

    return (res != 0);
}

/*******************************************************************
 *
 * READ/WRITE FUNCTIONS
 *
 *******************************************************************/

gmx_bool gmx_fio_writee_string(t_fileio *fio, const char *item,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    void    *it = (void*)item; /* ugh.. */
    gmx_fio_lock(fio);
    ret = do_xdr(fio, it, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_real(t_fileio *fio, real *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;

}

gmx_bool gmx_fio_doe_float(t_fileio *fio, float *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_double(t_fileio *fio, double *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_doe_gmx_bool(t_fileio *fio, gmx_bool *item,
                              const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;

    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        int itmp = 0;
        ret      = do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
        *item    = itmp;
    }
    else
    {
        int itmp = *item;
        ret      = do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int(t_fileio *fio, int *item,
                         const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int64(t_fileio *fio, gmx_int64_t *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioINT64, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_uchar(t_fileio *fio, unsigned char *item,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ushort(t_fileio *fio, unsigned short *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_rvec(t_fileio *fio, rvec *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ivec(t_fileio *fio, ivec *item,
                          const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_string(t_fileio *fio, char *item,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


/* Array reading & writing */

gmx_bool gmx_fio_ndoe_real(t_fileio *fio, real *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioREAL, desc,
                            srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioFLOAT, desc,
                            srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioDOUBLE, desc,
                            srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_gmx_bool(t_fileio *fio, gmx_bool *item, int n,
                               const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;

    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        if (fio->bRead)
        {
            int itmp = 0;
            ret      = ret && do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
            item[i]  = itmp;
        }
        else
        {
            int itmp = item[i];
            ret      = ret && do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioINT, desc,
                            srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioINT64, desc,
                            srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_uchar(t_fileio *fio, unsigned char *item, int n,
                            const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    ret = ret && do_xdr(fio, item, n, eioNUCHAR, desc,
                        srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioUSHORT, desc,
                            srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}



gmx_bool gmx_fio_ndoe_rvec(t_fileio *fio, rvec *item, int n,
                           const char *desc, const char *srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    ret = ret && do_xdr(fio, item, n, eioNRVEC, desc, srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioIVEC, desc,
                            srcfile, line);
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
        ret = ret && do_xdr(fio, &(item[i]), 1, eioSTRING, desc,
                            srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}
