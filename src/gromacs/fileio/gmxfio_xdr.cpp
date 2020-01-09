/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

#include "gmxfio_xdr.h"

#include <cstddef>
#include <cstdio>
#include <cstring>

#include <limits>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "gmxfio_impl.h"

/* Enumerated for data types in files */
enum
{
    eioREAL,
    eioFLOAT,
    eioDOUBLE,
    eioINT,
    eioINT32,
    eioINT64,
    eioUCHAR,
    eioCHAR,
    eioNCHAR,
    eioNUCHAR,
    eioUSHORT,
    eioRVEC,
    eioNRVEC,
    eioIVEC,
    eioSTRING,
    eioOPAQUE,
    eioNR
};

static const char* eioNames[eioNR] = { "REAL",  "FLOAT", "DOUBLE", "INT",    "INT32",  "INT64",
                                       "UCHAR", "CHAR",  "NCHAR",  "NUCHAR", "USHORT", "RVEC",
                                       "NRVEC", "IVEC",  "STRING", "OPAQUE" };

void gmx_fio_setprecision(t_fileio* fio, gmx_bool bDouble)
{
    gmx_fio_lock(fio);
    fio->bDouble = bDouble;
    gmx_fio_unlock(fio);
}

bool gmx_fio_is_double(t_fileio* fio)
{
    bool isDouble = false;
    gmx_fio_lock(fio);
    isDouble = fio->bDouble;
    gmx_fio_unlock(fio);
    return isDouble;
}

XDR* gmx_fio_getxdr(t_fileio* fio)
{
    XDR* ret = nullptr;
    gmx_fio_lock(fio);
    GMX_RELEASE_ASSERT(fio->xdr != nullptr, "Implementation error: NULL XDR pointers");
    ret = fio->xdr;
    gmx_fio_unlock(fio);
    return ret;
}

/* check the number of items given against the type */
static void gmx_fio_check_nitem(int eio, std::size_t nitem, const char* file, int line)
{
    if ((nitem != 1)
        && !((eio == eioNRVEC) || (eio == eioNUCHAR) || (eio == eioNCHAR) || (eio == eioOPAQUE)))
    {
        gmx_fatal(FARGS,
                  "nitem may differ from 1 only for %s, %s, %s or %s, not for %s"
                  "(%s, %d)",
                  eioNames[eioNUCHAR], eioNames[eioNRVEC], eioNames[eioNCHAR], eioNames[eioOPAQUE],
                  eioNames[eio], file, line);
    }
}

/* output a data type error. */
[[noreturn]] static void gmx_fio_fe(t_fileio* fio, int eio, const char* desc, const char* srcfile, int line)
{
    gmx_fatal(FARGS, "Trying to %s %s type %d (%s), src %s, line %d", fio->bRead ? "read" : "write",
              desc, eio, ((eio >= 0) && (eio < eioNR)) ? eioNames[eio] : "unknown", srcfile, line);
}

/* This is the part that reads xdr files.  */

static gmx_bool
do_xdr(t_fileio* fio, void* item, std::size_t nitem, int eio, const char* desc, const char* srcfile, int line)
{
    unsigned char  ucdum, *ucptr;
    char           cdum, *cptr;
    bool_t         res = 0;
    float          fvec[DIM];
    double         dvec[DIM];
    int            m, *iptr, idum;
    int32_t        s32dum;
    int64_t        s64dum;
    real*          ptr;
    unsigned short us;
    double         d = 0;
    float          f = 0;

    GMX_RELEASE_ASSERT(fio->xdr != nullptr, "Implementation error: NULL XDR pointers");
    gmx_fio_check_nitem(eio, nitem, srcfile, line);
    switch (eio)
    {
        case eioREAL:
            if (fio->bDouble)
            {
                if (item && !fio->bRead)
                {
                    d = *(static_cast<real*>(item));
                }
                res = xdr_double(fio->xdr, &d);
                if (item)
                {
                    *(static_cast<real*>(item)) = d;
                }
            }
            else
            {
                if (item && !fio->bRead)
                {
                    f = *(static_cast<real*>(item));
                }
                res = xdr_float(fio->xdr, &f);
                if (item)
                {
                    *(static_cast<real*>(item)) = f;
                }
            }
            break;
        case eioFLOAT:
            if (item && !fio->bRead)
            {
                f = *(static_cast<float*>(item));
            }
            res = xdr_float(fio->xdr, &f);
            if (item)
            {
                *(static_cast<float*>(item)) = f;
            }
            break;
        case eioDOUBLE:
            if (item && !fio->bRead)
            {
                d = *(static_cast<double*>(item));
            }
            res = xdr_double(fio->xdr, &d);
            if (item)
            {
                *(static_cast<double*>(item)) = d;
            }
            break;
        case eioINT:
            if (item && !fio->bRead)
            {
                idum = *static_cast<int*>(item);
            }
            res = xdr_int(fio->xdr, &idum);
            if (item)
            {
                *static_cast<int*>(item) = idum;
            }
            break;
        case eioINT32:
            if (item && !fio->bRead)
            {
                s32dum = *static_cast<int32_t*>(item);
            }
            res = xdr_int32(fio->xdr, &s32dum);
            if (item)
            {
                *static_cast<int32_t*>(item) = s32dum;
            }
            break;
        case eioINT64:
            if (item && !fio->bRead)
            {
                s64dum = *static_cast<int64_t*>(item);
            }
            res = xdr_int64(fio->xdr, &s64dum);
            if (item)
            {
                *static_cast<int64_t*>(item) = s64dum;
            }
            break;
        case eioUCHAR:
            if (item && !fio->bRead)
            {
                ucdum = *static_cast<unsigned char*>(item);
            }
            res = xdr_u_char(fio->xdr, &ucdum);
            if (item)
            {
                *static_cast<unsigned char*>(item) = ucdum;
            }
            break;
        case eioCHAR:
            if (item && !fio->bRead)
            {
                cdum = *static_cast<char*>(item);
            }
            res = xdr_char(fio->xdr, &cdum);
            if (item)
            {
                *static_cast<char*>(item) = cdum;
            }
            break;
        case eioNCHAR:
            cptr = static_cast<char*>(item);
            GMX_RELEASE_ASSERT(nitem < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                               "The XDR interface cannot handle array lengths > 2^31");
            res = xdr_vector(fio->xdr, cptr, static_cast<int>(nitem),
                             static_cast<unsigned int>(sizeof(char)),
                             reinterpret_cast<xdrproc_t>(xdr_char));
            break;
        case eioNUCHAR:
            ucptr = static_cast<unsigned char*>(item);
            GMX_RELEASE_ASSERT(nitem < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                               "The XDR interface cannot handle array lengths > 2^31");
            res = xdr_vector(fio->xdr, reinterpret_cast<char*>(ucptr), static_cast<int>(nitem),
                             static_cast<unsigned int>(sizeof(unsigned char)),
                             reinterpret_cast<xdrproc_t>(xdr_u_char));
            break;
        case eioUSHORT:
            if (item && !fio->bRead)
            {
                us = *static_cast<unsigned short*>(item);
            }
            res = xdr_u_short(fio->xdr, &us);
            if (item)
            {
                *static_cast<unsigned short*>(item) = us;
            }
            break;
        case eioRVEC:
            if (fio->bDouble)
            {
                if (item && !fio->bRead)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        dvec[m] = (static_cast<real*>(item))[m];
                    }
                }
                res = xdr_vector(fio->xdr, reinterpret_cast<char*>(dvec), DIM,
                                 static_cast<unsigned int>(sizeof(double)),
                                 reinterpret_cast<xdrproc_t>(xdr_double));
                if (item)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        (static_cast<real*>(item))[m] = dvec[m];
                    }
                }
            }
            else
            {
                if (item && !fio->bRead)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        fvec[m] = (static_cast<real*>(item))[m];
                    }
                }
                res = xdr_vector(fio->xdr, reinterpret_cast<char*>(fvec), DIM,
                                 static_cast<unsigned int>(sizeof(float)),
                                 reinterpret_cast<xdrproc_t>(xdr_float));
                if (item)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        (static_cast<real*>(item))[m] = fvec[m];
                    }
                }
            }
            break;
        case eioNRVEC:
            ptr = nullptr;
            res = 1;
            for (std::size_t j = 0; j < nitem && res; j++)
            {
                if (item)
                {
                    ptr = (static_cast<rvec*>(item))[j];
                }
                res = static_cast<bool_t>(do_xdr(fio, ptr, 1, eioRVEC, desc, srcfile, line));
            }
            break;
        case eioIVEC:
            iptr = static_cast<int*>(item);
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
            char* cptr;
            int   slen;

            if (item)
            {
                if (!fio->bRead)
                {
                    slen = strlen(static_cast<char*>(item)) + 1;
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
                gmx_fatal(FARGS,
                          "wrong string length %d for string %s"
                          " (source %s, line %d)",
                          slen, desc, srcfile, line);
            }
            if (!item && fio->bRead)
            {
                snew(cptr, slen);
            }
            else
            {
                cptr = static_cast<char*>(item);
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
        case eioOPAQUE:
        {
            if (item == nullptr && nitem > 0)
            {
                gmx_fatal(FARGS, "Null pointer provided for non-zero length XDR opaque data.");
            }

            if (nitem > 0)
            {
                // We need to support very large opaque data objects although the default
                // XDR interface only uses integers for the size field, since gromacs-2020
                // e.g. embeds the entire TPR body as a single such object, which would break all
                // TPR files larger than 2GB unless we handle it as a special case.
                // To avoid inserting extra padding, we calculate the chunk size as:
                // - The max value of a signed integer + 1
                // - Subtract 4 (the XDR object size) to get a size within the range of the signed int.
                const std::size_t maxChunk =
                        static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1 - 4;

                for (res = 1; res > 0 && nitem > 0;)
                {
                    std::size_t thisChunk = std::min(maxChunk, nitem);
                    res = xdr_opaque(fio->xdr, reinterpret_cast<char*>(item), thisChunk);
                    nitem -= thisChunk;
                }
            }
            else
            {
                res = 1;
            }
            break;
        }
        default: gmx_fio_fe(fio, eio, desc, srcfile, line);
    }

    return (res != 0);
}

/*******************************************************************
 *
 * READ/WRITE FUNCTIONS
 *
 *******************************************************************/

gmx_bool gmx_fio_writee_string(t_fileio* fio, const char* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    void*    it = const_cast<char*>(item); /* ugh.. */
    gmx_fio_lock(fio);
    ret = do_xdr(fio, it, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_real(t_fileio* fio, real* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioREAL, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_float(t_fileio* fio, float* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioFLOAT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_double(t_fileio* fio, double* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioDOUBLE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_doe_gmx_bool(t_fileio* fio, gmx_bool* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;

    gmx_fio_lock(fio);
    if (fio->bRead)
    {
        int itmp = 0;
        ret      = do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
        *item    = (itmp != 0);
    }
    else
    {
        int itmp = static_cast<int>(*item);
        ret      = do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int(t_fileio* fio, int* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioINT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int32(t_fileio* fio, int32_t* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioINT32, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_int64(t_fileio* fio, int64_t* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioINT64, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_uchar(t_fileio* fio, unsigned char* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_char(t_fileio* fio, char* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ushort(t_fileio* fio, unsigned short* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioUSHORT, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_rvec(t_fileio* fio, rvec* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_ivec(t_fileio* fio, ivec* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioIVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_string(t_fileio* fio, char* item, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, item, 1, eioSTRING, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_doe_opaque(t_fileio* fio, char* data, std::size_t size, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret;
    gmx_fio_lock(fio);
    ret = do_xdr(fio, data, size, eioOPAQUE, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

/* Array reading & writing */

gmx_bool gmx_fio_ndoe_real(t_fileio* fio, real* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioREAL, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_float(t_fileio* fio, float* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioFLOAT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_double(t_fileio* fio, double* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioDOUBLE, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_gmx_bool(t_fileio* fio, gmx_bool* item, int n, const char* desc, const char* srcfile, int line)
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
            item[i]  = (itmp != 0);
        }
        else
        {
            int itmp = static_cast<int>(item[i]);
            ret      = ret && do_xdr(fio, &itmp, 1, eioINT, desc, srcfile, line);
        }
    }
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_ndoe_int(t_fileio* fio, int* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioINT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_int64(t_fileio* fio, int64_t* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioINT64, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_uchar(t_fileio* fio, unsigned char* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    ret = ret && do_xdr(fio, item, n, eioNUCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}

gmx_bool gmx_fio_ndoe_char(t_fileio* fio, char* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    ret = ret && do_xdr(fio, item, n, eioNCHAR, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_ushort(t_fileio* fio, unsigned short* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioUSHORT, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_rvec(t_fileio* fio, rvec* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    gmx_fio_lock(fio);
    ret = ret && do_xdr(fio, item, n, eioNRVEC, desc, srcfile, line);
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_ivec(t_fileio* fio, ivec* item, int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioIVEC, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}


gmx_bool gmx_fio_ndoe_string(t_fileio* fio, char* item[], int n, const char* desc, const char* srcfile, int line)
{
    gmx_bool ret = TRUE;
    int      i;
    gmx_fio_lock(fio);
    for (i = 0; i < n; i++)
    {
        ret = ret && do_xdr(fio, &(item[i]), 1, eioSTRING, desc, srcfile, line);
    }
    gmx_fio_unlock(fio);
    return ret;
}

namespace gmx
{

FileIOXdrSerializer::FileIOXdrSerializer(t_fileio* fio) : fio_(fio)
{
    GMX_RELEASE_ASSERT(fio, "Need valid file io handle");
}

bool FileIOXdrSerializer::reading() const
{
    return fio_->bRead;
}

void FileIOXdrSerializer::doBool(bool* value)
{
    gmx_fio_do_gmx_bool(fio_, *value);
}

void FileIOXdrSerializer::doUChar(unsigned char* value)
{
    gmx_fio_do_uchar(fio_, *value);
}

void FileIOXdrSerializer::doChar(char* value)
{
    gmx_fio_do_char(fio_, *value);
}

void FileIOXdrSerializer::doUShort(unsigned short* value)
{
    gmx_fio_do_ushort(fio_, *value);
}

void FileIOXdrSerializer::doInt(int* value)
{
    gmx_fio_do_int(fio_, *value);
}

void FileIOXdrSerializer::doInt32(int32_t* value)
{
    gmx_fio_do_int32(fio_, *value);
}

void FileIOXdrSerializer::doInt64(int64_t* value)
{
    gmx_fio_do_int64(fio_, *value);
}

void FileIOXdrSerializer::doFloat(float* value)
{
    gmx_fio_do_float(fio_, *value);
}

void FileIOXdrSerializer::doDouble(double* value)
{
    gmx_fio_do_double(fio_, *value);
}

void FileIOXdrSerializer::doReal(real* value)
{
    gmx_fio_do_real(fio_, *value);
}

void FileIOXdrSerializer::doIvec(ivec* value)
{
    gmx_fio_do_ivec(fio_, *value);
}

void FileIOXdrSerializer::doRvec(rvec* value)
{
    gmx_fio_do_rvec(fio_, *value);
}

void FileIOXdrSerializer::doCharArray(char* values, int elements)
{
    gmx_fio_ndo_char(fio_, values, elements);
}

void FileIOXdrSerializer::doUCharArray(unsigned char* values, int elements)
{
    gmx_fio_ndo_uchar(fio_, values, elements);
}

void FileIOXdrSerializer::doRvecArray(rvec* values, int elements)
{
    gmx_fio_ndo_rvec(fio_, values, elements);
}

void FileIOXdrSerializer::doString(std::string* value)
{
    // TODO: Use an arbitrary length buffer (but that is not supported in
    // gmx_fio, either).
    char buf[STRLEN];
    if (!fio_->bRead)
    {
        std::strncpy(buf, value->c_str(), STRLEN);
        buf[STRLEN - 1] = 0;
    }
    gmx_fio_do_string(fio_, buf);
    if (fio_->bRead)
    {
        *value = buf;
    }
}

void FileIOXdrSerializer::doOpaque(char* data, std::size_t size)
{
    gmx_fio_do_opaque(fio_, data, size);
}

} // namespace gmx
