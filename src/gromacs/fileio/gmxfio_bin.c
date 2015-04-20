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
#include <string.h>

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

/* This is the part that reads dummy and ascii files.  */


static gmx_bool do_binread(t_fileio *fio, void *item, int nitem, int eio,
                           const char *desc, const char *srcfile, int line);
static gmx_bool do_binwrite(t_fileio *fio, const void *item, int nitem, int eio,
                            const char *desc, const char *srcfile, int line);


const t_iotype bin_iotype = {do_binread, do_binwrite};


static gmx_bool do_binwrite(t_fileio *fio, const void *item, int nitem, int eio,
                            const char *desc, const char *srcfile, int line)
{
    size_t size = 0, wsize;
    int    ssize;

    gmx_fio_check_nitem(eio, nitem, srcfile, line);
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
        case eioINT64:
            size = sizeof(gmx_int64_t);
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
    int    ssize;

    gmx_fio_check_nitem(eio, nitem, srcfile, line);
    switch (eio)
    {
        case eioREAL:
            if (fio->bDouble)
            {
                size = sizeof(double);
            }
            else
            {
                size = sizeof(float);
            }
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
        case eioINT64:
            size = sizeof(gmx_int64_t);
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
            {
                size = sizeof(double) * DIM;
            }
            else
            {
                size = sizeof(float) * DIM;
            }
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
    {
        rsize = fread(item, size, nitem, fio->fp);
    }
    else
    {
        /* Skip over it if we have a NULL pointer here */
        gmx_fseek(fio->fp, (gmx_off_t)(size*nitem), SEEK_CUR);
        rsize = nitem;
    }
    if ((rsize != nitem) && (fio->bDebug))
    {
        fprintf(stderr,
                "Error reading %s %s from file %s (source %s, line %d)\n",
                eioNames[eio], desc, fio->fn, srcfile, line);
    }

    return (rsize == nitem);
}
