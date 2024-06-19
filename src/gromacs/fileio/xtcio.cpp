/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "xtcio.h"

#include <cstdio>
#include <cstring>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/smalloc.h"

static int xdr_r2f(XDR* xdrs, real* r, gmx_bool gmx_unused bRead)
{
#if GMX_DOUBLE
    float f;
    int   ret;

    if (!bRead)
    {
        f = *r;
    }
    ret = xdr_float(xdrs, &f);
    if (bRead)
    {
        *r = f;
    }

    return ret;
#else
    return xdr_float(xdrs, static_cast<float*>(r));
#endif
}


t_fileio* open_xtc(const std::filesystem::path& fn, const char* mode)
{
    return gmx_fio_open(fn, mode);
}

void close_xtc(t_fileio* fio)
{
    gmx_fio_close(fio);
}

static void check_xtc_magic(int magic)
{
    if (magic != XTC_MAGIC && magic != XTC_NEW_MAGIC)
    {
        gmx_fatal(FARGS, "Magic Number Error in XTC file (read %d, should be %d or %d)", magic, XTC_MAGIC, XTC_NEW_MAGIC);
    }
}

static int xtc_check(const char* str, gmx_bool bResult, const char* file, int line)
{
    if (!bResult)
    {
        if (debug)
        {
            fprintf(debug,
                    "\nXTC error: read/write of %s failed, "
                    "source file %s, line %d\n",
                    str,
                    file,
                    line);
        }
        return 0;
    }
    return 1;
}

#define XTC_CHECK(s, b) xtc_check(s, b, __FILE__, __LINE__)

static int xtc_header(XDR* xd, int* magic, int* natoms, int64_t* step, real* time, gmx_bool bRead, gmx_bool* bOK)
{
    int result;

    if (xdr_int(xd, magic) == 0)
    {
        return 0;
    }
    result = XTC_CHECK("natoms", xdr_int(xd, natoms)); /* number of atoms */
    if (result)
    {
        /* Note that XTC wasn't defined to be extensible, so we can't
         * fix the fact that we used xdr_int for the step number,
         * which is defined to be signed and 32 bit. */
        int intStep = *step;
        result      = XTC_CHECK("step", xdr_int(xd, &intStep)); /* frame number    */
        *step       = intStep;
    }
    if (result)
    {
        result = XTC_CHECK("time", xdr_r2f(xd, time, bRead)); /* time */
    }
    *bOK = (result != 0);

    return result;
}

static int xtc_coord(XDR* xd, int* natoms, rvec* box, rvec* x, real* prec, int magic_number, gmx_bool bRead)
{
    int i, j, result;
#if GMX_DOUBLE
    float* ftmp;
    float  fprec;
#endif

    /* box */
    result = 1;
    for (i = 0; ((i < DIM) && result); i++)
    {
        for (j = 0; ((j < DIM) && result); j++)
        {
            result = XTC_CHECK("box", xdr_r2f(xd, &(box[i][j]), bRead));
        }
    }

    if (!result)
    {
        return result;
    }

#if GMX_DOUBLE
    /* allocate temp. single-precision array */
    snew(ftmp, static_cast<std::size_t>(*natoms) * DIM);

    /* Copy data to temp. array if writing */
    if (!bRead)
    {
        for (i = 0; (i < *natoms); i++)
        {
            ftmp[DIM * i + XX] = x[i][XX];
            ftmp[DIM * i + YY] = x[i][YY];
            ftmp[DIM * i + ZZ] = x[i][ZZ];
        }
        fprec = *prec;
    }
    result = XTC_CHECK("x", xdr3dfcoord(xd, ftmp, natoms, &fprec, magic_number));

    /* Copy from temp. array if reading */
    if (bRead)
    {
        for (i = 0; (i < *natoms); i++)
        {
            x[i][XX] = ftmp[DIM * i + XX];
            x[i][YY] = ftmp[DIM * i + YY];
            x[i][ZZ] = ftmp[DIM * i + ZZ];
        }
        *prec = fprec;
    }
    sfree(ftmp);
#else
    result = XTC_CHECK("x", xdr3dfcoord(xd, x[0], natoms, prec, magic_number));
#endif

    return result;
}


int write_xtc(t_fileio* fio, int natoms, int64_t step, real time, const rvec* box, const rvec* x, real prec)
{
    // By default we only write the new format for very large systems, but since the reading code
    // will adapt to whatever magic number is present in the header you could generate frames
    // for small systems that use the new format (which is useful for testing), and those should
    // be readable by normal implementations no matter how many atoms are present in the file.
    int      magic_number = (natoms > XTC_1995_MAX_NATOMS) ? XTC_NEW_MAGIC : XTC_MAGIC;
    XDR*     xd;
    gmx_bool bDum;
    int      bOK;

    if (!fio)
    {
        /* This means the fio object is not being used, e.g. because
           we are actually writing TNG output. We still have to return
           a pseudo-success value, to keep some callers happy. */
        return 1;
    }

    xd = gmx_fio_getxdr(fio);
    /* write magic number and xtc identidier */
    if (xtc_header(xd, &magic_number, &natoms, &step, &time, FALSE, &bDum) == 0)
    {
        return 0;
    }

    /* write data */
    bOK = xtc_coord(xd, &natoms, const_cast<rvec*>(box), const_cast<rvec*>(x), &prec, magic_number, FALSE); /* bOK will be 1 if writing went well */

    if (bOK)
    {
        if (gmx_fio_flush(fio) != 0)
        {
            bOK = 0;
        }
    }
    return bOK; /* 0 if bad, 1 if writing went well */
}

int read_first_xtc(t_fileio* fio, int* natoms, int64_t* step, real* time, matrix box, rvec** x, real* prec, gmx_bool* bOK)
{
    int  magic;
    XDR* xd;

    *bOK = TRUE;
    xd   = gmx_fio_getxdr(fio);

    /* read header and malloc x */
    if (!xtc_header(xd, &magic, natoms, step, time, TRUE, bOK))
    {
        return 0;
    }

    /* Check magic number */
    check_xtc_magic(magic);

    snew(*x, *natoms);

    *bOK = (xtc_coord(xd, natoms, box, *x, prec, magic, TRUE) != 0);

    return static_cast<int>(*bOK);
}

int read_next_xtc(t_fileio* fio, int natoms, int64_t* step, real* time, matrix box, rvec* x, real* prec, gmx_bool* bOK)
{
    int  magic;
    int  n;
    XDR* xd;

    *bOK = TRUE;
    xd   = gmx_fio_getxdr(fio);

    /* read header */
    if (!xtc_header(xd, &magic, &n, step, time, TRUE, bOK))
    {
        return 0;
    }

    /* Check magic number */
    check_xtc_magic(magic);

    if (n > natoms)
    {
        gmx_fatal(FARGS, "Frame contains more atoms (%d) than expected (%d)", n, natoms);
    }

    *bOK = (xtc_coord(xd, &natoms, box, x, prec, magic, TRUE) != 0);

    return static_cast<int>(*bOK);
}
