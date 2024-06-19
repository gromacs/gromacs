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

#include <cstdint>

#include "gromacs/fileio/xdrf.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

int xdr_real(XDR* xdrs, real* r)
{
#if GMX_DOUBLE
    float f;
    int   ret;

    f   = *r;
    ret = xdr_float(xdrs, &f);
    *r  = f;

    return ret;
#else
    return xdr_float(xdrs, static_cast<float*>(r));
#endif
}

int xdr3drcoord(XDR* xdrs, real* fp, int* size, real* precision, int magic_number)
{
#if GMX_DOUBLE
    float* ffp;
    float  fprec;
    int    i, ret, isize;

    isize = *size * DIM;
    if (isize <= 0)
    {
        gmx_fatal(FARGS, "Don't know what to malloc for ffp, isize = %d", isize);
    }

    snew(ffp, isize);

    for (i = 0; (i < isize); i++)
    {
        ffp[i] = fp[i];
    }
    fprec = *precision;
    ret   = xdr3dfcoord(xdrs, ffp, size, &fprec, magic_number);

    *precision = fprec;
    for (i = 0; (i < isize); i++)
    {
        fp[i] = ffp[i];
    }

    sfree(ffp);
    return ret;
#else
    return xdr3dfcoord(
            xdrs, reinterpret_cast<float*>(fp), size, reinterpret_cast<float*>(precision), magic_number);
#endif
}

int xdr_int32(XDR* xdrs, int32_t* i)
{
    // Note that this implementation assumes that an int is at least
    // 32 bits, which is not strictly required by the language, but
    // good enough in practice on 32- or 64-bit systems. GROMACS
    // requires 64-bit systems.
    static_assert(sizeof(int) >= 4, "XDR handling assumes that an int32_t can be stored in an int");
    int temporary = static_cast<int>(*i);
    int ret       = xdr_int(xdrs, &temporary);
    *i            = static_cast<int32_t>(temporary);

    return ret;
}

int xdr_int64(XDR* xdrs, int64_t* i)
{
    // Note that this implementation assumes that an int is at least
    // 32 bits, which is not strictly required by the language, but
    // good enough in practice on 32- or 64-bit systems. GROMACS
    // requires 64-bit systems.
    static_assert(2 * sizeof(int) >= 8,
                  "XDR handling assumes that an int64_t can be stored in two ints");

    static const uint64_t two_p32_m1 = 0xFFFFFFFF;

    uint64_t imaj64 = ((*i) >> 32) & two_p32_m1;
    uint64_t imin64 = (*i) & two_p32_m1;
    int      imaj   = static_cast<int>(imaj64);
    int      imin   = static_cast<int>(imin64);
    int      ret    = xdr_int(xdrs, &imaj);
    ret |= xdr_int(xdrs, &imin);

    *i = ((static_cast<uint64_t>(imaj) << 32) | (static_cast<uint64_t>(imin) & two_p32_m1));

    return ret;
}
