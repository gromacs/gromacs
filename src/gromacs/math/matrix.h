/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal
 * \file
 * \brief Declares special case of 3x3 matrix frequently used.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_MATRIX_H_
#define GMX_MATH_MATRIX_H_

#include <algorithm>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "multidimarray.h"

namespace gmx
{

//! Three-by-three matrix
using Matrix3x3 = MultiDimArray < std::array<real, DIM*DIM>, extents < DIM, DIM>>;

static inline void clear_mat(Matrix3x3::view_type m)
{
    std::fill(begin(m), end(m), 0);
}

static inline void mmul(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    dest(XX, XX) = a(XX, XX) * b(XX, XX) + a(XX, YY) * b(YY, XX) + a(XX, ZZ) * b(ZZ, XX);
    dest(YY, XX) = a(YY, XX) * b(XX, XX) + a(YY, YY) * b(YY, XX) + a(YY, ZZ) * b(ZZ, XX);
    dest(ZZ, XX) = a(ZZ, XX) * b(XX, XX) + a(ZZ, YY) * b(YY, XX) + a(ZZ, ZZ) * b(ZZ, XX);
    dest(XX, YY) = a(XX, XX) * b(XX, YY) + a(XX, YY) * b(YY, YY) + a(XX, ZZ) * b(ZZ, YY);
    dest(YY, YY) = a(YY, XX) * b(XX, YY) + a(YY, YY) * b(YY, YY) + a(YY, ZZ) * b(ZZ, YY);
    dest(ZZ, YY) = a(ZZ, XX) * b(XX, YY) + a(ZZ, YY) * b(YY, YY) + a(ZZ, ZZ) * b(ZZ, YY);
    dest(XX, ZZ) = a(XX, XX) * b(XX, ZZ) + a(XX, YY) * b(YY, ZZ) + a(XX, ZZ) * b(ZZ, ZZ);
    dest(YY, ZZ) = a(YY, XX) * b(XX, ZZ) + a(YY, YY) * b(YY, ZZ) + a(YY, ZZ) * b(ZZ, ZZ);
    dest(ZZ, ZZ) = a(ZZ, XX) * b(XX, ZZ) + a(ZZ, YY) * b(YY, ZZ) + a(ZZ, ZZ) * b(ZZ, ZZ);
}

static inline void transpose(Matrix3x3::const_view_type src, Matrix3x3::view_type dest)
{
    dest(XX, XX) = src(XX, XX);
    dest(YY, XX) = src(XX, YY);
    dest(ZZ, XX) = src(XX, ZZ);
    dest(XX, YY) = src(YY, XX);
    dest(YY, YY) = src(YY, YY);
    dest(ZZ, YY) = src(YY, ZZ);
    dest(XX, ZZ) = src(ZZ, XX);
    dest(YY, ZZ) = src(ZZ, YY);
    dest(ZZ, ZZ) = src(ZZ, ZZ);
}

static inline void tmmul(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    /* Computes dest=mmul(transpose(a),b,dest) - used in do_pr_pcoupl */
    dest(XX, XX) = a(XX, XX) * b(XX, XX) + a(YY, XX) * b(YY, XX) + a(ZZ, XX) * b(ZZ, XX);
    dest(XX, YY) = a(XX, XX) * b(XX, YY) + a(YY, XX) * b(YY, YY) + a(ZZ, XX) * b(ZZ, YY);
    dest(XX, ZZ) = a(XX, XX) * b(XX, ZZ) + a(YY, XX) * b(YY, ZZ) + a(ZZ, XX) * b(ZZ, ZZ);
    dest(YY, XX) = a(XX, YY) * b(XX, XX) + a(YY, YY) * b(YY, XX) + a(ZZ, YY) * b(ZZ, XX);
    dest(YY, YY) = a(XX, YY) * b(XX, YY) + a(YY, YY) * b(YY, YY) + a(ZZ, YY) * b(ZZ, YY);
    dest(YY, ZZ) = a(XX, YY) * b(XX, ZZ) + a(YY, YY) * b(YY, ZZ) + a(ZZ, YY) * b(ZZ, ZZ);
    dest(ZZ, XX) = a(XX, ZZ) * b(XX, XX) + a(YY, ZZ) * b(YY, XX) + a(ZZ, ZZ) * b(ZZ, XX);
    dest(ZZ, YY) = a(XX, ZZ) * b(XX, YY) + a(YY, ZZ) * b(YY, YY) + a(ZZ, ZZ) * b(ZZ, YY);
    dest(ZZ, ZZ) = a(XX, ZZ) * b(XX, ZZ) + a(YY, ZZ) * b(YY, ZZ) + a(ZZ, ZZ) * b(ZZ, ZZ);
}

static inline void mtmul(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    /* Computes dest=mmul(a,transpose(b),dest) - used in do_pr_pcoupl */
    dest(XX, XX) = a(XX, XX) * b(XX, XX) + a(XX, YY) * b(XX, YY) + a(XX, ZZ) * b(XX, ZZ);
    dest(XX, YY) = a(XX, XX) * b(YY, XX) + a(XX, YY) * b(YY, YY) + a(XX, ZZ) * b(YY, ZZ);
    dest(XX, ZZ) = a(XX, XX) * b(ZZ, XX) + a(XX, YY) * b(ZZ, YY) + a(XX, ZZ) * b(ZZ, ZZ);
    dest(YY, XX) = a(YY, XX) * b(XX, XX) + a(YY, YY) * b(XX, YY) + a(YY, ZZ) * b(XX, ZZ);
    dest(YY, YY) = a(YY, XX) * b(YY, XX) + a(YY, YY) * b(YY, YY) + a(YY, ZZ) * b(YY, ZZ);
    dest(YY, ZZ) = a(YY, XX) * b(ZZ, XX) + a(YY, YY) * b(ZZ, YY) + a(YY, ZZ) * b(ZZ, ZZ);
    dest(ZZ, XX) = a(ZZ, XX) * b(XX, XX) + a(ZZ, YY) * b(XX, YY) + a(ZZ, ZZ) * b(XX, ZZ);
    dest(ZZ, YY) = a(ZZ, XX) * b(YY, XX) + a(ZZ, YY) * b(YY, YY) + a(ZZ, ZZ) * b(YY, ZZ);
    dest(ZZ, ZZ) = a(ZZ, XX) * b(ZZ, XX) + a(ZZ, YY) * b(ZZ, YY) + a(ZZ, ZZ) * b(ZZ, ZZ);
}

static inline real det(Matrix3x3::const_view_type a)
{
    return (    a(XX, XX) * (a(YY, YY) * a(ZZ, ZZ) - a(ZZ, YY) * a(YY, ZZ))
                - a(YY, XX) * (a(XX, YY) * a(ZZ, ZZ) - a(ZZ, YY) * a(XX, ZZ))
                + a(ZZ, XX) * (a(XX, YY) * a(YY, ZZ) - a(YY, YY) * a(XX, ZZ)));
}

static inline void m_add(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    std::transform(begin(a), end(a), begin(b), begin(dest), std::plus<>());
}

static inline void m_sub(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    std::transform(begin(a), end(a), begin(b), begin(dest), std::minus<>());
}

static inline void msmul(Matrix3x3::const_view_type m1, real scalar, Matrix3x3::view_type dest)
{
    std::transform(begin(m1), end(m1), begin(dest), [scalar](real matrixEntry) { return matrixEntry * scalar; });
}

static inline void mvmul(Matrix3x3::const_view_type a, const rvec src, rvec dest)
{
    dest[XX] = a(XX, XX) * src[XX] + a(XX, YY) * src[YY] + a(XX, ZZ) * src[ZZ];
    dest[YY] = a(YY, XX) * src[XX] + a(YY, YY) * src[YY] + a(YY, ZZ) * src[ZZ];
    dest[ZZ] = a(ZZ, XX) * src[XX] + a(ZZ, YY) * src[YY] + a(ZZ, ZZ) * src[ZZ];
}

static inline void mmul_ur0(Matrix3x3::const_view_type a, Matrix3x3::const_view_type b, Matrix3x3::view_type dest)
{
    dest(XX, XX) = a(XX, XX) * b(XX, XX);
    dest(XX, YY) = 0.0;
    dest(XX, ZZ) = 0.0;
    dest(YY, XX) = a(YY, XX) * b(XX, XX) + a(YY, YY) * b(YY, XX);
    dest(YY, YY) = a(YY, YY) * b(YY, YY);
    dest(YY, ZZ) = 0.0;
    dest(ZZ, XX) = a(ZZ, XX) * b(XX, XX) + a(ZZ, YY) * b(YY, XX) + a(ZZ, ZZ) * b(ZZ, XX);
    dest(ZZ, YY) = a(ZZ, YY) * b(YY, YY) + a(ZZ, ZZ) * b(ZZ, YY);
    dest(ZZ, ZZ) = a(ZZ, ZZ) * b(ZZ, ZZ);
}

static inline void mvmul_ur0(Matrix3x3::const_view_type a, const rvec src, rvec dest)
{
    dest[ZZ] = a(ZZ, XX) * src[XX] + a(ZZ, YY) * src[YY] + a(ZZ, ZZ) * src[ZZ];
    dest[YY] = a(YY, XX) * src[XX] + a(YY, YY) * src[YY];
    dest[XX] = a(XX, XX) * src[XX];
}

static inline void tmvmul_ur0(Matrix3x3::const_view_type a, const rvec src, rvec dest)
{
    dest[XX] = a(XX, XX) * src[XX] + a(YY, XX) * src[YY] + a(ZZ, XX) * src[ZZ];
    dest[YY] = a(YY, YY) * src[YY] + a(ZZ, YY) * src[ZZ];
    dest[ZZ] = a(ZZ, ZZ) * src[ZZ];
}

} // namespace gmx

#endif
