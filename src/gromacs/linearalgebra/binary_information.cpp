/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
/*! \internal \file
 * \brief Implements functionality for printing information about the
 * linear algebra support in the currently running binary
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/linearalgebra/binary_information.h"

#if GMX_FFT_MKL
#    include <mkl.h>
#endif

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

namespace
{

// This function is duplicated in fft, keep both in sync.
std::string describeMkl()
{
#if GMX_FFT_MKL
    MKLVersion mklVersion;
    mkl_get_version(&mklVersion);
    auto description = gmx::formatString("Intel MKL version %d.%d.%d Build %s",
                                         mklVersion.MajorVersion,
                                         mklVersion.MinorVersion,
                                         mklVersion.UpdateVersion,
                                         mklVersion.Build);
    if (mklVersion.ProductStatus != std::string("Product"))
    {
        description += " ";
        description += mklVersion.ProductStatus;
    }
    return description;
#else
    return "Intel MKL";
#endif
}

// Turn A into a string literal without expanding macro definitions.
#define STRINGIZE_NO_EXPAND(A) #A

// Turn A into a string literal after macro-expanding it.
#define STRINGIZE(A) STRINGIZE_NO_EXPAND(A)

const char* sc_blasDescription   = STRINGIZE(GMX_DESCRIBE_BLAS);
const char* sc_lapackDescription = STRINGIZE(GMX_DESCRIBE_LAPACK);

} // namespace

namespace gmx
{

std::string blasDescription()
{
    // Describe the BLAS library. We generally don't know
    // much about what external library was detected, but we do in the
    // case of MKL so then it is reported.
    const bool descriptionContainsMkl = std::strstr(sc_blasDescription, "MKL") != nullptr;
    if (descriptionContainsMkl)
    {
        return describeMkl();
    }
    return sc_blasDescription;
}

std::string lapackDescription()
{
    // Describe the LAPACK library. We generally don't know
    // much about what external library was detected, but we do in the
    // case of MKL so then it is reported.
    const bool descriptionContainsMkl = std::strstr(sc_lapackDescription, "MKL") != nullptr;
    if (descriptionContainsMkl)
    {
        return describeMkl();
    }
    return sc_lapackDescription;
}

} // namespace gmx
