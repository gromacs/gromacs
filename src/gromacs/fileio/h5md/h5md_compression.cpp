/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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

/*! \brief Definitions of H5md compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_compression.h"

#if GMX_HAVE_SZ3
#    include <H5Z_SZ3.hpp>
#endif

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

void setSz3Compression(const hid_t creationPropertyList, const double absoluteError)
{
#if GMX_HAVE_SZ3
    GMX_H5MD_THROW_UPON_ERROR(absoluteError <= 0.0,
                              "Absolute error for SZ3 compression must be positive");
    GMX_H5MD_THROW_UPON_ERROR(
            H5Pset_sz3(creationPropertyList, H5Z_SZ3_ALGO_BIOMD, H5Z_SZ3_EB_ABS, absoluteError, 0, 0, 0) < 0,
            "Cannot set SZ3 compression filter");
#else
    GMX_THROW(NotImplementedError("GROMACS was compiled without SZ3 support"));
    GMX_UNUSED_VALUE(creationPropertyList);
    GMX_UNUSED_VALUE(absoluteError);
#endif
}

} // namespace gmx
