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
/*! \libinternal \file
 * \brief
 * Provides templated functions to map C++ types to MPI datatypes.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_MPITYPES_H
#define GMX_UTILITY_MPITYPES_H

#include "config.h"

#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

#if GMX_MPI

//! Returns the MPI data type corresponding to int
template<typename T>
std::enable_if_t<std::is_same_v<T, int>, MPI_Datatype> mpiType()
{
    return MPI_INT;
}

//! Returns the MPI data type corresponding to long
template<typename T>
std::enable_if_t<std::is_same_v<T, long>, MPI_Datatype> mpiType()
{
    return MPI_LONG;
}

//! Returns the MPI data type corresponding to int64_t
template<typename T>
std::enable_if_t<std::is_same_v<T, int64_t> && !std::is_same_v<int64_t, long>, MPI_Datatype> mpiType()
{
    return MPI_INT64_T;
}

//! Returns the MPI data type corresponding to float
template<typename T>
std::enable_if_t<std::is_same_v<T, float>, MPI_Datatype> mpiType()
{
    return MPI_FLOAT;
}

//! Returns the MPI data type corresponding to double
template<typename T>
std::enable_if_t<std::is_same_v<T, double>, MPI_Datatype> mpiType()
{
    return MPI_DOUBLE;
}

#endif // GMX_MPI

} // namespace gmx

#endif
