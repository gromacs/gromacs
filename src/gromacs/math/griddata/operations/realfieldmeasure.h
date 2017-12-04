/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#ifndef GMX_MATH_REALFIELDMEASURE_H
#define GMX_MATH_REALFIELDMEASURE_H

#include <vector>
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/vec.h"

namespace gmx
{

namespace realgriddatameasure
{
RVec center_of_mass(const GridDataReal3D &realfield);
}

namespace datavectormeasure
{
/*! \brief The sum of all grid values. */
real sum(const std::vector<real> &vecData_);
/*! \brief The minimum grid data value. */
real min(const std::vector<real> &vecData_);
/*! \brief The maximum grid data value. */
real max(const std::vector<real> &vecData_);
/*! \brief The mean grid data value. */
real mean(const std::vector<real> &vecData_);
/*! \brief
 * The size of the griddata in bytes.
 */
size_t data_size(const std::vector<real> &vecData_);

real sumOfSquareDeviation(const std::vector<real> &vecData_);

real normSquared(const std::vector<real> &vecData_);

real norm(const std::vector<real> &vecData_);
/*! \brief The variance of data values.
 *
 * Two pass algorithm, where mean is calculated fist, due to large expected
 * amount of data. (Typical data size=1,000,000)
 */
real var(const std::vector<real> &vecData_);
/*! \brief The root mean square deviation of grid data values. */
real rms(const std::vector<real> &vecData_);

real entropy(const std::vector<real> &realfield, real unitCellVolume);
}

}      // namespace gmx
#endif
