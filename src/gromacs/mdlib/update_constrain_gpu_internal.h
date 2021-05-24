/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Declares GPU implementations of backend-specific update-constraints functions.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_UPDATE_CONSTRAIN_GPU_INTERNAL_H
#define GMX_MDLIB_UPDATE_CONSTRAIN_GPU_INTERNAL_H

#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits.h"

class GpuEventSynchronizer;

namespace gmx
{

/*! \internal \brief Scaling matrix struct.
 *
 * \todo Should be generalized.
 */
struct ScalingMatrix
{
    ScalingMatrix(const matrix m) :
        xx(m[XX][XX]), yy(m[YY][YY]), zz(m[ZZ][ZZ]), yx(m[YY][XX]), zx(m[ZZ][XX]), zy(m[ZZ][YY])
    {
    }
    float xx, yy, zz, yx, zx, zy;
};

/*! \brief Launches positions of velocities scaling kernel.
 *
 * \param[in] numAtoms       Number of atoms in the system.
 * \param[in] d_coordinates  Device buffer with position or velocities to be scaled.
 * \param[in] mu             Scaling matrix.
 * \param[in] deviceStream   Stream to launch kernel in.
 */
void launchScaleCoordinatesKernel(int                  numAtoms,
                                  DeviceBuffer<Float3> d_coordinates,
                                  const ScalingMatrix& mu,
                                  const DeviceStream&  deviceStream);

} // namespace gmx

#endif // GMX_MDLIB_UPDATE_CONSTRAIN_GPU_INTERNAL_H
