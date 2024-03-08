/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 *
 * \brief Implements backend-specific functions of the update-constraints in SYCL.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/utility/gmxassert.h"

#include "update_constrain_gpu_internal.h"

//! \brief Class name for scaling kernel
class ScaleKernel;

namespace gmx
{

//! \brief Function returning the scaling kernel lambda.
static auto scaleKernel(Float3* gm_x, const ScalingMatrix& scalingMatrix)
{
    return [=](sycl::id<1> itemIdx) {
        Float3 x      = gm_x[itemIdx];
        x[0]          = scalingMatrix.xx * x[0] + scalingMatrix.yx * x[1] + scalingMatrix.zx * x[2];
        x[1]          = scalingMatrix.yy * x[1] + scalingMatrix.zy * x[2];
        x[2]          = scalingMatrix.zz * x[2];
        gm_x[itemIdx] = x;
    };
}

void launchScaleCoordinatesKernel(const int            numAtoms,
                                  DeviceBuffer<Float3> d_coordinates,
                                  const ScalingMatrix& mu,
                                  const DeviceStream&  deviceStream)
{
    const sycl::range<1> rangeAllAtoms(numAtoms);
    sycl::queue          queue = deviceStream.stream();

    sycl::event e = queue.submit([&](sycl::handler& cgh) {
        auto kernel = scaleKernel(d_coordinates.get_pointer(), mu);
        cgh.parallel_for<ScaleKernel>(rangeAllAtoms, kernel);
    });
    // TODO: Although this only happens on the pressure coupling steps, this synchronization
    //       can affect the performance if nstpcouple is small. See Issue #4018
    e.wait_and_throw();
}

} // namespace gmx
