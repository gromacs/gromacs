/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *  \brief Implements PME GPU spline calculation and charge spreading in SYCL.
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/syclutils.h"
#include "gromacs/math/vectypes.h"

#include "pme_gpu_internal.h"
#include "pme_gpu_types.h"

struct PmeGpuConstParams;
struct PmeGpuGridParams;

//! Contains most of the parameters used by the solve kernel
struct SolveKernelParams
{
    /*! \brief Ewald solving factor = (M_PI / pme->ewaldcoeff_q)^2 */
    float ewaldFactor;
    /*! \brief Real-space grid data dimensions. */
    gmx::IVec realGridSize;
    /*! \brief Fourier grid dimensions. This counts the complex numbers! */
    gmx::IVec complexGridSize;
    /*! \brief Fourier grid dimensions (padded). This counts the complex numbers! */
    gmx::IVec complexGridSizePadded;
    /*! \brief Offsets for X/Y/Z components of d_splineModuli */
    gmx::IVec splineValuesOffset;
    /*! \brief Reciprocal (inverted unit cell) box. */
    gmx::RVec recipBox[DIM];
    /*! \brief The unit cell volume for solving. */
    float boxVolume;
    /*! \brief Electrostatics coefficient = c_one4PiEps0 / pme->epsilon_r */
    float elFactor;
};

//! The kernel for PME solve
template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int gridIndex, int subGroupSize>
class PmeSolveKernel : public ISyclKernelFunctor
{
public:
    PmeSolveKernel();
    //! Sets the kernel arguments
    void setArg(size_t argIndex, void* arg) override;
    //! Launches the kernel with given \c config and \c deviceStream
    cl::sycl::event launch(const KernelLaunchConfig& config, const DeviceStream& deviceStream) override;

private:
    //! Kernel argument set by \c setArg()
    PmeGpuConstParams* constParams_ = nullptr;
    //! Kernel argument set by \c setArg()
    PmeGpuGridParams* gridParams_ = nullptr;
    //! Kernel argument set by \c setArg()
    SolveKernelParams solveKernelParams_;

    //! Called after each launch to ensure we set the arguments again properly
    void reset();
};
