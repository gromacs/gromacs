/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *  \brief Implements PME GPU spline calculation and charge spreading in SYCL.
 *  TODO: consider always pre-sorting particles (as in DD case).
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

#include "gromacs/gpu_utils/syclutils.h"

#include "pme_gpu_types_host.h"
#include "pme_grid.h"

struct PmeGpuGridParams;
struct PmeGpuAtomParams;
struct PmeGpuDynamicParams;

struct PmeGpuPipeliningParams
{
    int  pipelineAtomStart;
    int  pipelineAtomEnd;
    bool usePipeline;
};

template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
class PmeSplineAndSpreadKernel : public ISyclKernelFunctor
{
public:
    PmeSplineAndSpreadKernel();
    void setArg(size_t argIndex, void* arg) override;
    void launch(const KernelLaunchConfig& config, const DeviceStream& deviceStream) override;

private:
    PmeGpuGridParams*      gridParams_;
    PmeGpuAtomParams*      atomParams_;
    PmeGpuDynamicParams*   dynamicParams_;
    PmeGpuPipeliningParams pipeliningParams_;

    void reset();
};
