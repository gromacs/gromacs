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
 *  \brief Implements stub GPU 3D FFT routines for CPU-only builds
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Gaurav Garg <gaugarg@nvidia.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_impl.h"

#include "gromacs/gpu_utils/devicebuffer.h"

namespace gmx
{

Gpu3dFft::Impl::Impl() : complexGrid_(nullptr) {}

Gpu3dFft::Impl::Impl(bool performOutOfPlaceFFT) :
    performOutOfPlaceFFT_(performOutOfPlaceFFT), complexGrid_(nullptr)
{
}

void Gpu3dFft::Impl::allocateComplexGrid(const ivec           complexGridSizePadded,
                                         DeviceBuffer<float>* realGrid,
                                         DeviceBuffer<float>* complexGrid,
                                         const DeviceContext& context)
{
    if (performOutOfPlaceFFT_)
    {
        const int newComplexGridSize =
                complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ] * 2;

        allocateDeviceBuffer(complexGrid, newComplexGridSize, context);
    }
    else
    {
        *complexGrid = *realGrid;
    }

    complexGrid_ = *complexGrid;
}

void Gpu3dFft::Impl::deallocateComplexGrid()
{
    if (performOutOfPlaceFFT_)
    {
        freeDeviceBuffer(&complexGrid_);
    }
}

Gpu3dFft::Impl::~Impl() = default;

} // namespace gmx
