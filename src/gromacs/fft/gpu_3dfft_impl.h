/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2021, by the GROMACS development team, led by
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
 *  \brief Declares the GPU 3D FFT routines.
 *
 *  \author Gaurav Garg <gaugarg@nvidia.com>
 *  \ingroup module_fft
 */

#ifndef GMX_FFT_GPU_3DFFT_IMPL_H
#define GMX_FFT_GPU_3DFFT_IMPL_H

#include "gromacs/fft/fft.h"
#include "gromacs/fft/gpu_3dfft.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"


namespace gmx
{
/*! \internal \brief
 * Impl base class for all FFT backends
 */
class Gpu3dFft::Impl
{
public:
    //! Default constructor
    Impl() = default;

    /*! \brief
     * Constructs GPU FFT plans for performing 3D FFT on a PME grid.
     *
     * \param[in]  allocateGrids                True if fft grids are to be allocated, false if pre-allocated
     * \param[in]  comm                         MPI communicator, used with distributed-FFT backends
     * \param[in]  gridSizesInXForEachRank      Number of grid points used with each rank in X-dimension
     * \param[in]  gridSizesInYForEachRank      Number of grid points used with each rank in Y-dimension
     * \param[in]  nz                           Grid dimension in Z
     * \param[in]  performOutOfPlaceFFT         Whether the FFT will be performed out-of-place
     * \param[in]  context                      GPU context.
     * \param[in]  pmeStream                    GPU stream for PME.
     * \param[in,out]  realGridSize             Dimensions of the local real grid, out if allocateGrids=true
     * \param[in,out]  realGridSizePadded       Dimensions of the local real grid with padding, out if allocateGrids=true
     * \param[in,out]  complexGridSizePadded    Dimensions of the local complex grid with padding, out if allocateGrids=true
     * \param[in,out]  realGrid                 Device buffer of floats for the local real grid, out if allocateGrids=true
     * \param[in,out]  complexGrid              Device buffer of complex floats for the local complex grid, out if allocateGrids=true
     */
    Impl(bool                 allocateGrids,
         MPI_Comm             comm,
         ArrayRef<const int>  gridSizesInXForEachRank,
         ArrayRef<const int>  gridSizesInYForEachRank,
         int                  nz,
         bool                 performOutOfPlaceFFT,
         const DeviceContext& context,
         const DeviceStream&  pmeStream,
         ivec                 realGridSize,
         ivec                 realGridSizePadded,
         ivec                 complexGridSizePadded,
         DeviceBuffer<float>* realGrid,
         DeviceBuffer<float>* complexGrid);

    /*! \brief Default destructor */
    virtual ~Impl() = default;

    //! \copydoc Gpu3dFft::perform3dFft
    virtual void perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent) = 0;
};

} // namespace gmx

#endif
