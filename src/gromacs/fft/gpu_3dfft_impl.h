/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
    Impl();

    //! Constructor called from dervied classes
    Impl(bool performOutOfPlaceFFT);

    /*! \brief
     * Constructs GPU FFT plans for performing 3D FFT on a PME grid.
     *
     * \param[in]  allocateRealGrid          True if fft real-grid is to be allocated,
     *                                               false if pre-allocated
     * \param[in]  comm                      MPI communicator, used with distributed-FFT backends
     * \param[in]  gridSizesInXForEachRank   Number of grid points used with each rank in X-dimension
     * \param[in]  gridSizesInYForEachRank   Number of grid points used with each rank in Y-dimension
     * \param[in]  nz                        Grid dimension in Z
     * \param[in]  performOutOfPlaceFFT      Whether the FFT will be performed out-of-place
     * \param[in]  context                   GPU context
     * \param[in]  pmeStream                 GPU stream for PME
     * \param[in,out]  realGridSize          Dimensions of the local real grid, out if allocateRealGrid=true
     * \param[in,out]  realGridSizePadded    Dimensions of the local real grid with padding, out if allocateRealGrid=true
     * \param[in,out]  complexGridSizePadded Dimensions of the local complex grid with padding, out if allocateRealGrid=true
     * \param[in,out]  realGrid              Device buffer of floats for the local real grid, out if allocateRealGrid=true
     * \param[out]  complexGrid              Device buffer of complex floats for the local complex grid
     */
    Impl(bool                 allocateRealGrid,
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
    virtual ~Impl();

    //! \copydoc Gpu3dFft::perform3dFft
    virtual void perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent) = 0;

protected:
    //! Allocate and assign complexGrid
    void allocateComplexGrid(const ivec           complexGridSizePadded,
                             DeviceBuffer<float>* realGrid,
                             DeviceBuffer<float>* complexGrid,
                             const DeviceContext& context);

    //! free complexGrid
    void deallocateComplexGrid();

    /*! \brief A boolean which tells whether the complex and real grids are different or same. Currenty true. */
    bool performOutOfPlaceFFT_ = false;
    /*! \brief FFT complex grid */
    DeviceBuffer<float> complexGrid_;
};

} // namespace gmx

#endif
