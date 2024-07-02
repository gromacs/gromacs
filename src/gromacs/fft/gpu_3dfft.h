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
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Gaurav Garg <gaugarg@nvidia.com>
 *  \ingroup module_fft
 */

#ifndef GMX_FFT_GPU_3DFFT_H
#define GMX_FFT_GPU_3DFFT_H

#include <memory>

#include "gromacs/fft/fft.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/utility/gmxmpi.h"

class DeviceContext;
class DeviceStream;

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \internal \brief
 * Enum specifying all GPU FFT backends supported by GROMACS
 * Some of the backends support only single GPU, some only multi-node, multi-GPU
 */
enum class FftBackend
{
    Cufft,              //!< supports only single-GPU
    OclVkfft,           //!< supports only single-GPU
    Ocl,                //!< supports only single-GPU
    CuFFTMp,            //!< supports only multi-gpu
    HeFFTe_CUDA,        //!< supports only multi-gpu
    HeFFTe_Sycl_OneMkl, //!< supports only multi-gpu
    HeFFTe_Sycl_Rocfft, //!< supports only multi-gpu
    HeFFTe_Sycl_cuFFT,  //!< supports only multi-gpu
    SyclMkl,            //!< supports only single-GPU
    SyclOneMkl,         //!< supports only single-GPU
    SyclRocfft,         //!< supports only single-GPU
    SyclVkfft,          //!< supports only single-GPU
    SyclBbfft,          //!< supports only single-GPU
    Sycl,               //!< stubs for not supported configurations
    Hipfft,             //!< supports only single-GPU
    HipVkfft,           //!< supports only single-GPU
    Count
};

/*! \internal \brief
 * A 3D FFT class for performing R2C/C2R transforms
 */
class Gpu3dFft
{
public:
    /*! \brief
     * Construct 3D FFT object for given backend
     *
     * \param[in]  backend                   FFT backend to be instantiated
     * \param[in]  allocateRealGrid          True if fft real-grid is to be allocated,
     *                                          false if pre-allocated
     * \param[in]  comm                      MPI communicator, used with distributed-FFT backends
     * \param[in]  gridSizesInXForEachRank   Number of grid points used with each rank in X-dimension
     * \param[in]  gridSizesInYForEachRank   Number of grid points used with each rank in Y-dimension
     * \param[in]  nz                        Grid dimension in Z
     * \param[in]  performOutOfPlaceFFT      Whether the FFT will be performed out-of-place
     * \param[in] context                    GPU context.
     * \param[in]  pmeStream                 GPU stream for PME.
     * \param[in,out]  realGridSize          Dimensions of the local real grid, out if allocateRealGrid=true
     * \param[in,out]  realGridSizePadded    Dimensions of the local real grid with padding, out if allocateRealGrid=true
     * \param[in,out]  complexGridSizePadded Dimensions of the local complex grid with padding, out if allocateRealGrid=true
     * \param[in,out]  realGrid              Device buffer of floats for the local real grid, out if allocateRealGrid=true
     * \param[out] complexGrid               Device buffer of complex floats for the local complex grid
     */
    Gpu3dFft(FftBackend           backend,
             bool                 allocateRealGrid,
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

    /*! \brief Destroys the FFT plans. */
    ~Gpu3dFft();
    /*! \brief Performs the FFT transform in given direction
     *
     * \param[in]  dir           FFT transform direction specifier
     * \param[out] timingEvent   pointer to the timing event where timing data is recorded
     */
    void perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent);

private:
    class Impl;
    class ImplCuFft;
    class ImplCuFftMp;
    class ImplOclVkfft;
    class ImplOcl;
    class ImplSyclMkl;
    // oneMKL has an identical interface to MKL, so uses the ImplSyclMkl implemenation.
    class ImplSyclBbfft;
    class ImplSyclRocfft;
    class ImplSyclVkfft;
    class ImplSycl;
    class ImplHipFft;
    class ImplHipVkFft;

    template<typename backend_tag>
    class ImplHeFfte;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
