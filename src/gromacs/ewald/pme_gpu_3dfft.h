/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020, by the GROMACS development team.
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
 *  \brief Declares the 3D FFT class for PME.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_3DFFT_H
#define GMX_EWALD_PME_GPU_3DFFT_H

#include <memory>

#include "gromacs/fft/fft.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"

class DeviceContext;
class DeviceStream;
struct PmeGpu;

/*! \internal \brief
 * A 3D FFT class for performing R2C/C2R transforms
 * \todo Make this class actually parallel over multiple GPUs
 */
class GpuParallel3dFft
{
public:
    /*! \brief
     * Constructs GPU FFT plans for performing 3D FFT on a PME grid.
     *
     * \param[in]  realGridSize           Dimensions of the real grid
     * \param[in]  realGridSizePadded     Dimensions of the real grid with padding
     * \param[in]  complexGridSizePadded  Dimensions of the real grid with padding
     * \param[in]  useDecomposition       Whether PME decomposition will be used
     * \param[in]  performOutOfPlaceFFT   Whether the FFT will be performed out-of-place
     * \param[in]  context                GPU context.
     * \param[in]  pmeStream              GPU stream for PME.
     * \param[in]  realGrid               Device buffer of floats for the real grid
     * \param[in]  complexGrid            Device buffer of complex floats for the complex grid
     */
    GpuParallel3dFft(ivec                 realGridSize,
                     ivec                 realGridSizePadded,
                     ivec                 complexGridSizePadded,
                     bool                 useDecomposition,
                     bool                 performOutOfPlaceFFT,
                     const DeviceContext& context,
                     const DeviceStream&  pmeStream,
                     DeviceBuffer<float>  realGrid,
                     DeviceBuffer<float>  complexGrid);

    /*! \brief Destroys the FFT plans. */
    ~GpuParallel3dFft();
    /*! \brief Performs the FFT transform in given direction
     *
     * \param[in]  dir           FFT transform direction specifier
     * \param[out] timingEvent   pointer to the timing event where timing data is recorded
     */
    void perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

#endif
