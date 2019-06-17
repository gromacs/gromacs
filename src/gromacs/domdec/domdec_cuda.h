/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declaration of GPU halo exchange.
 *
 * \author Alan Gray <alang@nvidia.com>
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_DD_GPU_H
#define GMX_DD_GPU_H

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class DomdecCuda
{

    public:
        /*! \brief Creates Domdec GPU object
         *
         * \param [inout] dd      domdec structure
         */
        DomdecCuda(gmx_domdec_t gmx_unused *dd);
        ~DomdecCuda();

        /*! \brief
         * Initialization for GPU halo exchange of position buffer
         * \param [inout] dd      domdec structure, to be populated with index and param info for exchange
         * \param [in] box        box matrix required for shift information
         * \param [in] d_x_ptr   pointer to position buffer in GPU memory
         * \param [in] stream_ptr   pointer to CUDA stream to be used for init operations
         */
        void initHaloExchange(gmx_domdec_t gmx_unused *dd, matrix gmx_unused box,
                              rvec gmx_unused *d_x_ptr, void gmx_unused *stream_ptr);


        /*! \brief
         * GPU halo exchange of position buffer
         * \param [inout] d_x_ptr   pointer to position buffer in GPU memory
         * \param [in] streamNonLocal       CUDA stream to be used for buffer packing operation
         */
        void applyXHaloExchange(rvec gmx_unused *d_x_ptr, void gmx_unused *streamNonLocal);

        /*! \brief
         * GPU halo exchange of force buffer
         * \param [inout] d_f_ptr   pointer to force buffer in GPU memory
         * \param [inout] fshift    force shift array
         * \param [in] stream       CUDA stream to be used for buffer packing operation
         */
        void applyFHaloExchange(rvec gmx_unused *d_f_ptr, rvec gmx_unused *fshift, void gmx_unused *stream);


    private:
        class Impl;
        gmx::PrivateImplPointer<Impl> impl_;

};

} //namespace gmx

#endif
