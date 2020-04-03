/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief
 * Declares ClfftInitializer class, which initializes and
 * tears down the clFFT library resources in OpenCL builds,
 * and does nothing in other builds, and a factory function for it.
 *
 * clFFT itself is used in the OpenCL implementation of PME
 * for 3D R2C/C2R transforms. It is know to work with NVidia
 * OpenCL, AMD fglrx and AMDGPU-PRO drivers, and to not work with
 * AMD Rocm dev preview as of May 2018 (#2515).
 * TODO: find out compatibility with Intel once the rest of PME
 * gets there (#2516), or by building clFFT own tests.
 *
 * \inlibraryapi
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */
#ifndef GMX_GPU_UTILS_CLFFTINITIALIZER_H
#define GMX_GPU_UTILS_CLFFTINITIALIZER_H

#include <memory>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \libinternal
 * \brief Handle clFFT library init and tear down in RAII style also
 * with mutual exclusion.
 *
 * Only one thread per process needs to attempt to set up and tear
 * down the clFFT library, but this wrapper object ensures that it is
 * safe to do so more than once, or from any thread. It is the
 * responsibility of the caller to ensure that no use of the clFFT
 * library is made before making an object of this type, or after the
 * first such object is destroyed. If no more objects remain to be
 * destroyed, then it is safe to create another and resume clFFT work.
 *
 * \todo Consider making a composite object that also handles
 * on-demand compilation, managing lifetime of PME FFT kernel programs
 * to avoid exhausting resources and/or recompiling kernels previously
 * used.  See Issue #2535.
 *
 * \todo Consider implementing an appropriate flavor of the nifty
 * counter idiom so that both static initialization and
 * deinitialization can work in a fast, leak-free, and thread-safe way
 * without imposing constraints on the calling code.
 * See Issue #2535. */
class ClfftInitializer
{
public:
    /*! \brief Constructor
     *
     * This initializers the clFFT library if there is a
     * possibility of an FFT task on the device, and preserves it
     * until destruction. Any time required for this
     * initialization or tear down should not be accrued to
     * per-MD-step counters. */
    ClfftInitializer();
    //! Destructor
    ~ClfftInitializer();

    GMX_DISALLOW_COPY_AND_ASSIGN(ClfftInitializer);
};

} // namespace gmx

#endif
