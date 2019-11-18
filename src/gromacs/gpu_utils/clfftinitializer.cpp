/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief
 * Implements ClfftInitializer class.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "clfftinitializer.h"

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/mutex.h"
#include "gromacs/utility/stringutil.h"

#if GMX_GPU == GMX_GPU_OPENCL
#    include <clFFT.h>
#endif

namespace gmx
{

namespace
{

/*! \brief The clFFT library may only be initialized once per process,
 * and this is orchestrated by this shared value and mutex.
 *
 * This ensures that thread-MPI and OpenMP builds can't accidentally
 * initialize it more than once. */
//! @{
gmx_unused bool g_clfftInitialized = false;
gmx::Mutex      g_clfftMutex;
//! @}

} // namespace

ClfftInitializer::ClfftInitializer()
{
#if GMX_GPU == GMX_GPU_OPENCL
    gmx::lock_guard<gmx::Mutex> guard(g_clfftMutex);
    clfftSetupData              fftSetup;
    int                         initErrorCode = clfftInitSetupData(&fftSetup);
    if (initErrorCode != 0)
    {
        GMX_THROW(InternalError(formatString(
                "Failed to initialize the clFFT library, error code %d", initErrorCode)));
    }
    initErrorCode = clfftSetup(&fftSetup);
    if (initErrorCode != 0)
    {
        GMX_THROW(InternalError(formatString(
                "Failed to initialize the clFFT library, error code %d", initErrorCode)));
    }
    g_clfftInitialized = true;
#else
    GMX_UNUSED_VALUE(g_clfftInitialized);
#endif
}

ClfftInitializer::~ClfftInitializer()
{
#if GMX_GPU == GMX_GPU_OPENCL
    gmx::lock_guard<gmx::Mutex> guard(g_clfftMutex);
    if (g_clfftInitialized)
    {
        clfftTeardown();
        // TODO: log non-zero return values (errors)
    }
    g_clfftInitialized = false;
#endif
}

} // namespace gmx
