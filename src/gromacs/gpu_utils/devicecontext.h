/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 *  \brief Declare type for handling GPU device contexts with RAII style.
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_DEVICECONTEXT_H
#define GMX_GPU_UTILS_DEVICECONTEXT_H

#include <memory>

#include "gromacs/utility/mutex.h"

struct gmx_device_info_t;

namespace gmx
{
class MDLogger;

/*! \libinternal
 * \brief Implements RAII for handling device context lifetime. */
class DeviceContext
{
    public:
        //! Default constructor.
        DeviceContext();
        //! Creates the context.
        DeviceContext(const MDLogger &mdlog, gmx_device_info_t *deviceInfo);
        //! Frees the context, if any.
        ~DeviceContext();
        //! Getter for underlying device info.
        gmx_device_info_t *getDeviceInfo();
        //! Const getter for underlying device info.
        const gmx_device_info_t *getDeviceInfo() const;
    private:
        /*! \brief The underlying (non-owned) device info struct.
         *
         * \todo deviceInfo_ should be logically const, but currently
         * init_gpu modifies it to set up NVML support. This could
         * happen during the detection phase, and deviceInfo_ could
         * the become const. */
        gmx_device_info_t *deviceInfo_;
        //! Mutex to ensure thread-MPI shared contexts is free of races.
        static Mutex       mutex_;
};

/*! \libinternal \brief Convenience typedef. */
using DeviceContextPtr = std::shared_ptr<DeviceContext>;

} // namespace

#endif
