/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_NVML_CUH
#define GMX_GPU_UTILS_NVML_CUH

#include "config.h"

#if HAVE_NVML
#include <nvml.h>
#endif

#include <string>

namespace gmx
{

/*! \brief Exception class to throw when NVML errors occur.
 *
 * This might be used when the NVML API returns an error code, or when
 * something else is inconsistent (e.g. a CUDA device that should have
 * an NVML device ID doesn't have one).
 *
 * \todo This would be improved by replacing it with
 * gmx::InternalError (or subclassed from it) but we can't do that
 * until we can have C++11 host-side CUDA code. */
class NvmlException
{
    public:
        explicit NvmlException(const std::string &message) : message_(message) {
        }
    private:
        std::string message_;
};

/*! Manages interactions with NVML to set application clocks on GPUs to their maximum values.
 *
 * The setup() and changeClocks() methods should always be called
 * exactly once during initialization, and resetClocks() exactly once
 * during clean-up. No other method should be called after
 * resetClocks(). */
class NvmlManager
{
    public:
        /*! Determine whether changing GPU clocks is possible.
         *
         * \param[in]  prop     CUDA device properties
         * \param[out] message  If non-empty upon return, a message to report to the log file
         *
         * \throws    NvmlException   If an NVML API call failed, or a suitable device doesn't have an NVML device ID.
         *            std::bad_alloc  If memory allocation fails.
         *
         * \todo This method might be better as a factory method, but
         * first some more C++ conversion of GPU-related code needs to
         * take place.
         */
        void setup(const cudaDeviceProp &prop, std::string *message);
        //! Getter
        bool getClocksCanBeChanged() const { return clocksCanBeChanged_; };
        /*! Change GPU clocks (if possible).
         *
         * \param[in]  prop     CUDA device properties
         * \param[out] message  If non-empty upon return, a message to report to the log file
         *
         * \throws    NvmlException   If an NVML API call failed
         *            std::bad_alloc  If memory allocation fails.
         */
        void changeClocks(std::string *message);
        /*! Reset GPU clocks to the original values (if possible and appropriate).
         *
         * We reset them only if we changed them, and they haven't
         * been changed in the meantime. An NvmlManager is no longer
         * valid for use after calling this method.
         *
         * \throws    NvmlException  If an NVML API call failed
         */
        void resetClocks();
    private:
        //! Whether application clocks can be changed
        bool                clocksCanBeChanged_;
        //! Whether application clocks were actually changed
        bool                clocksWereChanged_;
        //! The original SM clock before we (perhaps) changed it
        unsigned int        originalSmClock_;
        //! The original memory clock before we (perhaps) changed it
        unsigned int        originalMemoryClock_;
        //! The SM clock we set (if any)
        unsigned int        ourSmClock_;
        //! The memory clock we set (if any)
        unsigned int        ourMemoryClock_;
        /*! \brief Descriptive name for this device, e.g. obtained from cudaDeviceProp name.
         *
         * \todo Until we can make a copy in a std::string, the caller
         * must ensure that the lifetime of prop used in setup()
         * exceeds that of this object. */
        const char         *name_;
#if HAVE_NVML
        //! The ID of this NVML device (not necessarily the same as the CUDA device ID)
        nvmlDevice_t        deviceId_;
#endif  /* HAVE_NVML */
};

} // namespace

#endif
