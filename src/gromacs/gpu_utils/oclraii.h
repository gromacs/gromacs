/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \libinternal \file
 * \brief Declare RAII helpers for OpenCL types, along with
 * supporting type traits.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_OCLRAII_H
#define GMX_GPU_UTILS_OCLRAII_H

#include "gromacs/gpu_utils/gmxopencl.h"

namespace gmx
{

/*! \libinternal \brief Stub for OpenCL type traits */
template<typename cl_type>
struct OpenClTraits;

/*! \libinternal \brief Implements common trait infrastructure for OpenCL types. */
template<typename cl_type>
struct OpenClTraitsBase
{
    //! Type of the function that will release a handle of this type.
    using ReleaserType = cl_int (*)(cl_type);
};

/*! \libinternal \brief Implements traits for cl_context. */
template<>
struct OpenClTraits<cl_context> : public OpenClTraitsBase<cl_context>
{
    //! Function that will release a handle of this type.
    static constexpr ReleaserType releaser = clReleaseContext;
};

/*! \libinternal \brief Implements traits for cl_command_queue. */
template<>
struct OpenClTraits<cl_command_queue> : public OpenClTraitsBase<cl_command_queue>
{
    //! Function that will release a handle of this type.
    static constexpr ReleaserType releaser = clReleaseCommandQueue;
};

/*! \libinternal \brief Implements traits for cl_program. */
template<>
struct OpenClTraits<cl_program> : public OpenClTraitsBase<cl_program>
{
    //! Function that will release a handle of this type.
    static constexpr ReleaserType releaser = clReleaseProgram;
};

/*! \libinternal \brief Implements traits for cl_kernel. */
template<>
struct OpenClTraits<cl_kernel> : public OpenClTraitsBase<cl_kernel>
{
    //! Function that will release a handle of this type.
    static constexpr ReleaserType releaser = clReleaseKernel;
};

/*! \libinternal \brief Wrapper of OpenCL type \c cl_type to implement RAII.
 *
 * Works by calling the releaser function associated with cl_type
 * by OpenClTraits.
 *
 * Simple copying and assignment are not supported, because there's no
 * need for that, and would require OpenCL API calls for deep copies
 * if they were needed. Move and move assignment are fine, however. */
template<typename cl_type>
class ClHandle
{
public:
    //! Constructor that takes an already created handle.
    explicit ClHandle(cl_type handle) : handle_(handle) {}
    //! Destructor that calls the releaser associated with cl_type.
    ~ClHandle() { OpenClTraits<cl_type>::releaser(handle_); }
    //! Deleted default constructor.
    ClHandle() = delete;
    //! Deleted assignment operator.
    ClHandle& operator=(const ClHandle&) = delete;
    //! Deleted copy constructor.
    ClHandle(const ClHandle&) = delete;
    //! Default move assignment operator.
    ClHandle& operator=(ClHandle&&) noexcept = default;
    //! Default copy constructor.
    ClHandle(ClHandle&&) noexcept = default;
    /*! \brief Convenience conversion operator so the wrapper type
     * can simply convert to the wrapped type. */
    operator cl_type() const { return handle_; }

private:
    //! The wrapped object.
    cl_type handle_;
};

//! Convenience declarations.
/*! @{ */
using ClContext      = ClHandle<cl_context>;
using ClCommandQueue = ClHandle<cl_command_queue>;
using ClProgram      = ClHandle<cl_program>;
using ClKernel       = ClHandle<cl_kernel>;
/*! @} */

} // namespace gmx

#endif
