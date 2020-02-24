/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
#ifndef GMX_EWALD_TEST_HARDWARE_CONTEXT_H
#define GMX_EWALD_TEST_HARDWARE_CONTEXT_H

/*! \internal \file
 * \brief
 * Describes test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_ewald
 */

#include <map>
#include <string>
#include <vector>

#include "gromacs/ewald/pme_gpu_program.h"
#include "gromacs/utility/gmxassert.h"

class DeviceContext;
struct DeviceInformation;
class DeviceStream;

namespace gmx
{
namespace test
{
//! Hardware code path being tested
enum class CodePath
{
    CPU,
    GPU
};

//! Return a string useful for human-readable messages describing a \c codePath.
const char* codePathToString(CodePath codePath);

/*! \internal \brief
 * A structure to describe a hardware context  that persists over the lifetime
 * of the test binary - an abstraction over PmeGpuProgram with a human-readable string.
 */
struct TestHardwareContext
{
    //! Hardware path for the code being tested.
    CodePath codePath_;
    //! Readable description
    std::string description_;
    //! Device context
    DeviceContext* deviceContext_ = nullptr;
    //! Device stream
    DeviceStream* deviceStream_ = nullptr;
    //! Persistent compiled GPU kernels for PME.
    PmeGpuProgramStorage program_;

public:
    //! Retuns the code path for this context.
    CodePath codePath() const { return codePath_; }
    //! Returns a human-readable context description line
    std::string description() const { return description_; }
    //! Returns the device info pointer
    const DeviceInformation* deviceInfo() const;
    //! Get the device context
    const DeviceContext* deviceContext() const;
    //! Get the device stream
    const DeviceStream* deviceStream() const;
    //! Returns the persistent PME GPU kernels
    const PmeGpuProgram* pmeGpuProgram() const { return program_.get(); }
    //! Constructs the context for CPU builds
    TestHardwareContext(CodePath codePath, const char* description);
    //! Constructs the context for GPU builds
    TestHardwareContext(CodePath codePath, const char* description, const DeviceInformation& deviceInfo);
    //! Destructor
    ~TestHardwareContext();
};

} // namespace test
} // namespace gmx
#endif
