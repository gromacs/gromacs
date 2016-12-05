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
#ifndef GMX_EWALD_TEST_HARDWARE_CONTEXTS_H
#define GMX_EWALD_TEST_HARDWARE_CONTEXTS_H

/*! \internal \file
 * \brief
 * Describes test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include <list>
#include <map>

#include <gtest/gtest.h>

#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
namespace test
{
//! Hardware code path being tested
enum class CodePath
{
    CPU,
    CUDA
};

/*! \internal \brief
 * A structure to describe a hardware context - an abstraction over
 * gmx_device_info_t with a human-readable string.
 * TODO: currently this does not know which CodePath it belongs too.
 * It probably should! That would save us one loop in all the PME tests.
 */
struct TestHardwareContext
{
    //! Readable description
    std::string        description_;
    //! Device information pointer
    gmx_device_info_t *deviceInfo_;

    public:
        //! Returns a human-readable context description line
        std::string getDescription() const{return description_; }
//! Returns the device info pointer
        gmx_device_info_t *getDeviceInfo() const{return deviceInfo_; }
        //! Constructs the context
        TestHardwareContext(const char *description, gmx_device_info_t *deviceInfo) : description_(description), deviceInfo_(deviceInfo){}
};

//! A list of hardware contexts
typedef std::list<TestHardwareContext> TestHardwareContexts;

/*! \internal \brief
 * This class performs one-time test initialization (enumerating the hardware)
 */
class PmeTestEnvironment : public ::testing::Environment
{
    private:
        //! General hardware info
        unique_cptr<gmx_hw_info_t, gmx_hardware_info_free> hardwareInfo_;
        //! Storage of hardware contexts
        std::map<CodePath, TestHardwareContexts>           hardwareContextsByMode_;
        //! Simple GPU initialization, allowing for PME to work on GPU
        void hardwareInit();

    public:
        //! Default
        ~PmeTestEnvironment() = default;
        //! This is called by GTest framework once to query the hardware
        void SetUp();
        //! Get available hardware contexts for given code path
        const TestHardwareContexts &getHardwareContexts(CodePath mode){return hardwareContextsByMode_.at(mode); }
};

//! The test environment
extern PmeTestEnvironment * const pmeEnv;
}
}
#endif
