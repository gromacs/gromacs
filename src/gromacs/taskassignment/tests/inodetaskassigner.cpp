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
/*! \internal \file
 * \brief
 * Tests for node task-assigner factory.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/inodetaskassigner.h"

#include <string>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/taskassignment/nonbondedoncpu.h"
#include "gromacs/taskassignment/nonbondedongpuauto.h"
#include "gromacs/taskassignment/nonbondedongpufromuser.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace test
{
namespace
{

//! Has the task assigner got the correct type?
template <typename AssignerType>
bool taskAssignerIs(INodeTaskAssigner *taskAssigner)
{
    auto rawPointer = dynamic_cast<AssignerType *>(taskAssigner);
    return rawPointer != nullptr;
}

//! Helper function to expect that the right kind of task assigner is made.
template <typename AssignerType, typename ... Args>
void makeNodeTaskAssignerOf(Args && ... args)
{
    std::unique_ptr<INodeTaskAssigner> taskAssigner;
    EXPECT_NO_THROW(taskAssigner = makeNodeTaskAssigner(std::forward<Args>(args) ...));
    EXPECT_TRUE(taskAssignerIs<AssignerType>(taskAssigner.get()));
}

//! Helper function to expect that the right kind of exception is thrown.
template <typename ExceptionType, typename ... Args>
void makeNodeTaskAssignerThrows(Args && ... args)
{
    std::unique_ptr<INodeTaskAssigner> taskAssigner;
    EXPECT_THROW(makeNodeTaskAssigner(std::forward<Args>(args) ...), ExceptionType);
}

// Tests the four different behaviours resulting from the 12 different
// code paths through makeNodeTaskAssigner (three -nb options,
// empty/filled -gpu_id, GMX_EMULATE_GPU on/off).
TEST(MakeNodeTaskAssigner, WorksAsExpected)
{
    // The factory merely passes the hardware info through to (some of) the
    // instantiated objects, so it is OK to leave it
    // default-initialized, and re-use it.
    gmx_hw_info_t dummyHardwareInfo;

    // Strings for user assignments of GPU IDs - only the state of "emptiness" is relevant for the factory.
    std::string anything("anything"), empty("");

    makeNodeTaskAssignerOf<NonbondedOnCpu>(            "cpu",  empty,    EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnCpu>(            "cpu",  empty,    EmulateGpuNonbonded::Yes, dummyHardwareInfo);
    makeNodeTaskAssignerThrows<InconsistentInputError>("cpu",  anything, EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerThrows<InconsistentInputError>("cpu",  anything, EmulateGpuNonbonded::Yes, dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnGpuAuto>(        "auto", empty,    EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnCpu>(            "auto", empty,    EmulateGpuNonbonded::Yes, dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnGpuFromUser>(    "auto", anything, EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerThrows<InconsistentInputError>("auto", anything, EmulateGpuNonbonded::Yes, dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnGpuAuto>(        "gpu",  empty,    EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerThrows<InconsistentInputError>("gpu",  empty,    EmulateGpuNonbonded::Yes, dummyHardwareInfo);
    makeNodeTaskAssignerOf<NonbondedOnGpuFromUser>(    "gpu",  anything, EmulateGpuNonbonded::No,  dummyHardwareInfo);
    makeNodeTaskAssignerThrows<InconsistentInputError>("gpu",  anything, EmulateGpuNonbonded::Yes, dummyHardwareInfo);
}

} // namespace
} // namespace
} // namespace
