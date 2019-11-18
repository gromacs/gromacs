/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "threadaffinitytest.h"

#include "config.h"

#include <memory>

#include <gmock/gmock.h>

#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
namespace test
{

MockThreadAffinityAccess::MockThreadAffinityAccess() : supported_(true)
{
    using ::testing::_;
    using ::testing::Return;
#ifndef __clang_analyzer__
    ON_CALL(*this, setCurrentThreadAffinityToCore(_)).WillByDefault(Return(true));
#endif
}

MockThreadAffinityAccess::~MockThreadAffinityAccess() {}


ThreadAffinityTestHelper::ThreadAffinityTestHelper()
{
    snew(cr_, 1);
    cr_->nnodes = gmx_node_num();
    cr_->nodeid = gmx_node_rank();
    cr_->duty   = DUTY_PP;
#if GMX_MPI
    cr_->mpi_comm_mysim = MPI_COMM_WORLD;
#endif
    hwOpt_.threadAffinity      = ThreadAffinity::Auto;
    hwOpt_.totNumThreadsIsAuto = false;
    physicalNodeId_            = 0;
}

ThreadAffinityTestHelper::~ThreadAffinityTestHelper()
{
    sfree(cr_);
}

void ThreadAffinityTestHelper::setLogicalProcessorCount(int logicalProcessorCount)
{
    hwTop_ = std::make_unique<HardwareTopology>(logicalProcessorCount);
}

} // namespace test
} // namespace gmx
