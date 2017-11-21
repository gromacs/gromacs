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
 * Tests for GPU memory status checker isHostMemoryPinned() being correct.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */
#include "gmxpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gputest.h"

namespace gmx
{

namespace test
{

namespace
{

//! Test fixture
using PinnedMemoryCheckerTest = GpuTest;

TEST_F(PinnedMemoryCheckerTest, DefaultContainerIsRecognized)
{
    std::vector<real> dummy(3, 1.5);
    EXPECT_FALSE(isHostMemoryPinned(dummy.data()));
}

TEST_F(PinnedMemoryCheckerTest, NonpinnedContainerIsRecognized)
{
    HostVector<real> dummy(3, 1.5);
    changePinningPolicy(&dummy, PinningPolicy::CannotBePinned);
    EXPECT_FALSE(isHostMemoryPinned(dummy.data()));
}

TEST_F(PinnedMemoryCheckerTest, PinnedContainerIsRecognized)
{
    HostVector<real> dummy(3, 1.5);
    changePinningPolicy(&dummy, PinningPolicy::CanBePinned);
    EXPECT_TRUE(isHostMemoryPinned(dummy.data()));
}

TEST_F(PinnedMemoryCheckerTest, DefaultCBufferIsRecognized)
{
    real *dummy;
    snew(dummy, 3);
    EXPECT_FALSE(isHostMemoryPinned(dummy));
    sfree(dummy);
}

TEST_F(PinnedMemoryCheckerTest, PinnedCBufferIsRecognized)
{
    real *dummy = nullptr;
    pmalloc((void **)&dummy, 3 * sizeof(real));
    EXPECT_TRUE(isHostMemoryPinned(dummy));
    pfree(dummy);
}

} // namespace
} // namespace
} // namespace
