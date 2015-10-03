/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief Tests for gmx::AlignedAllocator
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/alignedallocator.h"

#include <list>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

namespace gmx
{

TEST(AlignedAllocatorTest, AllocatorAlign)
{
    std::size_t mask;
    real *      p;

    // Alignment 1 element
    AlignedAllocator<real, sizeof(real)>      a1;
    mask = static_cast<std::size_t>(sizeof(real)-1);
    p    = a1.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);

    // Alignment 2 elements
    AlignedAllocator<real, 2*sizeof(real)>   a2;
    mask = static_cast<std::size_t>(2*sizeof(real)-1);
    p    = a2.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);

    // Alignment 4 elements
    AlignedAllocator<real, 4*sizeof(real)>   a4;
    mask = static_cast<std::size_t>(4*sizeof(real)-1);
    p    = a4.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);

    // Alignment 8 elements
    AlignedAllocator<real, 8*sizeof(real)>   a8;
    mask = static_cast<std::size_t>(8*sizeof(real)-1);
    p    = a8.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);

    // Alignment 16 elements
    AlignedAllocator<real, 16*sizeof(real)>   a16;
    mask = static_cast<std::size_t>(16*sizeof(real)-1);
    p    = a16.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);

    // Alignment 32 elements
    AlignedAllocator<real, 32*sizeof(real)>   a32;
    mask = static_cast<std::size_t>(32*sizeof(real)-1);
    p    = a32.allocate(1000);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(p) & mask);
    a1.deallocate(p, 1000);
}


TEST(AlignedAllocator, Vector)
{
    // Here we just use a single value for alignments that is so large that it is unlikely
    // to happen by chance with the system malloc() - 32 elements.
    std::size_t mask = static_cast<std::size_t>(32*sizeof(real)-1);

    std::vector<real, AlignedAllocator<real, 32*sizeof(real)> > v(10);
    EXPECT_EQ(0, reinterpret_cast<std::size_t>(v.data()) & mask);

    for (std::size_t i = 10000; i <= 100000; i += 10000)
    {
        v.resize(i);
        EXPECT_EQ(0, reinterpret_cast<std::size_t>(v.data()) & mask);
    }
}


}
