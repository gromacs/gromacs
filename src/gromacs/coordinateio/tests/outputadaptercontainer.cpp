/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*!\internal
 * \file
 * \brief
 * Tests for outputadaptercontainer.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "gromacs/coordinateio/outputadaptercontainer.h"

#include <algorithm>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/utility/exceptions.h"

#include "testmodule.h"

namespace gmx
{

namespace test
{

TEST(OutputAdapterContainer, MakeEmpty)
{
    OutputAdapterContainer container(CoordinateFileFlags::Base);
    EXPECT_TRUE(container.isEmpty());
}

TEST(OutputAdapterContainer, AddAdapter)
{
    OutputAdapterContainer container(CoordinateFileFlags::Base);
    container.addAdapter(std::make_unique<DummyOutputModule>(CoordinateFileFlags::Base),
                         CoordinateFileFlags::RequireNewFrameStartTime);
    EXPECT_FALSE(container.isEmpty());
}

TEST(OutputAdapterContainer, RejectBadAdapter)
{
    OutputAdapterContainer container(CoordinateFileFlags::Base);
    EXPECT_THROW(container.addAdapter(
                         std::make_unique<DummyOutputModule>(CoordinateFileFlags::RequireVelocityOutput),
                         CoordinateFileFlags::RequireVelocityOutput),
                 InconsistentInputError);
    EXPECT_TRUE(container.isEmpty());
}

TEST(OutputAdapterContainer, RejectDuplicateAdapter)
{
    OutputAdapterContainer container(CoordinateFileFlags::Base);
    EXPECT_NO_THROW(container.addAdapter(std::make_unique<DummyOutputModule>(CoordinateFileFlags::Base),
                                         CoordinateFileFlags::RequireNewFrameStartTime));
    EXPECT_FALSE(container.isEmpty());
    EXPECT_THROW(container.addAdapter(std::make_unique<DummyOutputModule>(CoordinateFileFlags::Base),
                                      CoordinateFileFlags::RequireNewFrameStartTime),
                 InternalError);
}

TEST(OutputAdapterContainer, AcceptMultipleAdapters)
{
    OutputAdapterContainer container(CoordinateFileFlags::Base);
    EXPECT_NO_THROW(container.addAdapter(std::make_unique<DummyOutputModule>(CoordinateFileFlags::Base),
                                         CoordinateFileFlags::RequireForceOutput));
    EXPECT_FALSE(container.isEmpty());
    EXPECT_NO_THROW(container.addAdapter(std::make_unique<DummyOutputModule>(CoordinateFileFlags::Base),
                                         CoordinateFileFlags::RequireVelocityOutput));
    EXPECT_FALSE(container.isEmpty());
}

} // namespace test

} // namespace gmx
