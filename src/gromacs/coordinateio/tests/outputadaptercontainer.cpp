/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Tests for outputadaptercontainer.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "gromacs/coordinateio/outputadaptercontainer.h"

#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/coordinateio/outputadapters/dummymodule.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace test
{

TEST(OutputAdapterContainer, MakeEmpty)
{
    OutputAdapterContainer container(efBaseOutputManager);
    EXPECT_TRUE(container.isEmpty());
}

TEST(OutputAdapterContainer, AddAdapter)
{
    OutputAdapterContainer container(efBaseOutputManager);
    container.addAdapter(
            compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                   efDummyModule));
    EXPECT_FALSE(container.isEmpty());

    for (const auto &adapter : container.getAdapters())
    {
        EXPECT_EQ(efDummyModule, adapter->getModuleIDFlag());
        EXPECT_EQ(efBaseOutputManager, adapter->getModuleRequirementFlag());
    }
}

TEST(OutputAdapterContainer, RejectBadAdapter)
{
    OutputAdapterContainer container(efBaseOutputManager);
    EXPECT_THROW(container.addAdapter(
                         compat::make_unique<DummyOutputModule>(efChangeVelocityModule,
                                                                efDummyModule)),
                 InconsistentInputError);
    EXPECT_TRUE(container.isEmpty());
}

TEST(OutputAdapterContainer, RejectDuplicateAdapter)
{
    OutputAdapterContainer container(efBaseOutputManager);
    EXPECT_NO_THROW(container.addAdapter(
                            compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                   efDummyModule)));
    EXPECT_FALSE(container.isEmpty());
    EXPECT_THROW(container.addAdapter(
                         compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                efDummyModule)),
                 InternalError);
}

TEST(OutputAdapterContainer, AcceptMultipleAdapters)
{
    OutputAdapterContainer container(efBaseOutputManager);
    EXPECT_NO_THROW(container.addAdapter(
                            compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                   efDummyModule)));
    EXPECT_FALSE(container.isEmpty());
    EXPECT_NO_THROW(container.addAdapter(
                            compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                   efChangeVelocityModule)));
    EXPECT_FALSE(container.isEmpty());
}

TEST(OutputAdapterContainer, TriggersOnBadOrder)
{
    OutputAdapterContainer container(efBaseOutputManager);
    EXPECT_NO_THROW(container.addAdapter(
                            compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                   efChangeCoordinateSelectionModule)));
    EXPECT_FALSE(container.isEmpty());
    EXPECT_THROW(container.addAdapter(
                         compat::make_unique<DummyOutputModule>(efBaseOutputManager,
                                                                efDummyModule)),
                 InternalError);
}

} // namespace test

} // namespace gmx
