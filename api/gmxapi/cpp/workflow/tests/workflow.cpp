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

#include "workflow.h"

#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "programs/mdrun/tests/moduletest.h"

#include "testingconfiguration.h"
#include "workflow_impl.h"

namespace gmxapi
{

namespace testing
{

namespace
{

//! Create a work spec, then the implementation graph, then the container
TEST_F(GmxApiTest, BuildApiWorkflowImpl)
{
    makeTprFile(100);
    // Create work spec
    auto node = std::make_unique<gmxapi::MDNodeSpecification>(runner_.tprFileName_);
    EXPECT_NE(node, nullptr);

    // Create key
    std::string key{ "MD" };
    key.append(runner_.tprFileName_);

    // Create graph (workflow implementation object)
    gmxapi::Workflow::Impl impl;
    impl[key] = std::move(node);
    EXPECT_EQ(impl.count(key), 1);
    EXPECT_EQ(impl.size(), 1);

    // Create workflow container
    EXPECT_NO_THROW(gmxapi::Workflow work{ std::move(impl) });
}

//! Create from create() method(s)
TEST_F(GmxApiTest, CreateApiWorkflow)
{
    makeTprFile(100);
    auto work = gmxapi::Workflow::create(runner_.tprFileName_);
    EXPECT_NE(work, nullptr);
}

} // end anonymous namespace

} // end namespace testing

} // end namespace gmxapi
