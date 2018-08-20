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
//
// Created by Eric Irrgang on 11/28/17.
//

#include <memory>

#include "testingconfiguration.h"
#include "workflow.h"
#include "workflow-impl.h"
#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"
#include "gmxapi/md/mdmodule.h"
#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"

namespace
{

const auto filename = gmxapi::testing::sample_tprfilename;

// Create a work spec, then the implementation graph, then the container
TEST(ApiWorkflowImpl, Build)
{
    // Create work spec
    auto node = gmx::compat::make_unique<gmxapi::MDNodeSpecification>(filename);
    ASSERT_NE(node, nullptr);

    // Create key
    std::string key {
        "MD"
    };
    key.append(filename);

    // Create graph (workflow implementation object)
    gmxapi::Workflow::Impl impl;
    impl[key] = std::move(node);
    ASSERT_EQ(impl.count(key), 1);
    ASSERT_EQ(impl.size(), 1);

    // Create workflow container
    gmxapi::Workflow work {
        std::move(impl)
    };
}

TEST(ApiWorkflow, Creation)
{
    // Create from create() method(s)
    auto work = gmxapi::Workflow::create(filename);
    ASSERT_NE(work, nullptr);
}

TEST(ApiWorkflow, Accessors)
{
    auto work = gmxapi::Workflow::create(filename);
//    work->addNode()
//    work->getNode()
}

} // end anonymous namespace
