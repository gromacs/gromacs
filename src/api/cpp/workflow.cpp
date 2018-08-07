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
// Created by Eric Irrgang on 11/27/17.
//

#include "workflow.h"

#include <cassert>

#include "workflow-impl.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/status.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{


NodeSpecification::~NodeSpecification() = default;


std::unique_ptr<NodeSpecification> MDNodeSpecification::clone()
{
    assert(!tprfilename_.empty());
    std::unique_ptr<NodeSpecification> node {
        nullptr
    };
    node = gmx::compat::make_unique<MDNodeSpecification>(tprfilename_);
    return node;
}

MDNodeSpecification::MDNodeSpecification(std::string filename) :
    tprfilename_ {std::move(filename)}
{
    assert(!tprfilename_.empty());
}

NodeSpecification::paramsType MDNodeSpecification::params() const noexcept
{
    return tprfilename_;
}

NodeKey Workflow::addNode(std::unique_ptr<NodeSpecification> &&spec) noexcept
{
    (void)spec;
    return {};
}

std::unique_ptr<Workflow> Workflow::create(const std::string &filename)
{
    std::string name {
        "MD"
    };
    auto           spec = gmx::compat::make_unique<MDNodeSpecification>(filename);
    Workflow::Impl graph;
    graph.emplace(std::make_pair(name, std::move(spec)));
    auto           workflow = gmx::compat::make_unique<Workflow>(std::move(graph));
    return workflow;
}

std::unique_ptr<NodeSpecification> Workflow::getNode(const NodeKey &key) const noexcept
{
    const Impl &graph = graph_;
    assert((graph.count(key) == 0) || (graph.count(key) == 1));
    auto const  iter = graph.find(key);
    // Can return a null NodeSpecification if key is not found...
//    if (iter == graph.end())
//    {
//        auto const error = WorkflowKeyError(std::move(std::string(key)));
//        throw error;
//    }
    std::unique_ptr<NodeSpecification> node {
        nullptr
    };
    if (iter != graph.end())
    {
        node = iter->second->clone();
    }
    return node;
}

Workflow::Workflow(Workflow::Impl &&impl) :
    graph_ {std::forward<Workflow::Impl>(impl)}
{}

Workflow::Impl::const_iterator
Workflow::cbegin() const
{
    return graph_.cbegin();
}

Workflow::Impl::const_iterator
Workflow::cend() const
{
    return graph_.cend();
}

Workflow::Impl::const_iterator
Workflow::begin() const
{
    return graph_.cbegin();
}

Workflow::Impl::const_iterator
Workflow::end() const
{
    return graph_.cend();
}

} // end namespace gmxapi
