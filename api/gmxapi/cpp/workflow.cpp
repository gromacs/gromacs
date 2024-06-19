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

#include "gromacs/utility/gmxassert.h"

#include "gmxapi/exceptions.h"

#include "workflow_impl.h"

namespace gmxapi
{

NodeSpecification::~NodeSpecification() = default;


std::unique_ptr<NodeSpecification> MDNodeSpecification::clone()
{
    GMX_ASSERT(!tprfilename_.empty(), "Need a non-empty filename string.");
    std::unique_ptr<NodeSpecification> node = nullptr;
    node                                    = std::make_unique<MDNodeSpecification>(tprfilename_);
    return node;
}

MDNodeSpecification::MDNodeSpecification(const std::string& filename) : tprfilename_{ filename }
{
    GMX_ASSERT(!tprfilename_.empty(), "Need a non-empty filename string.");
}

NodeSpecification::paramsType MDNodeSpecification::params() const noexcept
{
    return tprfilename_;
}

NodeKey Workflow::addNode(std::unique_ptr<NodeSpecification> spec)
{
    // TODO capture provided NodeSpecification.
    // Relates to gmxapi milestone 7, described at https://gitlab.com/gromacs/gromacs/-/issues/2585
    throw gmxapi::MissingImplementationError("Member function not yet implemented or used.");
    (void)spec;
    return {};
}

std::unique_ptr<Workflow> Workflow::create(const std::string& filename)
{
    const std::string name = "MD";
    auto              spec = std::make_unique<MDNodeSpecification>(filename);
    Workflow::Impl    graph;
    graph.emplace(std::make_pair(name, std::move(spec)));
    auto workflow = std::make_unique<Workflow>(std::move(graph));
    return workflow;
}

std::unique_ptr<NodeSpecification> Workflow::getNode(const NodeKey& key) const noexcept
{
    const Impl& graph = graph_;
    GMX_ASSERT((graph.count(key) == 0) || (graph.count(key) == 1),
               "Key should occur zero or one times.");
    auto const                         iter = graph.find(key);
    std::unique_ptr<NodeSpecification> node = nullptr;
    if (iter == graph.end())
    {
        // key not found. Return nullptr.
        // Alternatively, we could throw a WorkflowKeyError.
    }
    else
    {
        node = (*iter).second->clone();
    }
    return node;
}

Workflow::Workflow(Workflow::Impl&& impl) : graph_{ std::forward<Workflow::Impl>(impl) } {}

Workflow::Impl::const_iterator Workflow::cbegin() const
{
    return graph_.cbegin();
}

Workflow::Impl::const_iterator Workflow::cend() const
{
    return graph_.cend();
}

Workflow::Impl::const_iterator Workflow::begin() const
{
    return graph_.cbegin();
}

Workflow::Impl::const_iterator Workflow::end() const
{
    return graph_.cend();
}

} // end namespace gmxapi
