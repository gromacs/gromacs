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

#ifndef GMXAPI_WORKFLOW_H
#define GMXAPI_WORKFLOW_H

/*! \libinternal \file
 * \brief Declare public interface for Workflow and related infrastructure.
 *
 * \ingroup gmxapi
 */

#include <map>

#include <forward_list>

#include "gmxapi/gmxapi.h"

namespace gmxapi
{

/*!
 * \brief Uniquely identify a workflow node in the graph.
 *
 * The key probably needs a human-readable aspect, some machine-decipherable encoding of roles taken by the node,
 * and a hash to uniquely identify the output of the node (i.e. deterministic input parameters). It is probably not
 * necessary for nodes to refer to the consumers of their output by key, but they should abstractly refer to their
 * inputs by a key that is not dependent on a currently-running workflow.
 *
 * Requirements and roles:
 *
 * * serve as a key for use by other nodes to name their inputs
 * * encode workflow scheduling hints (TBD)
 * * provide robust assurance of reproducible results and restartability
 * * allow nodes to specify only their immediately dependent nodes (inwards directed edges)
 *
 * Workflow specifications need to be serializeable and portable across job restarts and porting to other computing
 * resources. The data graph manager and/or work scheduler need to be able to look at the inputs specified for a node
 * and be able to determine that the required node or its output is available. If a node is used as the input for
 * multiple other nodes, it should be clear how to avoid wasting resources when meeting the data requirement. If
 * similar looking nodes have different inputs or parameters, they must not be mistaken to be equivalent.
 *
 * Context-dependent aspects of the workflow specification cannot be included in a hash, then, but context-independent
 * aspects that affect the output of a node must be reflected.
 *
 * For example, an input filename should be included as identifying information, but the absolute path should not,
 * though path hints or conventions should be clear in the context. The filename is sufficient as a parameter with which
 * to construct the workflow node in an execution context, but is insufficient to uniquely identify the file since
 * several names get reused a lot. Some sort of checksum of the file should also be included so that the inputs of the
 * workflow at execution time can be checked against the inputs when the workflow was specified.
 *
 * Uniqueness of inputs could be more elaborate. For instance, a node may require the trajectory of a specific simulation
 * as input, but flexibly handle starting from an arbitrary step in that trajectory to allow check-pointed workflows.
 *
 * The workflow object can have a list of keys that can be instantiated with no input dependencies, the scheduler could
 * scan for keys that represent source nodes, or workflow containers could be turned into graphs through an additional
 * preprocessing or clustering phase, but it will be easiest if we assert a protocol such as a node is not instantiated
 * or activated until its inputs are ready.
 */
//class NodeKey final
//{
//    public:
//        std::string name();
//    private:
//        std::string name_;
//};
// Type alias until more elaborate implementation is needed.
using NodeKey = std::string;

// Forward declarations for definitions below.
class NodeKeyIterator;
class NodeSpecification;
class WorkflowKeyError;

/*!
 * \brief Recipe for a computational workflow.
 *
 * Provides a lightweight and portable container defining the nodes and edges in a workflow with enough information for
 * the workflow to be instantiated and run.
 *
 * \ingroup gmxapi
 */
class Workflow final
{
    public:
        using Impl = typename std::map < NodeKey, std::unique_ptr < NodeSpecification>>;

        /*! \brief Use create() to get Workflow objects.
         *
         * An empty workflow is not meaningful except to a builder, which does not
         * yet exist. Even a builder, though, will probably create the implementation
         * object directly and the Workflow object from that.
         */
        Workflow() = delete;

        /*!
         * \brief Construct from implementation
         *
         * \param impl
         */
        explicit Workflow(Impl &&impl);

        /*!
         * \brief Add a node to the workflow graph.
         *
         * The work specification must already have its inputs assigned to existing
         * nodes. This operation should only be permitted if it does not render a
         * valid workflow invalid.
         */
        NodeKey addNode(std::unique_ptr<NodeSpecification> &&spec) noexcept;

        /*!
         * \brief Get the node specification for a provided key.
         *
         * \param key Unique identifier for a node in the graph.
         * \return copy of the node specification.
         */
        std::unique_ptr<NodeSpecification> getNode(const gmxapi::NodeKey &key) const noexcept;

        /*!
         * \brief Get an iterator to the described node keys.
         *
         * The order in which the nodes are returned is unspecified. Only forward iterator is provided.
         *
         * \return iterator of keys for nodes that can be retrieved by the client if needed.
         *
         */
//        NodeKey::iterator keys() noexcept;
//        NodeKey::const_iterator keys() const noexcept;

        /*!
         * \brief Get an iterator to the node key--value pairs.
         *
         *
         * The order in which the nodes are returned is unspecified. Only forward iterator is provided.
         */
        Impl::const_iterator cbegin() const;
        Impl::const_iterator cend() const;
        // Allow range based for loop to work before C++17
        Impl::const_iterator begin() const;
        Impl::const_iterator end() const;

        /*!
         * \brief Create a new workflow.
         *
         * \param filename TPR filename accessible both to the client and library.
         * \return Ownership of a new Workflow instance.
         */
        static std::unique_ptr<Workflow> create(const std::string &filename);
    private:
        /*!
         * \brief Storage structure.
         */
        Impl graph_;

};

/*!
 * \brief Portable specification to define work and inform instantiation by the library.
 *
 * The GROMACS library creates the objects it needs to run as late as possible while
 * optimizing parallel resources at run time. The specifications provide a way for
 * client code to interact with the definition of the work to be performed while carrying
 * enough information for GROMACS to launch.
 *
 * Client input is translated into serializeable parameters sufficient to instantiate
 * the node at runtime.
 *
 * On the library side, the spec should have a pointer to a factory function for
 * the library object(s) it represents that is valid in the current Context. Thus,
 * when a workflow specification (and thus Node Specifications) are cloned to new
 * Contexts, the Contexts must resolve an appropriate function pointer or raise an
 * appropriate exception indicating the specified work is not possible on the targeted
 * execution context.
 *
 * Different node types will have different sorts of parameters and such. Should we
 * try to identify some sort of class hierarchy to aid control flow in looking for
 * parameters?
 */
class NodeSpecification
{
    public:
        /// Base class destructor.
        virtual ~NodeSpecification();

//        using paramsType = std::map<std::string, std::string>;
        using paramsType = std::string;

        /*!
         * \brief Get a copy of a node.
         *
         * \return ownership of a new node specification
         *
         * Future versions may use this function to translate a node spec from one
         * context to another.
         */
        virtual std::unique_ptr<NodeSpecification> clone() = 0;

        paramsType params_ {};

        virtual paramsType params() const noexcept = 0;

//        Node build() = 0;

};

}      //end namespace gmxapi

#endif //GMXAPI_WORKFLOW_H
