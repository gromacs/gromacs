/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#ifndef GMXAPI_CONTEXT_H
#define GMXAPI_CONTEXT_H
/*! \file
 * \brief Declares classes representing execution context.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <memory>
#include <string>
#include <vector>

namespace gmxapi
{

class Workflow;
class Session;

/*!
 * \brief Container for arguments passed to the simulation runner.
 *
 * This is part of an implementation that essentially wraps the command line
 * implementation at a high level.
 * \todo Modernize expression and handling of MD user options.
 */
using MDArgs = std::vector<std::string>;

/*!
 * \brief Context implementation abstract base class.
 *
 * Context Implementations depend on the execution environment, hardware resources, and
 * possibly other factors, and so are not constructed directly but by helper functions.
 * Their details are not exposed at the high level API.
 * \ingroup gmxapi
 */
class ContextImpl;

/*!
 * \brief Execution context.
 *
 * The execution context represents computing resources and zero, one, or more
 * workflows to execute. All API objects exist in some context, which determines
 * how the API objects interact underneath the hood.
 *
 * A proxy can be configured with information needed to initialize a runtime
 * environment capable of executing a work load, independently of defining the
 * work.
 * The actual execution
 * environment is not necessarily instantiated / configured until the work is
 * performed.
 * Thus, construction of a Context object does not necessarily imply
 * initialization of compute resources, but any active compute resources are
 * appropriately deinitialized when the object is destroyed. However, to allow
 * opportunities for exception handling, resources should be deinitialized when
 * and as work is completed by the Runner.
 *
 * Ultimately, the class in this header file should just define an interface.
 * Implementations for different execution environments will be provided by the
 * library or extension code and documented separately,
 * and it should be straight-forward to
 * write extension code for other execution environments.
 *
 * \todo Further encapsulate resource state transitions with Session objects.
 *
 * \ingroup gmxapi
 */
class Context
{
public:
    /*!
     * \brief Get a handle to a new default context object.
     */
    Context();
    ~Context();

    /*!
     * \brief Nearly trivial copy
     *
     * \{
     */
    Context(const Context&) = default;
    Context& operator=(const Context&) = default;
    //! \}

    /*!
     * \brief Allow move
     *
     * \{
     */
    Context(Context&&) = default;
    Context& operator=(Context&&) = default;
    //! \}

    /*!
     * \brief Construct by wrapping an implementation object.
     *
     * \param impl Ownership of Context definition and resources.
     */
    explicit Context(std::shared_ptr<ContextImpl> impl);

    /*!
     * \brief Set the simulation runtime arguments for this instance.
     *
     * \param mdArgs User-provided runtime parameters (such as for `gmx mdrun`)
     *
     * This is awkwardly named and due for some evolution, since most of the mdrun CLI options
     * pertain to the execution environment rather than the simulation parameters. For the first
     * implementation, we just map user arguments to the equivalent command-line substrings.
     */
    void setMDArgs(const MDArgs& mdArgs);

    /*!
     * \brief Launch a workflow in the current context, if possible.
     *
     * \param work Configured workflow to instantiate.
     * \return Ownership of a new session or nullptr if not possible.
     */
    std::shared_ptr<Session> launch(const Workflow& work);

private:
    /*!
     * \brief Private implementation
     *
     * Early implementation has not yet developed distinct handle classes
     * for different levels of access to the Context. In this implementation,
     * new Context handles to the same resources can be created (in library code)
     * by constructing from a smart pointer to an implementation object.
     *
     * \todo Consider client requirements and ownership contracts.
     * Whether/how API client might be required to create and hold a Context
     * and/or Session object for the length of a process.
     */
    std::shared_ptr<ContextImpl> impl_;
};

} // end namespace gmxapi

#endif // header guard
