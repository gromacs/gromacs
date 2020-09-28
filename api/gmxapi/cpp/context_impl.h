/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMXAPI_CONTEXT_IMPL_H
#define GMXAPI_CONTEXT_IMPL_H
/*! \file
 * \brief Declare gmxapi::ContextImpl
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <memory>
#include <string>

#include "gromacs/mdrun/legacymdrunoptions.h"
#include "gromacs/mdtypes/mdrunoptions.h"

#include "gmxapi/context.h"
#include "gmxapi/session.h"

namespace gmxapi
{

/*!
 * \brief Provide RAII management of communications resource state.
 *
 * To acquire an MpiContextManager is to have assurance that any external MPI
 * environment is ready to use. When the MpiContextManager is released or
 * goes out of scope, the destructor finalizes the resources.
 *
 * Note that thread-MPI chooses the number of ranks and constructs its
 * MPI communicator internally, so does not and is unlikely to ever
 * participate here.
 *
 * \todo There is no resource for logging or reporting errors during initialization
 *
 * \ingroup gmxapi
 */
class MpiContextManager
{
public:
    MpiContextManager();

    ~MpiContextManager();

    /*!
     * \brief Exclusive ownership of a scoped context means copying is impossible.
     *
     * \{
     */
    MpiContextManager(const MpiContextManager&) = delete;
    MpiContextManager& operator=(const MpiContextManager&) = delete;
    //! \}

    /*!
     * \brief Move semantics are non-trivial.
     *
     * \{
     */
    MpiContextManager(MpiContextManager&&) noexcept = delete;
    MpiContextManager& operator=(MpiContextManager&&) noexcept = delete;
    //! \}
};

/*!
 * \brief Context implementation.
 *
 * Execution contexts have a uniform interface specified by the API. Implementations for
 * particular execution environments can specialize / derive from this base.
 *
 * \todo Separate interface and implementation.
 *
 * \warning Definition and semantics depend on configure-time details. For example,
 *          MPI-enabled libraries always hold a valid MPI communicator via MpiContextManager,
 *          whereas tMPI and non-MPI builds hold a meaningless MpiContextManager.
 *
 * \todo Provide functions or traits for introspection.
 *
 * \ingroup gmxapi
 */
class ContextImpl final : public std::enable_shared_from_this<ContextImpl>
{
public:
    ~ContextImpl();

    /*!
     * \brief Factory function
     *
     * Since this class provides `shared_from_this`, we need to make sure
     * that it never exists without a shared_ptr owning it.
     *
     * If we can confirm `shared_from_this` is no longer necessary, implementation may change.
     * \todo: Use registration/deregistration of launched Sessions to log warnings on shutdown
     *        instead of letting Session keep ContextImpl alive.
     *
     * \return ownership of a new object
     */
    static std::shared_ptr<ContextImpl> create();

    /*!
     * \brief Copy disallowed because Session state would become ambiguous.
     *
     * The API implementation needs to unambiguously determine
     * which Sessions and Contexts are associated with each other.
     * \{
     */
    ContextImpl(const ContextImpl&) = delete;
    ContextImpl& operator=(const ContextImpl&) = delete;
    //! \}

    /*!
     * \brief Objects are not trivial to move.
     *
     * \todo Implement move semantics. Requires a moveable const MpiContextManager and
     *       LegacyMdrunOptions members.
     *
     * \{
     */
    ContextImpl(ContextImpl&&) = delete;
    ContextImpl& operator=(ContextImpl&&) = delete;
    //! \}

    /*!
     * \brief Translate the workflow to the execution context and launch.
     *
     * \param work workflow graph
     * \return ownership of a new session
     *
     * \todo This probably makes more sense as a free function, but we need to determine access policies.
     *
     * Session is returned with shared_ptr ownership so that Context
     * can hold a weak_ptr and because Session Resources handles
     * are still evolving.
     * \todo Hide lifetime management and ownership from handle object.
     * We can achieve the necessary aspects of this shared_ptr at a lower level of implementation.
     * Note also that returned value policies can be implemented in higher level wrappers to ensure
     * correct object lifetime scope. (See pybind, for instance.)
     */
    std::shared_ptr<Session> launch(const Workflow& work);

    /*!
     * \brief Retain the ability to find a launched session while it exists.
     *
     * The client owns the Session launched by a Context, but it is helpful
     * for the Context to know if it has an active Session associated with it.
     *
     * \todo Use registration/deregistration protocol instead.
     *       Requires logging facility.
     */
    std::weak_ptr<Session> session_;

    /*!
     * \brief mdrun command line arguments.
     *
     * Store arguments provided by the client and pass them when launching
     * a simulation runner. This allows client code to access the same
     * options as are available to mdrun on the command line while the API
     * evolves.
     */
    MDArgs mdArgs_;

    /*!
     * \brief Legacy option-handling and set up for mdrun.
     *
     * This object should not exist, but is necessary now to introduce
     * the API in a way that means CLI and API work similarly and do not
     * duplicate definitions e.g. of command-line options.
     */
    gmx::LegacyMdrunOptions options_;

    /*!
     * \brief Scoped MPI management.
     *
     * Part of the ContextImpl invariant establishes a point where MPI initialization status is
     * known.
     *
     * To ensure the MpiContextManager is initialized only once, we use a const member that must
     * be initialized at construction.
     */
    const MpiContextManager mpi_;

private:
    /*!
     * \brief Basic constructor.
     *
     * Don't use this. Use create() to get a shared pointer right away.
     * Otherwise, shared_from_this() is potentially dangerous.
     */
    explicit ContextImpl() noexcept(std::is_nothrow_constructible_v<gmx::LegacyMdrunOptions>);
};


class CommHandle
{
public:
    using commType = MPI_Comm;
    MPI_Comm communicator{ MPI_COMM_NULL };
};

} // end namespace gmxapi
#endif // GMXAPI_CONTEXT_IMPL_H
