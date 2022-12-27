/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#ifndef GMXAPI_MPI_H
#define GMXAPI_MPI_H

#include <mpi.h>

#include <cassert>

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/mpi/resourceassignment.h"

/*! \file
 * \brief Provide details of any MPI implementation used when building the library.
 *
 * If the gmxapi library was built with MPI-enabled GROMACS, client software
 * using a compatible MPI implementation may use this header to access additional
 * functionality.
 *
 * Client software should use the CMake infrastructure to ensure a compatible
 * MPI implementation. Use `find_package(MPI COMPONENTS ...)` and
 * `target_link_libraries(... MPI::MPI_CXX)`
 *
 * Use of an incompatible MPI implementation will produce linking errors.
 */

namespace gmxapi
{

/*!
 * \brief Convey resources assigned by the client to the new library Context.
 *
 * Only enabled if the required MPI library calls are available.
 */
template<typename CommT, decltype(MPI_Comm_rank(CommT(), nullptr)) = 0, decltype(MPI_Comm_size(CommT(), nullptr)) = 0>
class ResourceAssignmentImpl : public ResourceAssignment
{
public:
    ~ResourceAssignmentImpl() override = default;

    /*!
     * \brief Initialize from borrowed communicator.
     *
     * \param communicator Client-provided communicator relayed through assignResource()
     *
     * Create an abstract wrapper for client-provided values with which to initialize
     * simulation resources. When this wrapper is used, the client is responsible for providing
     * a valid communicator that will remain valid for the life of the consumer.
     *
     * \throws UsageError if client has not provided a valid communicator.
     */
    explicit ResourceAssignmentImpl(const CommT& communicator) : communicator_{ communicator }
    {
        if (communicator_ == MPI_COMM_NULL)
        {
            throw UsageError("Cannot assign a Null communicator.");
        }
        int flag = 0;
        MPI_Initialized(&flag);
        if (!flag)
        {
            throw UsageError("Client has not initialized MPI context.");
        }
    }

    [[nodiscard]] int size() const override
    {
        assert(communicator_ != MPI_COMM_NULL && "Class invariant broken: invalid communicator.");
        int size = 0;
        MPI_Comm_size(communicator_, &size);
        return size;
    }

    [[nodiscard]] int rank() const override
    {
        assert(communicator_ != MPI_COMM_NULL && "Class invariant broken: invalid communicator.");
        // This default value will never be read, but the compiler can't tell
        // that it is initialized by the MPI call.
        int rank = -1;
        MPI_Comm_rank(communicator_, &rank);
        return rank;
    }

    void applyCommunicator(CommHandle* dst) const override
    {
        assert(communicator_ != MPI_COMM_NULL && "Class invariant broken: invalid communicator.");
        offerComm(communicator_, dst);
    }

    CommT communicator_;
};

/*!
 * \brief Template header utility for connecting to MPI implementations.
 *
 * The client provides a communicator for work executed within the scope of the Context.
 * Client remains responsible for freeing the communicator and finalizing the MPI environment.
 *
 * To use this helper function, the client software build environment must be
 * configured for an MPI implementation compatible with the target GROMACS library.
 *
 * If the library and client software are built with compatible MPI implementations,
 * the template will capture the communicator in an opaque container that can be
 * passed to the ContextImpl creation function. Otherwise, the opaque container
 * will carry a suitable NULL object. The template definition is necessarily
 * provided by a generated template header, since the definition depends on the
 * GROMACS installation.
 *
 * \tparam CommT Deduced type of client-provided communicator.
 * \param communicator MPI (sub)communicator linking collaborating processes in the new Context.
 * \return Opaque resource handle usable when setting up execution context.
 *
 * \see createContext(std::unique_ptr<ResourceAssignment> resources)
 *
 * The communicator resource type is a template parameter because MPI_Comm is commonly a C
 * typedef that varies between implementations, so we do not want to couple our
 * API to it, but we cannot forward-declare it.
 *
 * See also https://gitlab.com/gromacs/gromacs/-/issues/3650
 *
 * \throws UsageError if the provided resource is not usable.
 */
template<typename CommT>
std::unique_ptr<ResourceAssignment> assignResource(CommT communicator)
{
    if (communicator == MPI_COMM_NULL)
    {
        throw UsageError("Cannot assign a Null communicator.");
    }
    return std::make_unique<ResourceAssignmentImpl<CommT>>(communicator);
}

// Also ref https://stackoverflow.com/questions/49259704/pybind11-possible-to-use-mpi4py

} // end namespace gmxapi

#endif // GMXAPI_MPI_H
