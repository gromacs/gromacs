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

#ifndef GMXAPI_RESOURCEASSIGNMENT_H
#define GMXAPI_RESOURCEASSIGNMENT_H

/*! \file
 * \brief Provide build-specific overloads for client-MPI-dependent stuff.
 *
 * Define the interfaces that a client must implement to generate MpiContextManager for the library.
 * Client code should use the gmxapi_mpi.h template header to generate code supporting
 * the required interface. (Client code should not need to include this header directly.)
 *
 * \note This is a generated header. Some definitions are determined when the GROMACS library
 * build is configured.
 *
 * If the library is built with tMPI, CommHandle is empty and offerComm is not defined.
 *
 * If the library is built with an MPI library, CommHandle holds a native MPI_Comm object and
 * the library-provided offerComm is used by the client MultiProcessingResources implementation
 * to support passing a copy of the communicator.
 *
 * Note: The communicator is not `dup`ed in this call. If the library later duplicates
 * the offered communicator, the library will be responsible for freeing the duplicate.
 * However, the caller is responsible for keeping any MPI environment valid while the library is
 * in use. For details, \see MpiContextManager.
 *
 * \author "M. Eric Irrgang <ericirrgang@gmail.com>"
 */

#include <functional>
#include <memory>

// The interface in this header is determined when the GROMACS library build is configured.
#cmakedefine01 GMX_LIB_MPI
#if GMX_LIB_MPI
#    include <mpi.h>
#endif

namespace gmxapi
{

// Forward declaration for library resources.
// CommHandle is opaque to the public API, but allows extensions to implement
// required library interfaces by composing the behavior of helpers like offerComm().
// The abstraction is important because the library implementations depend on
// the options with which GROMACS was built.
class CommHandle;

/*!
 * \brief Convey the resources that the client has directed the library to use within a Context.
 *
 * The object itself does not convey ownership of resources. However, the interfaces
 * for declaring assigned resources must have well-specified (and documented)
 * semantics for resource ownership. See assignResource().
 * (Note: the initial implementation only allows for assignment of an MPI Communicator.)
 *
 * The library and client share this interface definition. The implementation is
 * provided by client code with the aid of a template header.
 * Client developers should refer instead to gmxapi_mpi.h.
 */
class ResourceAssignment
{
public:
    virtual ~ResourceAssignment();
    [[nodiscard]] virtual int size() const = 0;
    [[nodiscard]] virtual int rank() const = 0;
    virtual void              applyCommunicator(class CommHandle* dst) const;
};

/*! \brief Offer the client communicator to the library.
 *
 * Helper function allowing clients to provide the MPI communicator for the library.
 * \param src client-provided communicator.
 * \param dst library recipient, abstracted to hide library MPI type details.
 *
 */
template<typename CommT>
void offerComm([[maybe_unused]] CommT src, [[maybe_unused]] CommHandle* dst)
{
    // Default overload does nothing.
}

#if GMX_LIB_MPI
/*! \brief Offer the client communicator to the library.
 *
 * Helper function allowing clients to provide communicator to library.
 * \param src communicator offered by client
 * \param dst opaque pointer to library resource destination
 *
 * This function is only available in MPI-enabled GROMACS.
 * Used indirectly in template helpers for client code to implement library interfaces.
 * Clients should use the assignResource() higher level function in explicit code,
 * which will evaluable appropriately for the target GROMACS library.
 *
 * \todo Provide a stub for docs even when not available.
 */
void offerComm(MPI_Comm src, CommHandle* dst);
#endif

} // end namespace gmxapi

#endif // GMXAPI_RESOURCEASSIGNMENT_H
