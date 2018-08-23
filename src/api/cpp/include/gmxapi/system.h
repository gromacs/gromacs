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
#ifndef GMXAPI_SYSTEM_H
#define GMXAPI_SYSTEM_H
/*! \file
 * \brief Declare container for molecular systems
 *
 * \ingroup gmxapi
 */
#include <memory>

#include "gmxapi/gmxapi.h"

namespace gmxapi
{

// Forward declaration for a return type defined elsewhere.
class Session;

/// Container for molecular model and simulation parameters.
/*!
 * \brief Deprecated: A wrapper for gmx::Mdrunner
 *
 * It was not intended to end up as such, but gmxapi::System has ended up as basically a
 * wrapper for the gmx::Mdrunner simulator object. As it is, this class does not fit the
 * current gmxapi paradigm and will be removed, reworked, or renamed soon.
 *
 * # Protocol
 *
 * As of gmxapi 0.0.6, a simulation is configured and launched as follows.
 *
 * 1. Caller gets a System handle with gmxapi::fromTprFile().
 * 2. Caller optionally attaches additional MD Modules with getSpec()->addModule(std::shared_ptr<gmxapi::MDModule> module). See gmxapi::MDHolder
 * 3. Caller gets a runnable object by passing a Context to System::launch()
 *
 * During launch() configured gmxapi::MDModules are attached to the simulator, which is
 * then run by calling run() on the object returned by launch().
 *
 * \ingroup gmxapi
 */
class System final
{
    public:
        /*! \brief Private implementation class.
         *
         * System::Impl does not have a public interface and is only exposed in opaque pointers.
         */
        class Impl;

        /// A blank system object is possible, but not yet useful.
        System();
        /// No copy.
        /*! The semantics of copying a System are ambiguous, so disallow implicit
         * copy. Some sort of prototype or clone idiom is probably useful, but
         * needs to explicitly identify any expensive operations.
         */
        System(const System &)            = delete;
        /// No copy.
        System              &operator=(const System &) = delete;

        /// Allow move.
        System(System &&) noexcept;
        /// Allow move.
        System &operator=(System &&) noexcept;

        /*!
         * \brief Create by taking ownership of an implementation object.
         *
         * \param implementation
         */
        explicit System(std::unique_ptr<Impl> &&implementation);

        /// \cond internal
        /// Destructor defined later to allow unique_ptr members of partially-defined types.
        ~System();
        /// \endcond

        std::shared_ptr<Session> launch(std::shared_ptr<Context> context);
        /// \}

        /*!
         * \brief Get the status of the last API call involving this system.
         *
         * \return copy of the most recent status.
         */
        Status status();

//        /// Get a handle to system atoms.
//        std::unique_ptr<Atoms> atoms();

    private:
        /*!
         * \brief Opaque pointer to implementation.
         */
        std::unique_ptr<Impl> impl_;
};


/// Defines an MD workflow from a TPR file.
/*! The TPR file has sufficient information to fully specify an MD run, though
 * various parameters are implicit until the work is launched. The TPR filename
 * provided must refer to identical TPR files at the API client and at the
 * master rank of the execution host.
 *
 * \param filename Filesystem path of TPR file.
 * \returns gmxapi::System object with the specified workflow.
 * \ingroup gmxapi
 */
std::unique_ptr<gmxapi::System> fromTprFile(std::string filename);

}      // end namespace gmxapi

#endif // include guard
