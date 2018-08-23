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
#ifndef GMXAPI_SYSTEM_IMPL_H
#define GMXAPI_SYSTEM_IMPL_H

/*! \file
 * \brief Declare implementation details for gmxapi::System.
 *
 * \ingroup gmxapi
 */

#include <string>

#include "gmxapi/status.h"
#include "gmxapi/system.h"

namespace gmxapi
{

class Context;
class Workflow;

/*!
 * \brief Private implementation for gmxapi::System
 *
 * \ingroup gmxapi
 */
class System::Impl final
{
    public:
        /*! \cond */
        Impl();
        ~Impl();

        Impl(Impl &&) noexcept            = default;
        Impl &operator=(Impl &&) noexcept = default;
        /*! \endcond */

        /*!
         * \brief Initialize from a work description.
         *
         * \param workflow Simulation work to perform.
         */
        explicit Impl(std::unique_ptr<gmxapi::Workflow> &&workflow) noexcept;

        /*!
         * \brief Get the status of the last operation.
         *
         * Force resolution of any pending operations and return the status to
         * the client.
         *
         * \return success if the last operation on the system completed without problems.
         */
        Status status() const;

        std::shared_ptr<Session> launch(std::shared_ptr<Context> context);

        Status setRestraint(std::shared_ptr<gmxapi::MDModule> module);
        std::shared_ptr<MDWorkSpec> getSpec();

    private:
        //! Execution context.
        std::shared_ptr<Context>            context_;
        //! Description of simulation work.
        std::shared_ptr<Workflow>           workflow_;
        // \todo merge Workflow and MDWorkSpec
        std::shared_ptr<gmxapi::MDWorkSpec> spec_;
        //! Cached Status object.
        std::unique_ptr<Status>             status_;
};

}      // end namespace gmxapi

#endif // header guard
