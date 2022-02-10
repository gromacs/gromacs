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
#ifndef GMXAPI_SYSTEM_IMPL_H
#define GMXAPI_SYSTEM_IMPL_H

/*! \file
 * \brief Declare implementation details for gmxapi::System.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <string>

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
    ~Impl();

    Impl(Impl&& /*unused*/) noexcept;
    Impl& operator=(Impl&& source) noexcept;
    /*! \endcond */

    /*!
     * \brief Initialize from a work description.
     *
     * \param workflow Simulation work to perform.
     */
    explicit Impl(std::unique_ptr<gmxapi::Workflow> workflow) noexcept;

    /*!
     * \brief Launch the configured simulation.
     *
     * \param context Runtime execution context in which to run simulation.
     * \return Ownership of a new simulation session.
     *
     * The session is returned as a shared pointer so that the Context can
     * maintain a weak reference to it via std::weak_ptr.
     */
    std::shared_ptr<Session> launch(const std::shared_ptr<Context>& context);

    //! Description of simulation work.
    std::shared_ptr<Workflow> workflow_;

    /*!
     * \brief Specified simulation work.
     *
     * \todo merge Workflow and MDWorkSpec
     */
    std::shared_ptr<gmxapi::MDWorkSpec> spec_;
};

} // end namespace gmxapi

#endif // header guard
