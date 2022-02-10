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

#ifndef GMXAPI_MDMODULE_H
#define GMXAPI_MDMODULE_H
/*! \file
 * \brief Declare MDModule class for interacting with GROMACS library.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_md
 */

#include <memory>

#include "gmxapi/gromacsfwd.h"

namespace gmxapi
{

/*!
 * \brief Base class for computational components of MD containers.
 *
 * Instances of this class and its descendents provide member functions that return
 * GROMACS library objects that are not defined in this API, but which should be
 * defined in the public part of the GROMACS library API. Forward declarations of the library
 * classes are in gmxapi/gromacsfwd.h so that basic API clients only need to compile
 * and link against the gmxapi target. To extend the API may require GROMACS library
 * headers and possibly linking against `libgromacs`. Refer to the GROMACS developer
 * documentation for details.
 *
 * gmxapi::MDModule participates in the binding protocol described for gmxapi::Session.
 *
 * For portability, a shared_ptr to MDModule can be passed within the
 * gmxapi::MDHolder wrapper declared in gmxapi.h
 *
 * \ingroup gmxapi_md
 */
class MDModule
{
public:
    virtual ~MDModule();

    /*!
     * \brief Get a user-friendly identifier for an instance.
     *
     * \return see documentation for specific MDModule implementations.
     */
    virtual const char* name() const;

    /*!
     * \brief Allows module to provide a restraint implementation.
     *
     * To implement a restraint, override this function.
     * \return shared ownership of a restraint implementation or nullptr if not implemented.
     *
     * With future maturation, this interface will presumably be revised to something
     * more abstract, though I'm not sure what form that would take. We will probably
     * still need to have a general set of possible module types defined with the API,
     * in which case it does make sense to have clearly typed dispatching, and
     * `bool hasRestraint = module->getRestraint() != nullptr;` might be the simplest thing.
     *
     * Implementing a restraint is explained in the GROMACS developer documentation,
     * which is currently built separately from the GMXAPI documentation.
     * Also, refer to the sample plugin in a repository hosted in the same
     * place this git repository is found.
     */
    virtual std::shared_ptr<::gmx::IRestraintPotential> getRestraint();
};


} // end namespace gmxapi

#endif // GMXAPI_MDMODULE_H
