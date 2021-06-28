/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2021, by the GROMACS development team, led by
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
/*! \file
 * \brief Declarations for Context wrappers.
 *
 * \ingroup module_python
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#ifndef GMXPY_PYCONTEXT_H
#define GMXPY_PYCONTEXT_H

#include "pybind11/pybind11.h"

#include "gmxapi/context.h"
#include "gmxapi/md.h"

namespace gmxpy
{

using gmxapi::MDArgs;


/*!
 * \brief Wrapper for gmxapi::Context
 *
 * Proxies gmxapi::Context methods and includes additions not yet provided by
 * by upstream library.
 */
class PyContext
{
public:
    PyContext();
    void                             setMDArgs(const MDArgs& mdArgs);
    std::shared_ptr<gmxapi::Session> launch(const gmxapi::Workflow& work);
    std::shared_ptr<gmxapi::Context> get() const;

    void addMDModule(const pybind11::object& forceProvider) const;

    /*!
     * \brief Borrow shared ownership of the System's container of associated modules.
     *
     * Used with gmxapi::MDHolder to add MD Modules to the simulation to be run.
     *
     * \return handle to be passed to gmxapi::MDHolder
     *
     */
    std::shared_ptr<gmxapi::MDWorkSpec> getSpec() const;

private:
    std::shared_ptr<gmxapi::Context>    context_;
    std::shared_ptr<gmxapi::MDWorkSpec> workNodes_;
};


} // end namespace gmxpy

#endif // GMXPY_PYCONTEXT_H
