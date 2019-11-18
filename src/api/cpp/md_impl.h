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
#ifndef GMXAPI_MD_IMPL_H
#define GMXAPI_MD_IMPL_H
/*! \file
 * \brief Declarations for molecular dynamics API implementation details.
 *
 * \ingroup gmxapi
 */

#include <memory>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"

namespace gmxapi
{

class MDWorkSpec;

/*!
 * \brief Implementation class to hide guts of MDHolder
 *
 * Holds the gmxapi interface for an object that can help instantiate the gmx::MdRunner
 */
class MDHolder::Impl
{
public:
    /*!
     * \brief Construct by capturing a messaging object.
     *
     * \param spec operations specified for a workflow and the means to instantiate them.
     */
    explicit Impl(std::shared_ptr<MDWorkSpec>&& spec);

    /*!
     * \brief Shared ownership of the gmxapi object used for higher level message passing.
     */
    std::shared_ptr<MDWorkSpec> spec_{ nullptr };
};

} // namespace gmxapi

#endif // header guard
