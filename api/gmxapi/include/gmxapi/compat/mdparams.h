/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

/*! \file
 * \brief Compatibility header for simulation parameters.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#ifndef GMXAPICOMPAT_MDPARAMS_H
#define GMXAPICOMPAT_MDPARAMS_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gmxapi/gmxapi.h"
#include "gmxapi/gmxapicompat.h"

namespace gmxapicompat
{

// Forward declaration for private implementation class for GmxMdParams
class GmxMdParamsImpl;

class GmxMdParams
{
public:
    GmxMdParams();
    ~GmxMdParams();
    GmxMdParams(const GmxMdParams&) = delete;
    GmxMdParams& operator=(const GmxMdParams&) = delete;
    GmxMdParams(GmxMdParams&& /*unused*/) noexcept;
    GmxMdParams& operator=(GmxMdParams&& /*unused*/) noexcept;

    explicit GmxMdParams(std::unique_ptr<GmxMdParamsImpl>&& impl);

    std::unique_ptr<GmxMdParamsImpl> params_;
};

/*!
 * \brief Get the list of parameter key words.
 *
 * \param params molecular simulation parameters object reference.
 * \return A new vector of parameter names.
 *
 * \note The returned data is a copy. Modifying the return value has no affect on
 * the original object inspected.
 */
std::vector<std::string> keys(const GmxMdParams& params);

} // end namespace gmxapicompat

#endif // GMXAPICOMPAT_MDPARAMS_H
