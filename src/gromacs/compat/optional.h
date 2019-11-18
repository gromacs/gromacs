/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Provides C++14-compatible implementation of std::optional.
 *
 * This imports the implementation found in src/external from
 * https://github.com/martinmoene/optional-lite.git.
 *
 * There is no Doxygen for this code, but it is intended to conform to
 * that of std::optional, so look in the usual C++17 documentation for
 * std::optional for that.
 *
 * \todo Remove when requiring C++17, which has a standardized version
 * of std::optional.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_compat
 * \inlibraryapi
 */
#ifndef GMX_COMPAT_OPTIONAL_H
#define GMX_COMPAT_OPTIONAL_H

#include <nonstd/optional.hpp>

namespace gmx
{
namespace compat
{

using nonstd::bad_optional_access;
using nonstd::in_place;
using nonstd::make_optional;
using nonstd::nullopt;
using nonstd::optional;

} // namespace compat
} // namespace gmx

#endif
