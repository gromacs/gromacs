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
//
// Created by Eric Irrgang on 11/9/17.
//

#ifndef GMXAPI_GROMACSFWD_H
#define GMXAPI_GROMACSFWD_H

/*! \ingroup gmxapi
 * \file
 * \brief Provide forward declarations for symbols in the GROMACS public headers.
 *
 * Basic API clients only need to compile
 * and link against the gmxapi target, but some gmxapi classes use opaque pointers to
 * library classes that are forward-declared here.
 *
 * We don't want to include ::gmx headers if we don't have to, but we need to declare
 * some things in the ::gmx namespace somewhere. These are forward declarations for
 * opaque pointers in libgromacs for client code building against libgmxapi.
 * Client code that is
 * more entwined with libgromacs can include headers from there.
 *
 * For maximal compatibility with other libgmxapi clients (such as third-party
 * Python modules), client code should use the wrappers and protocols in the
 * gmxapi.h header. Note that there is a separate CMake target to build the full
 * developer documentation for gmxapi.
 *
 * Refer to GMXAPI developer docs for the protocols that map gmxapi interfaces to
 * GROMACS library interfaces.
 * Refer to the GROMACS developer
 * documentation for details on library interfaces forward-declared in this header.
 *
 * \todo It would be nice to include links to the documentation for these classes, too.
 */

namespace gmx
{

class IRestraintPotential;
class TpxState;

}      // end namespace gmx

#endif //GMXAPI_GROMACSFWD_H
