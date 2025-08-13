/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

/*! \libinternal \file
 * \brief Defines factory structure for FMM MDModule class.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#ifndef GMX_FMM_MDMODULE_H
#define GMX_FMM_MDMODULE_H

#include <memory>
#include <string_view>

namespace gmx
{

class IMDModule;

/*! \internal
 *  \brief Information about the FMM module.
 *
 * Provides name and method to create a FMM module.
 */
struct FmmModuleInfo
{
    /*! \brief
     * Creates a module for computing forces using FMM.
     *
     */
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static constexpr std::string_view sc_name = "fmm";
};

} // namespace gmx

#endif // GMX_FMM_MDMODULE_H
