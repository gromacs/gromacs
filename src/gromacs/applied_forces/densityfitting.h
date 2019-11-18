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
#ifndef GMX_APPLIED_FORCES_DENSITYFITTING_H
#define GMX_APPLIED_FORCES_DENSITYFITTING_H

#include <memory>
#include <string>

namespace gmx
{

class IMDModule;
struct MdModulesNotifier;

/*! \libinternal \brief Information about the density fitting module.
 *
 * Provides name and method to create a density fitting module.
 */
struct DensityFittingModuleInfo
{
    /*! \brief
     * Creates a module for applying forces to fit a given density.
     *
     * Fitting an all-atom structure into an experimental cryo-EM density map is a
     * typical application.
     */
    static std::unique_ptr<IMDModule> create(MdModulesNotifier* notifier);
    //! The name of the module
    static const std::string name_;
};

} // namespace gmx

#endif
