/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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


/*! \libinternal
 * \defgroup module_fmm Fast Multipole Method (FMM)
 * \ingroup group_mdrun
 *
 * \brief
 * Allows to calculate electrostatics using the Fast Multipole Method (FMM).
 *
 * Relies on the fmsolvr git submodule (or, potentially another library that provides
 * FMM functionality) to be available in src/fmm/fmsolvr.
 *
 */

/*! \libinternal \file
 *
 * \brief
 * This file contains the necessary function declarations for
 * using the FMM to calculate electrostatic energies and forces.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_fmm
 */
#ifndef GMX_FMM_FMM_H
#define GMX_FMM_FMM_H

//#define GMX_WITH_FMM  // This controls whether or not the FMM code is compiled in

#include <memory>


namespace gmx
{

/*! \brief Identification string for the FMM options section in .mdp input.
 *
 * The same identifier string needs to be used in Fmm::initMdpTransform (fmm.cpp)
 * and makeModuleOptions() (mdmodules.cpp), therefore we use a global variable
 * defined here.
 */
extern const char *nameOfFmmOptionsSection; //= "fast-multipole-method";


class IMDModule;

/*! \brief
 * Creates a module allowing Coulomb force calculations with the FMM.
 *
 */
std::unique_ptr<IMDModule> createFastMultipoleModule();

}      // namespace gmx

#endif // GMX_FMM_FMM_H
