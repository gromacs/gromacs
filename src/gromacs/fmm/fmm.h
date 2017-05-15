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
 * Prepares \Gromacs to use a Fast Multipole Method (FMM) to evaluate Coulomb interactions.
 *
 * The code in fmm.h and fmm.cpp handles the .mdp input parameter parsing (like depth and
 * multipole order). An underlying FMM can insert itself into the calculateForces() method
 * in fmm-impl.h
 */

/*! \libinternal \file
 *
 * \brief
 * This file provides basic infrastructure allowing \Gromacs to interface to a FMM.
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

#include "gromacs/utility/real.h"

namespace gmx
{

/*! \internal
 * \brief Input parameters needed for FMM methods, defined in the .mdp file
 */
struct fmmInputParameters
{
    int  depth;      /**< tree depth                                                  */
    int  order;      /**< order of the multipole expansion                            */
    int  separation; /**< separation criterion                                        */
    real precision;  /**< Maximum accepted relative error in Coulomb potential energy */

    // default values of FMM input parameters:
    fmmInputParameters() : depth(-1), order(-1), separation(-1), precision(0.001) { }
};

class IMDModule;

/*! \brief
 * Creates a module allowing \Gromacs to link to a FMM implementation.
 *
 */
std::unique_ptr<IMDModule> createFastMultipoleModule();

}      // namespace gmx

#endif // GMX_FMM_FMM_H
