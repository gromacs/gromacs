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
/*! \internal
 * \brief Defines the dispatch function for the .mdp integrator field.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "legacysimulator.h"

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{
//! \brief Run the correct integrator function.
void LegacySimulator::run()
{
    switch (inputrec->eI)
    {
        case eiMD:
        case eiBD:
        case eiSD1:
        case eiVV:
        case eiVVAK:
            if (!EI_DYNAMICS(inputrec->eI))
            {
                GMX_THROW(APIError(
                        "do_md integrator would be called for a non-dynamical integrator"));
            }
            if (doRerun)
            {
                do_rerun();
            }
            else
            {
                do_md();
            }
            break;
        case eiMimic:
            if (doRerun)
            {
                do_rerun();
            }
            else
            {
                do_mimic();
            }
            break;
        case eiSteep: do_steep(); break;
        case eiCG: do_cg(); break;
        case eiNM: do_nm(); break;
        case eiLBFGS: do_lbfgs(); break;
        case eiTPI:
        case eiTPIC:
            if (!EI_TPI(inputrec->eI))
            {
                GMX_THROW(APIError("do_tpi integrator would be called for a non-TPI integrator"));
            }
            do_tpi();
            break;
        case eiSD2_REMOVED: GMX_THROW(NotImplementedError("SD2 integrator has been removed"));
        default: GMX_THROW(APIError("Non existing integrator selected"));
    }
}

} // namespace gmx
