/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "rf_util.h"

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

void calc_rffac(FILE* fplog, real eps_r, real eps_rf, real Rc, real* krf, real* crf)
{
    /* eps == 0 signals infinite dielectric */
    if (eps_rf == 0)
    {
        *krf = 1 / (2 * Rc * Rc * Rc);
    }
    else
    {
        *krf = (eps_rf - eps_r) / (2 * eps_rf + eps_r) / (Rc * Rc * Rc);
    }
    *crf = 1 / Rc + *krf * Rc * Rc;

    if (fplog)
    {
        fprintf(fplog,
                "%s:\n"
                "epsRF = %g, rc = %g, krf = %g, crf = %g, epsfac = %g\n",
                enumValueToString(CoulombInteractionType::RF),
                eps_rf,
                Rc,
                *krf,
                *crf,
                gmx::c_one4PiEps0 / eps_r);
        if (*krf > 0)
        {
            // Make sure we don't lose resolution in pow() by casting real arg to double
            real rmin = gmx::invcbrt(static_cast<double>(*krf * 2.0));
            fprintf(fplog, "The electrostatics potential has its minimum at r = %g\n", rmin);
        }
    }
}
