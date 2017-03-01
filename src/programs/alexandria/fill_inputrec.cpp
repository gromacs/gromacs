/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "fill_inputrec.h"

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/smalloc.h"

void fill_inputrec(t_inputrec *ir)
{
    ir->bAdress          = false;
    ir->cutoff_scheme    = ecutsGROUP;
    ir->tabext           = 0; /* nm */
    ir->ePBC             = epbcNONE;
    ir->ns_type          = ensSIMPLE;
    ir->epsilon_r        = 1;
    ir->vdwtype          = evdwCUT;
    ir->coulombtype      = eelCUT;
    ir->coulomb_modifier = eintmodNONE;
    ir->eDispCorr        = edispcNO;
    ir->vdw_modifier     = eintmodNONE;
    ir->niter            = 100;
    ir->em_stepsize      = 1e-2; // nm
    ir->em_tol           = 1e-2;
    ir->opts.ngener      = 1;
    snew(ir->fepvals, 1);
    snew(ir->opts.egp_flags, 1);
}
