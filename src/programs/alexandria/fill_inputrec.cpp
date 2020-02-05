/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "gmxpre.h"

#include "fill_inputrec.h"

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/keyvaluetree.h"
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
    ir->em_tol           = 1e-4;
    ir->opts.ngener      = 1;
    snew(ir->fepvals, 1);
    snew(ir->opts.egp_flags, 1);
    ir->params           = new gmx::KeyValueTreeObject();
}
