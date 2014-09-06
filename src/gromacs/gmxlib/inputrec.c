/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/legacyheaders/inputrec.h"

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/utility/fatalerror.h"


/* The minimum number of integration steps required for reasonably accurate
 * integration of first and second order coupling algorithms.
 */
const int nstmin_berendsen_tcouple =  5;
const int nstmin_berendsen_pcouple = 10;
const int nstmin_harmonic          = 20;

static int nst_wanted(const t_inputrec *ir)
{
    if (ir->nstlist > 0)
    {
        return ir->nstlist;
    }
    else
    {
        return 10;
    }
}

int ir_optimal_nstcalcenergy(const t_inputrec *ir)
{
    return nst_wanted(ir);
}

int tcouple_min_integration_steps(int etc)
{
    int n;

    switch (etc)
    {
        case etcNO:
            n = 0;
            break;
        case etcBERENDSEN:
        case etcYES:
            n = nstmin_berendsen_tcouple;
            break;
        case etcVRESCALE:
            /* V-rescale supports instantaneous rescaling */
            n = 0;
            break;
        case etcNOSEHOOVER:
            n = nstmin_harmonic;
            break;
        case etcANDERSEN:
        case etcANDERSENMASSIVE:
            n = 1;
            break;
        default:
            gmx_incons("Unknown etc value");
            n = 0;
    }

    return n;
}

int ir_optimal_nsttcouple(const t_inputrec *ir)
{
    int  nmin, nwanted, n;
    real tau_min;
    int  g;

    nmin = tcouple_min_integration_steps(ir->etc);

    nwanted = nst_wanted(ir);

    tau_min = 1e20;
    if (ir->etc != etcNO)
    {
        for (g = 0; g < ir->opts.ngtc; g++)
        {
            if (ir->opts.tau_t[g] > 0)
            {
                tau_min = min(tau_min, ir->opts.tau_t[g]);
            }
        }
    }

    if (nmin == 0 || ir->delta_t*nwanted <= tau_min)
    {
        n = nwanted;
    }
    else
    {
        n = (int)(tau_min/(ir->delta_t*nmin) + 0.001);
        if (n < 1)
        {
            n = 1;
        }
        while (nwanted % n != 0)
        {
            n--;
        }
    }

    return n;
}

int pcouple_min_integration_steps(int epc)
{
    int n;

    switch (epc)
    {
        case epcNO:
            n = 0;
            break;
        case etcBERENDSEN:
        case epcISOTROPIC:
            n = nstmin_berendsen_pcouple;
            break;
        case epcPARRINELLORAHMAN:
        case epcMTTK:
            n = nstmin_harmonic;
            break;
        default:
            gmx_incons("Unknown epc value");
            n = 0;
    }

    return n;
}

int ir_optimal_nstpcouple(const t_inputrec *ir)
{
    int  nmin, nwanted, n;

    nmin = pcouple_min_integration_steps(ir->epc);

    nwanted = nst_wanted(ir);

    if (nmin == 0 || ir->delta_t*nwanted <= ir->tau_p)
    {
        n = nwanted;
    }
    else
    {
        n = (int)(ir->tau_p/(ir->delta_t*nmin) + 0.001);
        if (n < 1)
        {
            n = 1;
        }
        while (nwanted % n != 0)
        {
            n--;
        }
    }

    return n;
}

gmx_bool ir_coulomb_switched(const t_inputrec *ir)
{
    return (ir->coulombtype == eelSWITCH ||
            ir->coulombtype == eelSHIFT ||
            ir->coulombtype == eelENCADSHIFT ||
            ir->coulombtype == eelPMESWITCH ||
            ir->coulombtype == eelPMEUSERSWITCH ||
            ir->coulomb_modifier == eintmodPOTSWITCH ||
            ir->coulomb_modifier == eintmodFORCESWITCH);
}

gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir->cutoff_scheme == ecutsVERLET ||
            ir_coulomb_switched(ir) || ir->coulomb_modifier != eintmodNONE ||
            ir->coulombtype == eelRF_ZERO);
}

gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir_coulomb_is_zero_at_cutoff(ir) || ir->coulombtype == eelUSER || ir->coulombtype == eelPMEUSER);
}

gmx_bool ir_vdw_switched(const t_inputrec *ir)
{
    return (ir->vdwtype == evdwSWITCH ||
            ir->vdwtype == evdwSHIFT ||
            ir->vdwtype == evdwENCADSHIFT ||
            ir->vdw_modifier == eintmodPOTSWITCH ||
            ir->vdw_modifier == eintmodFORCESWITCH);
}

gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir->cutoff_scheme == ecutsVERLET ||
            ir_vdw_switched(ir) || ir->vdw_modifier != eintmodNONE);
}

gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir_vdw_is_zero_at_cutoff(ir) || ir->vdwtype == evdwUSER);
}
