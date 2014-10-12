/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */

#include "gmxpre.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"

static void reset_pmeonly_counters(gmx_wallcycle_t wcycle,
                                   gmx_walltime_accounting_t walltime_accounting,
                                   t_nrnb *nrnb, t_inputrec *ir,
                                   gmx_int64_t step)
{
    /* Reset all the counters related to performance over the run */
    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    init_nrnb(nrnb);
    if (ir->nsteps >= 0)
    {
        /* ir->nsteps is not used here, but we update it for consistency */
        ir->nsteps -= step - ir->init_step;
    }
    ir->init_step = step;
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_start(walltime_accounting);
}


static void gmx_pmeonly_switch(int *npmedata, struct gmx_pme_t ***pmedata,
                               ivec grid_size,
                               t_commrec *cr, t_inputrec *ir,
                               struct gmx_pme_t **pme_ret)
{
    int               ind;
    struct gmx_pme_t *pme = NULL;

    ind = 0;
    while (ind < *npmedata)
    {
        pme = (*pmedata)[ind];
        if (pme->nkx == grid_size[XX] &&
            pme->nky == grid_size[YY] &&
            pme->nkz == grid_size[ZZ])
        {
            *pme_ret = pme;

            return;
        }

        ind++;
    }

    (*npmedata)++;
    srenew(*pmedata, *npmedata);

    /* Generate a new PME data structure, copying part of the old pointers */
    gmx_pme_reinit(&((*pmedata)[ind]), cr, pme, ir, grid_size);

    *pme_ret = (*pmedata)[ind];
}

int gmx_pmeonly(struct gmx_pme_t *pme,
                t_commrec *cr,    t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                real ewaldcoeff_q, real ewaldcoeff_lj,
                t_inputrec *ir)
{
    int                npmedata;
    struct gmx_pme_t **pmedata;
    gmx_pme_pp_t       pme_pp;
    int                ret;
    int                natoms;
    matrix             box;
    rvec              *x_pp       = NULL, *f_pp = NULL;
    real              *chargeA    = NULL, *chargeB = NULL;
    real              *c6A        = NULL, *c6B = NULL;
    real              *sigmaA     = NULL, *sigmaB = NULL;
    real               lambda_q   = 0;
    real               lambda_lj  = 0;
    int                maxshift_x = 0, maxshift_y = 0;
    real               energy_q, energy_lj, dvdlambda_q, dvdlambda_lj;
    matrix             vir_q, vir_lj;
    float              cycles;
    int                count;
    gmx_bool           bEnerVir;
    int                pme_flags;
    gmx_int64_t        step;
    ivec               grid_switch;

    /* This data will only use with PME tuning, i.e. switching PME grids */
    npmedata = 1;
    snew(pmedata, npmedata);
    pmedata[0] = pme;

    pme_pp = gmx_pme_pp_init(cr);

    init_nrnb(mynrnb);

    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* The reason for having a loop here is PME grid tuning/switching */
        do
        {
            /* Domain decomposition */
            ret = gmx_pme_recv_coeffs_coords(pme_pp,
                                             &natoms,
                                             &chargeA, &chargeB,
                                             &c6A, &c6B,
                                             &sigmaA, &sigmaB,
                                             box, &x_pp, &f_pp,
                                             &maxshift_x, &maxshift_y,
                                             &pme->bFEP_q, &pme->bFEP_lj,
                                             &lambda_q, &lambda_lj,
                                             &bEnerVir,
                                             &pme_flags,
                                             &step,
                                             grid_switch, &ewaldcoeff_q, &ewaldcoeff_lj);

            if (ret == pmerecvqxSWITCHGRID)
            {
                /* Switch the PME grid to grid_switch */
                gmx_pmeonly_switch(&npmedata, &pmedata, grid_switch, cr, ir, &pme);
            }

            if (ret == pmerecvqxRESETCOUNTERS)
            {
                /* Reset the cycle and flop counters */
                reset_pmeonly_counters(wcycle, walltime_accounting, mynrnb, ir, step);
            }
        }
        while (ret == pmerecvqxSWITCHGRID || ret == pmerecvqxRESETCOUNTERS);

        if (ret == pmerecvqxFINISH)
        {
            /* We should stop: break out of the loop */
            break;
        }

        if (count == 0)
        {
            wallcycle_start(wcycle, ewcRUN);
            walltime_accounting_start(walltime_accounting);
        }

        wallcycle_start(wcycle, ewcPMEMESH);

        dvdlambda_q  = 0;
        dvdlambda_lj = 0;
        clear_mat(vir_q);
        clear_mat(vir_lj);

        gmx_pme_do(pme, 0, natoms, x_pp, f_pp,
                   chargeA, chargeB, c6A, c6B, sigmaA, sigmaB, box,
                   cr, maxshift_x, maxshift_y, mynrnb, wcycle,
                   vir_q, ewaldcoeff_q, vir_lj, ewaldcoeff_lj,
                   &energy_q, &energy_lj, lambda_q, lambda_lj, &dvdlambda_q, &dvdlambda_lj,
                   pme_flags | GMX_PME_DO_ALL_F | (bEnerVir ? GMX_PME_CALC_ENER_VIR : 0));

        cycles = wallcycle_stop(wcycle, ewcPMEMESH);

        gmx_pme_send_force_vir_ener(pme_pp,
                                    f_pp, vir_q, energy_q, vir_lj, energy_lj,
                                    dvdlambda_q, dvdlambda_lj, cycles);

        count++;
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);

    walltime_accounting_end(walltime_accounting);

    return 0;
}
