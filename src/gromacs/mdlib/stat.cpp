/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include <stdio.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/rbin.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/vcm.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct gmx_global_stat
{
    t_bin *rb;
    int   *itc0;
    int   *itc1;
} t_gmx_global_stat;

gmx_global_stat_t global_stat_init(t_inputrec *ir)
{
    gmx_global_stat_t gs;

    snew(gs, 1);

    gs->rb = mk_bin();
    snew(gs->itc0, ir->opts.ngtc);
    snew(gs->itc1, ir->opts.ngtc);

    return gs;
}

void global_stat_destroy(gmx_global_stat_t gs)
{
    destroy_bin(gs->rb);
    sfree(gs->itc0);
    sfree(gs->itc1);
    sfree(gs);
}

static int filter_enerdterm(real *afrom, gmx_bool bToBuffer, real *ato,
                            gmx_bool bTemp, gmx_bool bPres, gmx_bool bEner)
{
    int i, to, from;

    from = 0;
    to   = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (bToBuffer)
        {
            from = i;
        }
        else
        {
            to = i;
        }
        switch (i)
        {
            case F_EKIN:
            case F_TEMP:
            case F_DKDL:
                if (bTemp)
                {
                    ato[to++] = afrom[from++];
                }
                break;
            case F_PRES:
            case F_PDISPCORR:
                if (bPres)
                {
                    ato[to++] = afrom[from++];
                }
                break;
            default:
                if (bEner)
                {
                    ato[to++] = afrom[from++];
                }
                break;
        }
    }

    return to;
}

void global_stat(FILE *fplog, gmx_global_stat_t gs,
                 t_commrec *cr, gmx_enerdata_t *enerd,
                 tensor fvir, tensor svir, rvec mu_tot,
                 t_inputrec *inputrec,
                 gmx_ekindata_t *ekind, gmx_constr_t constr,
                 t_vcm *vcm,
                 int nsig, real *sig,
                 gmx_mtop_t *top_global, t_state *state_local,
                 gmx_bool bSumEkinhOld, int flags)
/* instead of current system, gmx_booleans for summing virial, kinetic energy, and other terms */
{
    t_bin     *rb;
    int       *itc0, *itc1;
    int        ie    = 0, ifv = 0, isv = 0, irmsd = 0, imu = 0;
    int        idedl = 0, idvdll = 0, idvdlnl = 0, iepl = 0, icm = 0, imass = 0, ica = 0, inb = 0;
    int        isig  = -1;
    int        icj   = -1, ici = -1, icx = -1;
    int        inn[egNR];
    real       copyenerd[F_NRE];
    int        nener, j;
    real      *rmsd_data = NULL;
    double     nb;
    gmx_bool   bVV, bTemp, bEner, bPres, bConstrVir, bEkinAveVel, bReadEkin;

    bVV           = EI_VV(inputrec->eI);
    bTemp         = flags & CGLO_TEMPERATURE;
    bEner         = flags & CGLO_ENERGY;
    bPres         = (flags & CGLO_PRESSURE);
    bConstrVir    = (flags & CGLO_CONSTRAINT);
    bEkinAveVel   = (inputrec->eI == eiVV || (inputrec->eI == eiVVAK && bPres));
    bReadEkin     = (flags & CGLO_READEKIN);

    rb   = gs->rb;
    itc0 = gs->itc0;
    itc1 = gs->itc1;


    reset_bin(rb);
    /* This routine copies all the data to be summed to one big buffer
     * using the t_bin struct.
     */

    /* First, we neeed to identify which enerd->term should be
       communicated.  Temperature and pressure terms should only be
       communicated and summed when they need to be, to avoid repeating
       the sums and overcounting. */

    nener = filter_enerdterm(enerd->term, TRUE, copyenerd, bTemp, bPres, bEner);

    /* First, the data that needs to be communicated with velocity verlet every time
       This is just the constraint virial.*/
    if (bConstrVir)
    {
        isv = add_binr(rb, DIM*DIM, svir[0]);
        where();
    }

/* We need the force virial and the kinetic energy for the first time through with velocity verlet */
    if (bTemp || !bVV)
    {
        if (ekind)
        {
            for (j = 0; (j < inputrec->opts.ngtc); j++)
            {
                if (bSumEkinhOld)
                {
                    itc0[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinh_old[0]);
                }
                if (bEkinAveVel && !bReadEkin)
                {
                    itc1[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinf[0]);
                }
                else if (!bReadEkin)
                {
                    itc1[j] = add_binr(rb, DIM*DIM, ekind->tcstat[j].ekinh[0]);
                }
            }
            /* these probably need to be put into one of these categories */
            where();
            idedl = add_binr(rb, 1, &(ekind->dekindl));
            where();
            ica   = add_binr(rb, 1, &(ekind->cosacc.mvcos));
            where();
        }
    }
    where();

    if (bPres || !bVV)
    {
        ifv = add_binr(rb, DIM*DIM, fvir[0]);
    }


    if (bEner)
    {
        where();
        ie  = add_binr(rb, nener, copyenerd);
        where();
        if (constr)
        {
            rmsd_data = constr_rmsd_data(constr);
            if (rmsd_data)
            {
                irmsd = add_binr(rb, inputrec->eI == eiSD2 ? 3 : 2, rmsd_data);
            }
        }
        if (!NEED_MUTOT(*inputrec))
        {
            imu = add_binr(rb, DIM, mu_tot);
            where();
        }

        for (j = 0; (j < egNR); j++)
        {
            inn[j] = add_binr(rb, enerd->grpp.nener, enerd->grpp.ener[j]);
        }
        where();
        if (inputrec->efep != efepNO)
        {
            idvdll  = add_bind(rb, efptNR, enerd->dvdl_lin);
            idvdlnl = add_bind(rb, efptNR, enerd->dvdl_nonlin);
            if (enerd->n_lambda > 0)
            {
                iepl = add_bind(rb, enerd->n_lambda, enerd->enerpart_lambda);
            }
        }
    }

    if (vcm)
    {
        icm   = add_binr(rb, DIM*vcm->nr, vcm->group_p[0]);
        where();
        imass = add_binr(rb, vcm->nr, vcm->group_mass);
        where();
        if (vcm->mode == ecmANGULAR)
        {
            icj   = add_binr(rb, DIM*vcm->nr, vcm->group_j[0]);
            where();
            icx   = add_binr(rb, DIM*vcm->nr, vcm->group_x[0]);
            where();
            ici   = add_binr(rb, DIM*DIM*vcm->nr, vcm->group_i[0][0]);
            where();
        }
    }

    if (DOMAINDECOMP(cr))
    {
        nb  = cr->dd->nbonded_local;
        inb = add_bind(rb, 1, &nb);
    }
    where();
    if (nsig > 0)
    {
        isig = add_binr(rb, nsig, sig);
    }

    /* Global sum it all */
    if (debug)
    {
        fprintf(debug, "Summing %d energies\n", rb->maxreal);
    }
    sum_bin(rb, cr);
    where();

    /* Extract all the data locally */

    if (bConstrVir)
    {
        extract_binr(rb, isv, DIM*DIM, svir[0]);
    }

    /* We need the force virial and the kinetic energy for the first time through with velocity verlet */
    if (bTemp || !bVV)
    {
        if (ekind)
        {
            for (j = 0; (j < inputrec->opts.ngtc); j++)
            {
                if (bSumEkinhOld)
                {
                    extract_binr(rb, itc0[j], DIM*DIM, ekind->tcstat[j].ekinh_old[0]);
                }
                if (bEkinAveVel && !bReadEkin)
                {
                    extract_binr(rb, itc1[j], DIM*DIM, ekind->tcstat[j].ekinf[0]);
                }
                else if (!bReadEkin)
                {
                    extract_binr(rb, itc1[j], DIM*DIM, ekind->tcstat[j].ekinh[0]);
                }
            }
            extract_binr(rb, idedl, 1, &(ekind->dekindl));
            extract_binr(rb, ica, 1, &(ekind->cosacc.mvcos));
            where();
        }
    }
    if (bPres || !bVV)
    {
        extract_binr(rb, ifv, DIM*DIM, fvir[0]);
    }

    if (bEner)
    {
        extract_binr(rb, ie, nener, copyenerd);
        if (rmsd_data)
        {
            extract_binr(rb, irmsd, inputrec->eI == eiSD2 ? 3 : 2, rmsd_data);
        }
        if (!NEED_MUTOT(*inputrec))
        {
            extract_binr(rb, imu, DIM, mu_tot);
        }

        for (j = 0; (j < egNR); j++)
        {
            extract_binr(rb, inn[j], enerd->grpp.nener, enerd->grpp.ener[j]);
        }
        if (inputrec->efep != efepNO)
        {
            extract_bind(rb, idvdll, efptNR, enerd->dvdl_lin);
            extract_bind(rb, idvdlnl, efptNR, enerd->dvdl_nonlin);
            if (enerd->n_lambda > 0)
            {
                extract_bind(rb, iepl, enerd->n_lambda, enerd->enerpart_lambda);
            }
        }
        if (DOMAINDECOMP(cr))
        {
            extract_bind(rb, inb, 1, &nb);
            if ((int)(nb + 0.5) != cr->dd->nbonded_global)
            {
                dd_print_missing_interactions(fplog, cr, (int)(nb + 0.5), top_global, state_local);
            }
        }
        where();

        filter_enerdterm(copyenerd, FALSE, enerd->term, bTemp, bPres, bEner);
    }

    if (vcm)
    {
        extract_binr(rb, icm, DIM*vcm->nr, vcm->group_p[0]);
        where();
        extract_binr(rb, imass, vcm->nr, vcm->group_mass);
        where();
        if (vcm->mode == ecmANGULAR)
        {
            extract_binr(rb, icj, DIM*vcm->nr, vcm->group_j[0]);
            where();
            extract_binr(rb, icx, DIM*vcm->nr, vcm->group_x[0]);
            where();
            extract_binr(rb, ici, DIM*DIM*vcm->nr, vcm->group_i[0][0]);
            where();
        }
    }

    if (nsig > 0)
    {
        extract_binr(rb, isig, nsig, sig);
    }
    where();
}

int do_per_step(gmx_int64_t step, gmx_int64_t nstep)
{
    if (nstep != 0)
    {
        return ((step % nstep) == 0);
    }
    else
    {
        return 0;
    }
}
