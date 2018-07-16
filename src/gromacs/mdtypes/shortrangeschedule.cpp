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
/*! \internal \file
 * \brief
 * Implements short range schedules declared in
 * iforceprovider.h.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \ingroup module_mdtypes
 */

#include "gmxpre.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/math/vec.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/mdlib/calcmu.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/topology/topology.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/nbnxn_grid.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "forcerec.h"
#include "iforceschedule.h"
#include "state.h"
#include "mdatom.h"
#include "forceoutput.h"

using namespace gmx;

void SRCPUSingleNodeSchedule::computeStep(int flags)
{
    int                 cg1, i, j;
    double              mu[2*DIM];
    gmx_bool            bStateChanged, bNS, bFillGrid, bCalcCGCM;
    gmx_bool            bDoForces, bEmulGPU;
    rvec                vzero, box_diag;

    nonbonded_verlet_t *nbv = fr_->nbv;

    bStateChanged = ((flags & GMX_FORCE_STATECHANGED) != 0);
    bNS           = ((flags & GMX_FORCE_NS) != 0);
    bFillGrid     = (bNS && bStateChanged);
    bCalcCGCM     = bFillGrid;
    bDoForces     = ((flags & GMX_FORCE_FORCES) != 0);
    bEmulGPU      = (fr_->nbv->emulateGpu == EmulateGpuNonbonded::Yes);

    const int  start  = 0;
    const int  homenr = mdatoms_->homenr;

    int64_t    step = *(sp_->step_);

    auto       lambda = state_->lambda.data();

    matrix   * box = &(state_->box);

    history_t* hist = &(state_->hist);

    clear_mat(*vir_force_);

    gmx::ArrayRef<gmx::RVec> x = state_->x;

    auto                    *t = sp_->t_;

    cg1 = top_->cgs.nr;

    if (fr_->n_tpi > 0)
    {
        cg1--;
    }

    if (bStateChanged)
    {
        update_forcerec(fr_, state_->box);

        if (inputrecNeedMutot(inputrec_))
        {
            /* Calculate total (local) dipole moment in a temporary common array.
             * This makes it possible to sum them over nodes faster.
             */
            calc_mu(start, homenr,
                    x, mdatoms_->chargeA, mdatoms_->chargeB, mdatoms_->nChargePerturbed,
                    mu, mu+DIM);
        }
    }

    if (fr_->ePBC != epbcNONE)
    {
        /* Compute shift vectors every step,
         * because of pressure coupling or box deformation!
         */
        if ((flags & GMX_FORCE_DYNAMICBOX) && bStateChanged)
        {
            calc_shifts(*box, fr_->shift_vec);
        }

        if (bCalcCGCM)
        {
            put_atoms_in_box_omp(fr_->ePBC, *box, x.subArray(0, homenr));
            inc_nrnb(nrnb_, eNR_SHIFTX, homenr);
        }
        else if (EI_ENERGY_MINIMIZATION(inputrec_->eI) && graph_)
        {
            unshift_self(graph_, *box, as_rvec_array(x.data()));
        }
    }

    nbnxn_atomdata_copy_shiftvec((flags & GMX_FORCE_DYNAMICBOX) != 0,
                                 fr_->shift_vec, nbv->nbat);

    /* do gridding for pair search */
    if (bNS)
    {
        if (graph_ && bStateChanged)
        {
            /* Calculate intramolecular shift vectors to make molecules whole */
            mk_mshift(log_, graph_, fr_->ePBC, *box, as_rvec_array(x.data()));
        }

        clear_rvec(vzero);

        box_diag[XX] = (*box)[XX][XX];
        box_diag[YY] = (*box)[YY][YY];
        box_diag[ZZ] = (*box)[ZZ][ZZ];

        wallcycle_start(wcycle_, ewcNS);

        wallcycle_sub_start(wcycle_, ewcsNBS_GRID_LOCAL);
        nbnxn_put_on_grid(nbv->nbs.get(), fr_->ePBC, *box,
                          0, vzero, box_diag,
                          nullptr, 0, mdatoms_->homenr, -1,
                          fr_->cginfo, x,
                          0, nullptr,
                          nbv->grp[eintLocal].kernel_type,
                          nbv->nbat);
        wallcycle_sub_stop(wcycle_, ewcsNBS_GRID_LOCAL);


        nbnxn_atomdata_set(nbv->nbat, nbv->nbs.get(), mdatoms_, fr_->cginfo);

        wallcycle_stop(wcycle_, ewcNS);
    }

    /* do local pair search */
    if (bNS)
    {
        wallcycle_start_nocount(wcycle_, ewcNS);
        wallcycle_sub_start(wcycle_, ewcsNBS_SEARCH_LOCAL);
        nbnxn_make_pairlist(nbv->nbs.get(), nbv->nbat,
                            &top_->excls,
                            nbv->listParams->rlistOuter,
                            nbv->min_ci_balanced,
                            &nbv->grp[eintLocal].nbl_lists,
                            eintLocal,
                            nbv->grp[eintLocal].kernel_type,
                            nrnb_);
        nbv->grp[eintLocal].nbl_lists.outerListCreationStep = step;

        if (nbv->listParams->useDynamicPruning)
        {
            nbnxnPrepareListForDynamicPruning(&nbv->grp[eintLocal].nbl_lists);
        }
        wallcycle_sub_stop(wcycle_, ewcsNBS_SEARCH_LOCAL);

        wallcycle_stop(wcycle_, ewcNS);
    }
    else
    {
        nbnxn_atomdata_copy_x_to_nbat_x(nbv->nbs.get(), eatLocal, FALSE, as_rvec_array(x.data()),
                                        nbv->nbat, wcycle_);
    }

    if (bStateChanged && inputrecNeedMutot(inputrec_))
    {

        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                fr_->mu_tot[i][j] = mu[i*DIM + j];
            }
        }
    }
    if (fr_->efep == efepNO)
    {
        copy_rvec(fr_->mu_tot[0], *mu_tot_);
    }
    else
    {
        for (j = 0; j < DIM; j++)
        {
            *mu_tot_[j] = (1.0 - lambda[efptCOUL])*fr_->mu_tot[0][j] + lambda[efptCOUL]*fr_->mu_tot[1][j];
        }
    }

    /* Reset energies */
    reset_enerdata(enerd_);
    clear_rvecs(SHIFTS, fr_->fshift);

    if (inputrec_->bRot)
    {
        wallcycle_start(wcycle_, ewcROT);
        do_rotation(cr_, enforcedRotation_, *box, as_rvec_array(x.data()), *t, step, bNS);
        wallcycle_stop(wcycle_, ewcROT);
    }

    /* Temporary solution until all routines take PaddedRVecVector */
    auto        force_temp = force_->arrayRefWithPadding();
    rvec *const f          = as_rvec_array(force_temp.unpaddedArrayRef().data());

    /* Start the force cycle counter.
     * Note that a different counter is used for dynamic load balancing.
     */
    wallcycle_start(wcycle_, ewcFORCE);

    gmx::ArrayRef<gmx::RVec> forceRef = force_temp.unpaddedArrayRef();
    if (bDoForces)
    {
        /* If we need to compute the virial, we might need a separate
         * force buffer for algorithms for which the virial is calculated
         * directly, such as PME.
         */
        if ((flags & GMX_FORCE_VIRIAL) && fr_->haveDirectVirialContributions)
        {
            forceRef = *fr_->forceBufferForDirectVirialContributions;

            /* TODO: update comment
             * We only compute forces on local atoms. Note that vsites can
             * spread to non-local atoms, but that part of the buffer is
             * cleared separately in the vsite spreading code.
             */
            clear_rvecs_omp(forceRef.size(), as_rvec_array(forceRef.data()));
        }
        /* Clear the short- and long-range forces */
        clear_rvecs_omp(fr_->natoms_force_constr, f);
    }

    /* forceWithVirial uses the local atom range only */
    gmx::ForceWithVirial forceWithVirial(forceRef, (flags & GMX_FORCE_VIRIAL) != 0);

    if (inputrec_->bPull && pull_have_constraint(inputrec_->pull_work))
    {
        clear_pull_forces(inputrec_->pull_work);
    }

    /* We calculate the non-bonded forces, when done on the CPU, here.
     * We do this before calling do_force_lowlevel, because in that
     * function, the listed forces are calculated before PME, which
     * does communication.  With this order, non-bonded and listed
     * force calculation imbalance can be balanced out by the domain
     * decomposition load balancing.
     */

    if (!bEmulGPU)
    {
        do_nb_verlet(fr_, fr_->ic, enerd_, flags, eintLocal, enbvClearFYes,
                     step, nrnb_, wcycle_);
    }

    if (fr_->efep != efepNO)
    {
        /* Calculate the local and non-local free energy interactions here.
         * Happens here on the CPU both with and without GPU.
         */
        if (fr_->nbv->grp[eintLocal].nbl_lists.nbl_fep[0]->nrj > 0)
        {
            do_nb_verlet_fep(&fr_->nbv->grp[eintLocal].nbl_lists,
                             fr_, as_rvec_array(x.data()), f, mdatoms_,
                             inputrec_->fepvals, lambda,
                             enerd_, flags, nrnb_, wcycle_);
        }
    }

    if (!bEmulGPU)
    {
        int aloc;


        aloc = eintLocal;

        /* Add all the non-bonded force to the normal force array.
         * This can be split into a local and a non-local part when overlapping
         * communication with calculation with domain decomposition.
         */
        wallcycle_stop(wcycle_, ewcFORCE);

        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs.get(), eatAll, nbv->nbat, f, wcycle_);

        wallcycle_start_nocount(wcycle_, ewcFORCE);

        /* if there are multiple fshift output buffers reduce them */
        if ((flags & GMX_FORCE_VIRIAL) &&
            nbv->grp[aloc].nbl_lists.nnbl > 1)
        {
            /* This is not in a subcounter because it takes a
               negligible and constant-sized amount of time */
            nbnxn_atomdata_add_nbat_fshift_to_fshift(nbv->nbat,
                                                     fr_->fshift);
        }
    }

    /* update QMMMrec, if necessary */
    if (fr_->bQMMM)
    {
        update_QMMMrec(cr_, fr_, as_rvec_array(x.data()), mdatoms_, *box);
    }

    /* Compute the bonded and non-bonded energies and optionally forces */
    do_force_lowlevel(fr_, inputrec_, &(top_->idef),
                      cr_, ms_, nrnb_, wcycle_, mdatoms_,
                      as_rvec_array(x.data()), hist, f, &forceWithVirial, enerd_, fcd_,
                      *box, inputrec_->fepvals, lambda, graph_, &(top_->excls), fr_->mu_tot,
                      flags, nullptr);

    wallcycle_stop(wcycle_, ewcFORCE);

    computeSpecialForces(log_, cr_, inputrec_, awh_, enforcedRotation_,
                         step, *t, wcycle_,
                         fr_->forceProviders, *box, x, mdatoms_, lambda,
                         flags, &forceWithVirial, enerd_,
                         ed_, bNS);

    if (bEmulGPU)
    {
        // NOTE: emulation kernel is not included in the balancing region,
        // but emulation mode does not target performance anyway
        wallcycle_start_nocount(wcycle_, ewcFORCE);
        do_nb_verlet(fr_, fr_->ic, enerd_, flags, eintLocal,
                     DOMAINDECOMP(cr_) ? enbvClearFNo : enbvClearFYes,
                     step, nrnb_, wcycle_);
        wallcycle_stop(wcycle_, ewcFORCE);
    }

    /* Do the nonbonded GPU (or emulation) force buffer reduction
     * on the non-alternating path. */
    if (bEmulGPU)
    {
        nbnxn_atomdata_add_nbat_f_to_f(nbv->nbs.get(), eatLocal,
                                       nbv->nbat, f, wcycle_);
    }

    if (bDoForces)
    {
        /* If we have NoVirSum forces, but we do not calculate the virial,
         * we sum fr->f_novirsum=f later.
         */
        if (vsite_ && !(fr_->haveDirectVirialContributions && !(flags & GMX_FORCE_VIRIAL)))
        {
            spread_vsite_f(vsite_, as_rvec_array(x.data()), f, fr_->fshift, FALSE, nullptr, nrnb_,
                           &top_->idef, fr_->ePBC, fr_->bMolPBC, graph_, *box, cr_, wcycle_);
        }

        if (flags & GMX_FORCE_VIRIAL)
        {
            /* Calculation of the virial must be done after vsites! */
            calc_virial(0, mdatoms_->homenr, as_rvec_array(x.data()), f,
                        *vir_force_, graph_, *box, nrnb_, fr_, inputrec_->ePBC);
        }
    }

    if (bDoForces)
    {
        post_process_forces(cr_, step, nrnb_, wcycle_,
                            top_, *box, as_rvec_array(x.data()), f, &forceWithVirial,
                            *vir_force_, mdatoms_, graph_, fr_, vsite_,
                            flags);
    }

    if (flags & GMX_FORCE_ENERGY)
    {
        /* Sum the potential energy terms from group contributions */
        sum_epot(&(enerd_->grpp), enerd_->term);

        if (!EI_TPI(inputrec_->eI))
        {
            checkPotentialEnergyValidity(step, *enerd_, *inputrec_);
        }
    }
}
