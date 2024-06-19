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

#include "force.h"

#include <cassert>
#include <cmath>
#include <cstring>

#include <array>
#include <filesystem>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/ewald.h"
#include "gromacs/ewald/long_range_correction.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

using gmx::ArrayRef;
using gmx::RVec;

static void clearEwaldThreadOutput(ewald_corr_thread_t* ewc_t)
{
    ewc_t->Vcorr_q                                        = 0;
    ewc_t->Vcorr_lj                                       = 0;
    ewc_t->dvdl[FreeEnergyPerturbationCouplingType::Coul] = 0;
    ewc_t->dvdl[FreeEnergyPerturbationCouplingType::Vdw]  = 0;
    clear_mat(ewc_t->vir_q);
    clear_mat(ewc_t->vir_lj);
}

static void reduceEwaldThreadOuput(int nthreads, gmx::ArrayRef<ewald_corr_thread_t> ewc_t)
{
    ewald_corr_thread_t& dest = ewc_t[0];

    for (int t = 1; t < nthreads; t++)
    {
        dest.Vcorr_q += ewc_t[t].Vcorr_q;
        dest.Vcorr_lj += ewc_t[t].Vcorr_lj;
        dest.dvdl[FreeEnergyPerturbationCouplingType::Coul] +=
                ewc_t[t].dvdl[FreeEnergyPerturbationCouplingType::Coul];
        dest.dvdl[FreeEnergyPerturbationCouplingType::Vdw] +=
                ewc_t[t].dvdl[FreeEnergyPerturbationCouplingType::Vdw];
        m_add(dest.vir_q, ewc_t[t].vir_q, dest.vir_q);
        m_add(dest.vir_lj, ewc_t[t].vir_lj, dest.vir_lj);
    }
}

CpuPpLongRangeNonbondeds::CpuPpLongRangeNonbondeds(int                         numberOfTestPaticles,
                                                   real                        ewaldCoeffQ,
                                                   real                        epsilonR,
                                                   gmx::ArrayRef<const double> chargeC6Sum,
                                                   CoulombInteractionType      eeltype,
                                                   VanDerWaalsType             vdwtype,
                                                   const t_inputrec&           inputrec,
                                                   t_nrnb*                     nrnb,
                                                   gmx_wallcycle*              wcycle,
                                                   FILE*                       fplog) :
    numTpiAtoms_(numberOfTestPaticles),
    ewaldCoeffQ_(ewaldCoeffQ),
    epsilonR_(epsilonR),
    chargeC6Sum_(chargeC6Sum),
    coulombInteractionType_(eeltype),
    vanDerWaalsType_(vdwtype),
    ewaldGeometry_(inputrec.ewald_geometry),
    epsilonSurface_(inputrec.epsilon_surface),
    haveEwaldSurfaceTerm_(haveEwaldSurfaceContribution(inputrec)),
    wallEwaldZfac_(inputrec.wall_ewald_zfac),
    havePbcXY2Walls_(inputrecPbcXY2Walls(&inputrec)),
    freeEnergyPerturbationType_(inputrec.efep),
    nrnb_(nrnb),
    wcycle_(wcycle)
{
    GMX_ASSERT(epsilonR == inputrec.epsilon_r,
               "Forcerec and inputrec relative dielectrics don't match");

    // Use the thread count for the bonded module since reducing CPU-side
    // non-bonded contributions does not currently have its own thread count.
    outputPerThread_.resize(gmx_omp_nthreads_get(ModuleMultiThread::Bonded));

    if (inputrec.coulombtype == CoulombInteractionType::Ewald)
    {
        ewaldTable_ = std::make_unique<gmx_ewald_tab_t>(inputrec, fplog);
    }
}

CpuPpLongRangeNonbondeds::~CpuPpLongRangeNonbondeds() = default;

void CpuPpLongRangeNonbondeds::updateAfterPartition(const t_mdatoms& md)
{
    homenr_        = md.homenr;
    havePerturbed_ = md.nChargePerturbed != 0;
    chargeA_       = md.chargeA;
    chargeB_       = md.chargeB;
    sqrt_c6A_      = md.sqrt_c6A;
    sqrt_c6B_      = md.sqrt_c6B;
    sigmaA_        = md.sigmaA;
    sigmaB_        = md.sigmaB;
}

void CpuPpLongRangeNonbondeds::calculate(gmx_pme_t*                     pmedata,
                                         const t_commrec*               commrec,
                                         gmx::ArrayRef<const RVec>      coordinates,
                                         gmx::ForceWithVirial*          forceWithVirial,
                                         gmx_enerdata_t*                enerd,
                                         const matrix                   box,
                                         gmx::ArrayRef<const real>      lambda,
                                         gmx::ArrayRef<const gmx::RVec> mu_tot,
                                         const gmx::StepWorkload&       stepWork,
                                         const DDBalanceRegionHandler&  ddBalanceRegionHandler)
{
    const bool computePmeOnCpu = (usingPme(coulombInteractionType_) || usingLJPme(vanDerWaalsType_))
                                 && thisRankHasDuty(commrec, DUTY_PME)
                                 && (pme_run_mode(pmedata) == PmeRunMode::CPU);

    /* Do long-range electrostatics and/or LJ-PME
     * and compute PME surface terms when necessary.
     */
    if ((computePmeOnCpu || coulombInteractionType_ == CoulombInteractionType::Ewald
         || haveEwaldSurfaceTerm_ || chargeC6Sum_[0] != 0 || chargeC6Sum_[1] != 0)
        && stepWork.computeNonbondedForces)
    {
        real Vlr_q = 0, Vlr_lj = 0;
        /* We reduce all virial, dV/dlambda and energy contributions, except
         * for the reciprocal energies (Vlr_q, Vlr_lj) into the same struct.
         */
        ewald_corr_thread_t& ewaldOutput = outputPerThread_[0];
        clearEwaldThreadOutput(&ewaldOutput);

        if (usingPmeOrEwald(coulombInteractionType_) || usingLJPme(vanDerWaalsType_))
        {
            /* Calculate the Ewald surface force and energy contributions, when necessary */
            if (haveEwaldSurfaceTerm_)
            {
                wallcycle_sub_start(wcycle_, WallCycleSubCounter::EwaldCorrection);

                int nthreads = gmx::ssize(outputPerThread_);
#pragma omp parallel for num_threads(nthreads) schedule(static)
                for (int t = 0; t < nthreads; t++)
                {
                    try
                    {
                        ewald_corr_thread_t& ewc_t = outputPerThread_[t];
                        if (t > 0)
                        {
                            clearEwaldThreadOutput(&ewc_t);
                        }

                        /* Threading is only supported with the Verlet cut-off
                         * scheme and then only single particle forces (no
                         * exclusion forces) are calculated, so we can store
                         * the forces in the normal, single forceWithVirial->force_ array.
                         */
                        ewald_LRcorrection(
                                homenr_,
                                commrec,
                                nthreads,
                                t,
                                epsilonR_,
                                chargeC6Sum_,
                                ewaldGeometry_,
                                epsilonSurface_,
                                havePbcXY2Walls_,
                                wallEwaldZfac_,
                                chargeA_,
                                chargeB_,
                                havePerturbed_,
                                coordinates,
                                box,
                                mu_tot,
                                forceWithVirial->force_,
                                &ewc_t.Vcorr_q,
                                lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                                &ewc_t.dvdl[FreeEnergyPerturbationCouplingType::Coul]);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                if (nthreads > 1)
                {
                    reduceEwaldThreadOuput(nthreads, outputPerThread_);
                }
                wallcycle_sub_stop(wcycle_, WallCycleSubCounter::EwaldCorrection);
            }

            if (usingPmeOrEwald(coulombInteractionType_) && numTpiAtoms_ == 0)
            {
                /* This is not in a subcounter because it takes a
                   negligible and constant-sized amount of time */
                ewaldOutput.Vcorr_q += ewald_charge_correction(
                        commrec,
                        epsilonR_,
                        ewaldCoeffQ_,
                        chargeC6Sum_,
                        lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                        box,
                        &ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Coul],
                        ewaldOutput.vir_q);
            }

            if (computePmeOnCpu)
            {
                /* Do reciprocal PME for Coulomb and/or LJ. */
                assert(numTpiAtoms_ >= 0);
                if (numTpiAtoms_ == 0 || stepWork.stateChanged)
                {
                    /* With domain decomposition we close the CPU side load
                     * balancing region here, because PME does global
                     * communication that acts as a global barrier.
                     */
                    ddBalanceRegionHandler.closeAfterForceComputationCpu();

                    wallcycle_start(wcycle_, WallCycleCounter::PmeMesh);
                    int status = gmx_pme_do(
                            pmedata,
                            gmx::constArrayRefFromArray(coordinates.data(), homenr_ - numTpiAtoms_),
                            forceWithVirial->force_,
                            chargeA_,
                            chargeB_,
                            sqrt_c6A_,
                            sqrt_c6B_,
                            sigmaA_,
                            sigmaB_,
                            box,
                            commrec,
                            haveDDAtomOrdering(*commrec) ? dd_pme_maxshift_x(*commrec->dd) : 0,
                            haveDDAtomOrdering(*commrec) ? dd_pme_maxshift_y(*commrec->dd) : 0,
                            nrnb_,
                            wcycle_,
                            ewaldOutput.vir_q,
                            ewaldOutput.vir_lj,
                            &Vlr_q,
                            &Vlr_lj,
                            lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                            lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                            &ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Coul],
                            &ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Vdw],
                            stepWork);
                    wallcycle_stop(wcycle_, WallCycleCounter::PmeMesh);
                    if (status != 0)
                    {
                        gmx_fatal(FARGS, "Error %d in reciprocal PME routine", status);
                    }

                    /* We should try to do as little computation after
                     * this as possible, because parallel PME synchronizes
                     * the nodes, so we want all load imbalance of the
                     * rest of the force calculation to be before the PME
                     * call.  DD load balancing is done on the whole time
                     * of the force call (without PME).
                     */
                }
                if (numTpiAtoms_ > 0)
                {
                    /* Determine the PME grid energy of the test molecule
                     * with the PME grid potential of the other charges.
                     */
                    Vlr_q = gmx_pme_calc_energy(
                            pmedata,
                            coordinates.subArray(homenr_ - numTpiAtoms_, numTpiAtoms_),
                            chargeA_.subArray(homenr_ - numTpiAtoms_, numTpiAtoms_));
                }
            }
        }

        if (coulombInteractionType_ == CoulombInteractionType::Ewald)
        {
            Vlr_q = do_ewald(havePbcXY2Walls_,
                             wallEwaldZfac_,
                             epsilonR_,
                             freeEnergyPerturbationType_,
                             coordinates,
                             forceWithVirial->force_,
                             chargeA_,
                             chargeB_,
                             box,
                             commrec,
                             homenr_,
                             ewaldOutput.vir_q,
                             ewaldCoeffQ_,
                             lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                             &ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Coul],
                             ewaldTable_.get());
        }

        /* Note that with separate PME nodes we get the real energies later */
        // TODO it would be simpler if we just accumulated a single
        // long-range virial contribution.
        forceWithVirial->addVirialContribution(ewaldOutput.vir_q);
        forceWithVirial->addVirialContribution(ewaldOutput.vir_lj);
        enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Coul] +=
                ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Coul];
        enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Vdw] +=
                ewaldOutput.dvdl[FreeEnergyPerturbationCouplingType::Vdw];
        enerd->term[F_COUL_RECIP] = Vlr_q + ewaldOutput.Vcorr_q;
        enerd->term[F_LJ_RECIP]   = Vlr_lj + ewaldOutput.Vcorr_lj;

        if (debug)
        {
            fprintf(debug,
                    "Vlr_q = %g, Vcorr_q = %g, Vlr_corr_q = %g\n",
                    Vlr_q,
                    ewaldOutput.Vcorr_q,
                    enerd->term[F_COUL_RECIP]);
            pr_rvecs(debug, 0, "vir_el_recip after corr", ewaldOutput.vir_q, DIM);
            fprintf(debug,
                    "Vlr_lj: %g, Vcorr_lj = %g, Vlr_corr_lj = %g\n",
                    Vlr_lj,
                    ewaldOutput.Vcorr_lj,
                    enerd->term[F_LJ_RECIP]);
            pr_rvecs(debug, 0, "vir_lj_recip after corr", ewaldOutput.vir_lj, DIM);
        }
    }

    if (debug)
    {
        print_nrnb(debug, nrnb_);
    }
}
