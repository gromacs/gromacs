/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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

#include "mdatoms.h"

#include <cmath>

#include <memory>

#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#define ALMOST_ZERO 1e-30

namespace gmx
{

MDAtoms::MDAtoms() : mdatoms_(nullptr) {}

MDAtoms::~MDAtoms()
{
    if (mdatoms_ == nullptr)
    {
        return;
    }
    sfree(mdatoms_->massA);
    sfree(mdatoms_->massB);
    sfree(mdatoms_->massT);
    gmx::AlignedAllocationPolicy::free(mdatoms_->invmass);
    sfree(mdatoms_->invMassPerDim);
    sfree(mdatoms_->typeA);
    sfree(mdatoms_->typeB);
    /* mdatoms->chargeA and mdatoms->chargeB point at chargeA_.data()
     * and chargeB_.data() respectively. They get freed automatically. */
    sfree(mdatoms_->sqrt_c6A);
    sfree(mdatoms_->sigmaA);
    sfree(mdatoms_->sigma3A);
    sfree(mdatoms_->sqrt_c6B);
    sfree(mdatoms_->sigmaB);
    sfree(mdatoms_->sigma3B);
    sfree(mdatoms_->ptype);
    sfree(mdatoms_->cTC);
    sfree(mdatoms_->cENER);
    sfree(mdatoms_->cACC);
    sfree(mdatoms_->cFREEZE);
    sfree(mdatoms_->cVCM);
    sfree(mdatoms_->cORF);
    sfree(mdatoms_->bPerturbed);
    sfree(mdatoms_->cU1);
    sfree(mdatoms_->cU2);
}

void MDAtoms::resizeChargeA(const int newSize)
{
    chargeA_.resizeWithPadding(newSize);
    mdatoms_->chargeA = chargeA_.data();
}

void MDAtoms::resizeChargeB(const int newSize)
{
    chargeB_.resizeWithPadding(newSize);
    mdatoms_->chargeB = chargeB_.data();
}

void MDAtoms::reserveChargeA(const int newCapacity)
{
    chargeA_.reserveWithPadding(newCapacity);
    mdatoms_->chargeA = chargeA_.data();
}

void MDAtoms::reserveChargeB(const int newCapacity)
{
    chargeB_.reserveWithPadding(newCapacity);
    mdatoms_->chargeB = chargeB_.data();
}

std::unique_ptr<MDAtoms> makeMDAtoms(FILE* fp, const gmx_mtop_t& mtop, const t_inputrec& ir, const bool rankHasPmeGpuTask)
{
    auto mdAtoms = std::make_unique<MDAtoms>();
    // GPU transfers may want to use a suitable pinning mode.
    if (rankHasPmeGpuTask)
    {
        changePinningPolicy(&mdAtoms->chargeA_, pme_get_pinning_policy());
        changePinningPolicy(&mdAtoms->chargeB_, pme_get_pinning_policy());
    }
    t_mdatoms* md;
    snew(md, 1);
    mdAtoms->mdatoms_.reset(md);

    md->nenergrp = mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size();
    md->bVCMgrps = FALSE;
    for (int i = 0; i < mtop.natoms; i++)
    {
        if (getGroupType(mtop.groups, SimulationAtomGroupType::MassCenterVelocityRemoval, i) > 0)
        {
            md->bVCMgrps = TRUE;
        }
    }

    /* Determine the total system mass and perturbed atom counts */
    double totalMassA = 0.0;
    double totalMassB = 0.0;

    md->haveVsites                  = FALSE;
    gmx_mtop_atomloop_block_t aloop = gmx_mtop_atomloop_block_init(mtop);
    const t_atom*             atom;
    int                       nmol;
    while (gmx_mtop_atomloop_block_next(aloop, &atom, &nmol))
    {
        totalMassA += nmol * atom->m;
        totalMassB += nmol * atom->mB;

        if (atom->ptype == ParticleType::VSite)
        {
            md->haveVsites = TRUE;
        }

        if (ir.efep != FreeEnergyPerturbationType::No && PERTURBED(*atom))
        {
            md->nPerturbed++;
            if (atom->mB != atom->m)
            {
                md->nMassPerturbed += nmol;
            }
            if (atom->qB != atom->q)
            {
                md->nChargePerturbed += nmol;
            }
            if (atom->typeB != atom->type)
            {
                md->nTypePerturbed += nmol;
            }
        }
    }

    md->tmassA = totalMassA;
    md->tmassB = totalMassB;

    if (ir.efep != FreeEnergyPerturbationType::No && fp)
    {
        fprintf(fp,
                "There are %d atoms and %d charges for free energy perturbation\n",
                md->nPerturbed,
                md->nChargePerturbed);
    }

    md->havePartiallyFrozenAtoms = FALSE;
    for (int g = 0; g < ir.opts.ngfrz; g++)
    {
        for (int d = YY; d < DIM; d++)
        {
            if (ir.opts.nFreeze[g][d] != ir.opts.nFreeze[g][XX])
            {
                md->havePartiallyFrozenAtoms = TRUE;
            }
        }
    }

    md->bOrires = (gmx_mtop_ftype_count(mtop, F_ORIRES) != 0);

    return mdAtoms;
}

} // namespace gmx

void atoms2md(const gmx_mtop_t&  mtop,
              const t_inputrec&  inputrec,
              int                nindex,
              gmx::ArrayRef<int> index,
              int                homenr,
              gmx::MDAtoms*      mdAtoms)
{
    gmx_bool         bLJPME;
    const t_grpopts* opts;
    int nthreads     gmx_unused;

    bLJPME = EVDW_PME(inputrec.vdwtype);

    opts = &inputrec.opts;

    const SimulationGroups& groups = mtop.groups;

    auto* md = mdAtoms->mdatoms();
    /* nindex>=0 indicates DD where we use an index */
    if (nindex >= 0)
    {
        md->nr = nindex;
    }
    else
    {
        md->nr = mtop.natoms;
    }

    if (md->nr > md->nalloc)
    {
        md->nalloc = over_alloc_dd(md->nr);

        if (md->nMassPerturbed)
        {
            srenew(md->massA, md->nalloc);
            srenew(md->massB, md->nalloc);
        }
        srenew(md->massT, md->nalloc);
        /* The SIMD version of the integrator needs this aligned and padded.
         * The padding needs to be with zeros, which we set later below.
         */
        gmx::AlignedAllocationPolicy::free(md->invmass);
        md->invmass = new (gmx::AlignedAllocationPolicy::malloc(
                (md->nalloc + GMX_REAL_MAX_SIMD_WIDTH) * sizeof(*md->invmass))) real;
        srenew(md->invMassPerDim, md->nalloc);
        // TODO eventually we will have vectors and just resize
        // everything, but for now the semantics of md->nalloc being
        // the capacity are preserved by keeping vectors within
        // mdAtoms having the same properties as the other arrays.
        mdAtoms->reserveChargeA(md->nalloc);
        mdAtoms->resizeChargeA(md->nr);
        if (md->nPerturbed > 0)
        {
            mdAtoms->reserveChargeB(md->nalloc);
            mdAtoms->resizeChargeB(md->nr);
        }
        srenew(md->typeA, md->nalloc);
        if (md->nPerturbed)
        {
            srenew(md->typeB, md->nalloc);
        }
        if (bLJPME)
        {
            srenew(md->sqrt_c6A, md->nalloc);
            srenew(md->sigmaA, md->nalloc);
            srenew(md->sigma3A, md->nalloc);
            if (md->nPerturbed)
            {
                srenew(md->sqrt_c6B, md->nalloc);
                srenew(md->sigmaB, md->nalloc);
                srenew(md->sigma3B, md->nalloc);
            }
        }
        srenew(md->ptype, md->nalloc);
        if (opts->ngtc > 1)
        {
            srenew(md->cTC, md->nalloc);
            /* We always copy cTC with domain decomposition */
        }
        srenew(md->cENER, md->nalloc);
        if (inputrec.useConstantAcceleration)
        {
            srenew(md->cACC, md->nalloc);
        }
        if (inputrecFrozenAtoms(&inputrec))
        {
            srenew(md->cFREEZE, md->nalloc);
        }
        if (md->bVCMgrps)
        {
            srenew(md->cVCM, md->nalloc);
        }
        if (md->bOrires)
        {
            srenew(md->cORF, md->nalloc);
        }
        if (md->nPerturbed)
        {
            srenew(md->bPerturbed, md->nalloc);
        }

        /* Note that these user t_mdatoms array pointers are NULL
         * when there is only one group present.
         * Therefore, when adding code, the user should use something like:
         * gprnrU1 = (md->cU1==NULL ? 0 : md->cU1[localatindex])
         */
        if (!mtop.groups.groupNumbers[SimulationAtomGroupType::User1].empty())
        {
            srenew(md->cU1, md->nalloc);
        }
        if (!mtop.groups.groupNumbers[SimulationAtomGroupType::User2].empty())
        {
            srenew(md->cU2, md->nalloc);
        }
    }

    int molb = 0;

    nthreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(nthreads) schedule(static) firstprivate(molb)
    for (int i = 0; i < md->nr; i++)
    {
        try
        {
            int  g, ag;
            real mA, mB, fac;
            real c6, c12;

            if (index.empty())
            {
                ag = i;
            }
            else
            {
                ag = index[i];
            }
            const t_atom& atom = mtopGetAtomParameters(mtop, ag, &molb);

            if (md->cFREEZE)
            {
                md->cFREEZE[i] = getGroupType(groups, SimulationAtomGroupType::Freeze, ag);
            }
            if (EI_ENERGY_MINIMIZATION(inputrec.eI))
            {
                /* Displacement is proportional to F, masses used for constraints */
                mA = 1.0;
                mB = 1.0;
            }
            else if (inputrec.eI == IntegrationAlgorithm::BD)
            {
                /* With BD the physical masses are irrelevant.
                 * To keep the code simple we use most of the normal MD code path
                 * for BD. Thus for constraining the masses should be proportional
                 * to the friction coefficient. We set the absolute value such that
                 * m/2<(dx/dt)^2> = m/2*2kT/fric*dt = kT/2 => m=fric*dt/2
                 * Then if we set the (meaningless) velocity to v=dx/dt, we get the
                 * correct kinetic energy and temperature using the usual code path.
                 * Thus with BD v*dt will give the displacement and the reported
                 * temperature can signal bad integration (too large time step).
                 */
                if (inputrec.bd_fric > 0)
                {
                    mA = 0.5 * inputrec.bd_fric * inputrec.delta_t;
                    mB = 0.5 * inputrec.bd_fric * inputrec.delta_t;
                }
                else
                {
                    /* The friction coefficient is mass/tau_t */
                    fac = inputrec.delta_t
                          / opts->tau_t[md->cTC ? groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag] : 0];
                    mA = 0.5 * atom.m * fac;
                    mB = 0.5 * atom.mB * fac;
                }
            }
            else
            {
                mA = atom.m;
                mB = atom.mB;
            }
            if (md->nMassPerturbed)
            {
                md->massA[i] = mA;
                md->massB[i] = mB;
            }
            md->massT[i] = mA;

            if (mA == 0.0)
            {
                md->invmass[i]           = 0;
                md->invMassPerDim[i][XX] = 0;
                md->invMassPerDim[i][YY] = 0;
                md->invMassPerDim[i][ZZ] = 0;
            }
            else if (md->cFREEZE)
            {
                g = md->cFREEZE[i];
                GMX_ASSERT(opts->nFreeze != nullptr, "Must have freeze groups to initialize masses");
                if (opts->nFreeze[g][XX] && opts->nFreeze[g][YY] && opts->nFreeze[g][ZZ])
                {
                    /* Set the mass of completely frozen particles to ALMOST_ZERO
                     * iso 0 to avoid div by zero in lincs or shake.
                     */
                    md->invmass[i] = ALMOST_ZERO;
                }
                else
                {
                    /* Note: Partially frozen particles use the normal invmass.
                     * If such particles are constrained, the frozen dimensions
                     * should not be updated with the constrained coordinates.
                     */
                    md->invmass[i] = 1.0 / mA;
                }
                for (int d = 0; d < DIM; d++)
                {
                    md->invMassPerDim[i][d] = (opts->nFreeze[g][d] ? 0 : 1.0 / mA);
                }
            }
            else
            {
                md->invmass[i] = 1.0 / mA;
                for (int d = 0; d < DIM; d++)
                {
                    md->invMassPerDim[i][d] = 1.0 / mA;
                }
            }

            md->chargeA[i] = atom.q;
            md->typeA[i]   = atom.type;
            if (bLJPME)
            {
                c6  = mtop.ffparams.iparams[atom.type * (mtop.ffparams.atnr + 1)].lj.c6;
                c12 = mtop.ffparams.iparams[atom.type * (mtop.ffparams.atnr + 1)].lj.c12;
                md->sqrt_c6A[i] = std::sqrt(c6);
                if (c6 == 0.0 || c12 == 0)
                {
                    md->sigmaA[i] = 1.0;
                }
                else
                {
                    md->sigmaA[i] = gmx::sixthroot(c12 / c6);
                }
                md->sigma3A[i] = 1 / (md->sigmaA[i] * md->sigmaA[i] * md->sigmaA[i]);
            }
            if (md->nPerturbed)
            {
                md->bPerturbed[i] = PERTURBED(atom);
                md->chargeB[i]    = atom.qB;
                md->typeB[i]      = atom.typeB;
                if (bLJPME)
                {
                    c6  = mtop.ffparams.iparams[atom.typeB * (mtop.ffparams.atnr + 1)].lj.c6;
                    c12 = mtop.ffparams.iparams[atom.typeB * (mtop.ffparams.atnr + 1)].lj.c12;
                    md->sqrt_c6B[i] = std::sqrt(c6);
                    if (c6 == 0.0 || c12 == 0)
                    {
                        md->sigmaB[i] = 1.0;
                    }
                    else
                    {
                        md->sigmaB[i] = gmx::sixthroot(c12 / c6);
                    }
                    md->sigma3B[i] = 1 / (md->sigmaB[i] * md->sigmaB[i] * md->sigmaB[i]);
                }
            }
            md->ptype[i] = atom.ptype;
            if (md->cTC)
            {
                md->cTC[i] = groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag];
            }
            md->cENER[i] = getGroupType(groups, SimulationAtomGroupType::EnergyOutput, ag);
            if (md->cACC)
            {
                md->cACC[i] = groups.groupNumbers[SimulationAtomGroupType::Acceleration][ag];
            }
            if (md->cVCM)
            {
                md->cVCM[i] = groups.groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval][ag];
            }
            if (md->cORF)
            {
                md->cORF[i] = getGroupType(groups, SimulationAtomGroupType::OrientationRestraintsFit, ag);
            }

            if (md->cU1)
            {
                md->cU1[i] = groups.groupNumbers[SimulationAtomGroupType::User1][ag];
            }
            if (md->cU2)
            {
                md->cU2[i] = groups.groupNumbers[SimulationAtomGroupType::User2][ag];
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    if (md->nr > 0)
    {
        /* Pad invmass with 0 so a SIMD MD update does not change v and x */
        for (int i = md->nr; i < md->nr + GMX_REAL_MAX_SIMD_WIDTH; i++)
        {
            md->invmass[i] = 0;
        }
    }

    md->homenr = homenr;
    /* We set mass, invmass, invMassPerDim and tmass for lambda=0.
     * For free-energy runs, these should be updated using update_mdatoms().
     */
    md->tmass  = md->tmassA;
    md->lambda = 0;
}

void update_mdatoms(t_mdatoms* md, real lambda)
{
    if (md->nMassPerturbed && lambda != md->lambda)
    {
        real L1 = 1 - lambda;

        /* Update masses of perturbed atoms for the change in lambda */
        int gmx_unused nthreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0; i < md->nr; i++)
        {
            if (md->bPerturbed[i])
            {
                md->massT[i] = L1 * md->massA[i] + lambda * md->massB[i];
                /* Atoms with invmass 0 or ALMOST_ZERO are massless or frozen
                 * and their invmass does not depend on lambda.
                 */
                if (md->invmass[i] > 1.1 * ALMOST_ZERO)
                {
                    md->invmass[i] = 1.0 / md->massT[i];
                    for (int d = 0; d < DIM; d++)
                    {
                        if (md->invMassPerDim[i][d] > 1.1 * ALMOST_ZERO)
                        {
                            md->invMassPerDim[i][d] = md->invmass[i];
                        }
                    }
                }
            }
        }

        /* Update the system mass for the change in lambda */
        md->tmass = L1 * md->tmassA + lambda * md->tmassB;
    }

    md->lambda = lambda;
}
