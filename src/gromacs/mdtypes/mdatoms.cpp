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

#include "mdatoms.h"

#include <cmath>

#include <algorithm>
#include <memory>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/booltype.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vectypes.h"

#define ALMOST_ZERO 1e-30

namespace gmx
{

MDAtoms::MDAtoms(const bool                 rankHasPmeGpuTask,
                 const bool                 useGpuForUpdate,
                 const DeviceStreamManager* deviceStreamManager) :
    // GPU transfers may want to use a suitable pinning mode.
    invmass(makeHostAllocationPolicy(useGpuForUpdate, deviceStreamManager)),
    chargeA(makeHostAllocationPolicy(rankHasPmeGpuTask, deviceStreamManager)),
    chargeB(makeHostAllocationPolicy(rankHasPmeGpuTask, deviceStreamManager)),
    cTC(makeHostAllocationPolicy(useGpuForUpdate, deviceStreamManager))
{
    if (rankHasPmeGpuTask || useGpuForUpdate)
    {
        GMX_RELEASE_ASSERT(deviceStreamManager != nullptr,
                           "Must have device stream manager when there is a GPU task");
    }
}


std::unique_ptr<MDAtoms> makeMDAtoms(FILE*                      fp,
                                     const gmx_mtop_t&          mtop,
                                     const t_inputrec&          ir,
                                     const bool                 rankHasPmeGpuTask,
                                     const bool                 useGpuForUpdate,
                                     const DeviceStreamManager* deviceStreamManager)
{
    auto mdAtoms = std::make_unique<MDAtoms>(rankHasPmeGpuTask, useGpuForUpdate, deviceStreamManager);

    mdAtoms->bVCMgrps = FALSE;
    for (int i = 0; i < mtop.natoms; i++)
    {
        if (getGroupType(mtop.groups, SimulationAtomGroupType::MassCenterVelocityRemoval, i) > 0)
        {
            mdAtoms->bVCMgrps = TRUE;
        }
    }

    /* Determine the total system mass and perturbed atom counts */
    double totalMassA = 0.0;
    double totalMassB = 0.0;

    mdAtoms->haveVsites             = FALSE;
    gmx_mtop_atomloop_block_t aloop = gmx_mtop_atomloop_block_init(mtop);
    const t_atom*             atom;
    int                       nmol;

    mdAtoms->nPerturbed       = 0;
    mdAtoms->nMassPerturbed   = 0;
    mdAtoms->nChargePerturbed = 0;
    mdAtoms->nTypePerturbed   = 0;
    while (gmx_mtop_atomloop_block_next(aloop, &atom, &nmol))
    {
        totalMassA += nmol * atom->m;
        totalMassB += nmol * atom->mB;

        if (atom->ptype == ParticleType::VSite)
        {
            mdAtoms->haveVsites = TRUE;
        }

        if (ir.efep != FreeEnergyPerturbationType::No && PERTURBED(*atom))
        {
            mdAtoms->nPerturbed++;
            if (atom->mB != atom->m)
            {
                mdAtoms->nMassPerturbed += nmol;
            }
            if (atom->qB != atom->q)
            {
                mdAtoms->nChargePerturbed += nmol;
            }
            if (atom->typeB != atom->type)
            {
                mdAtoms->nTypePerturbed += nmol;
            }
        }
    }

    mdAtoms->tmassA = totalMassA;
    mdAtoms->tmassB = totalMassB;

    if (ir.efep != FreeEnergyPerturbationType::No && fp)
    {
        fprintf(fp,
                "There are %d atoms and %d charges for free energy perturbation\n",
                mdAtoms->nPerturbed,
                mdAtoms->nChargePerturbed);
    }

    mdAtoms->havePartiallyFrozenAtoms = FALSE;
    for (int g = 0; g < ir.opts.ngfrz; g++)
    {
        for (int d = YY; d < DIM; d++)
        {
            if (ir.opts.nFreeze[g][d] != ir.opts.nFreeze[g][XX])
            {
                mdAtoms->havePartiallyFrozenAtoms = TRUE;
            }
        }
    }

    mdAtoms->bOrires = (gmx_mtop_ftype_count(mtop, InteractionFunction::OrientationRestraints) != 0);

    return mdAtoms;
}

} // namespace gmx

/* Only the GPU update and constraints copy invmass and cTC to the GPU, asynchronously.
 * atoms2md() can reallocate them, because the search step first waits for the coordinates
 * on the host, so those copies are done. */
void atoms2md(const gmx_mtop_t&        mtop,
              const t_inputrec&        inputrec,
              int                      nindex,
              gmx::ArrayRef<const int> index,
              int                      numHomeAtoms,
              gmx::MDAtoms*            mdAtoms)
{
    gmx_bool         bLJPME;
    const t_grpopts* opts;

    bLJPME = usingLJPme(inputrec.vdwtype);

    opts = &inputrec.opts;

    const SimulationGroups& groups = mtop.groups;

    // When using DD, nindex (>= 0) indicates the size of the map
    // from local to global atom indices.  MDAtoms needs to allocate
    // space for home atoms and ghost atoms for force, constraint, and
    // virtual-site operations.
    const int numTotalAtoms = (nindex >= 0) ? nindex : mtop.natoms;

    if (mdAtoms->nMassPerturbed)
    {
        mdAtoms->massA.resize(numTotalAtoms);
        mdAtoms->massB.resize(numTotalAtoms);
    }
    mdAtoms->massT.resize(numTotalAtoms);
    mdAtoms->invmass.resizeWithPadding(numTotalAtoms);
    mdAtoms->invMassPerDim.resize(numTotalAtoms);
    mdAtoms->chargeA.resizeWithPadding(numTotalAtoms);
    if (mdAtoms->nPerturbed > 0)
    {
        mdAtoms->chargeB.resizeWithPadding(numTotalAtoms);
    }
    mdAtoms->typeA.resize(numTotalAtoms);
    if (mdAtoms->nPerturbed)
    {
        mdAtoms->typeB.resize(numTotalAtoms);
    }
    if (bLJPME)
    {
        mdAtoms->sqrt_c6A.resize(numTotalAtoms);
        mdAtoms->sigmaA.resize(numTotalAtoms);
        mdAtoms->sigma3A.resize(numTotalAtoms);
        if (mdAtoms->nPerturbed)
        {
            mdAtoms->sqrt_c6B.resize(numTotalAtoms);
            mdAtoms->sigmaB.resize(numTotalAtoms);
            mdAtoms->sigma3B.resize(numTotalAtoms);
        }
    }
    mdAtoms->ptype.resize(numTotalAtoms);
    if (opts->ngtc > 1)
    {
        mdAtoms->cTC.resize(numTotalAtoms);
        /* We always copy cTC with domain decomposition */
    }
    mdAtoms->cENER.resize(numTotalAtoms);
    if (inputrec.useConstantAcceleration)
    {
        mdAtoms->cACC.resize(numTotalAtoms);
    }
    if (inputrecFrozenAtoms(&inputrec))
    {
        mdAtoms->cFREEZE.resize(numTotalAtoms);
    }
    if (mdAtoms->bVCMgrps)
    {
        mdAtoms->cVCM.resize(numTotalAtoms);
    }
    if (mdAtoms->bOrires)
    {
        mdAtoms->cORF.resize(numTotalAtoms);
    }
    if (mdAtoms->nPerturbed)
    {
        mdAtoms->bPerturbed.resize(numTotalAtoms);
    }

    // Note that these user groups are empty
    // when there is only one group present.
    // Therefore, when adding code, the user should use something like:
    // gprnrU1 = (mdAtoms->cU1.empty() ? 0 : mdAtoms->cU1[localatindex])
    if (!mtop.groups.groupNumbers[SimulationAtomGroupType::User1].empty())
    {
        mdAtoms->cU1.resize(numTotalAtoms);
    }
    if (!mtop.groups.groupNumbers[SimulationAtomGroupType::User2].empty())
    {
        mdAtoms->cU2.resize(numTotalAtoms);
    }

    MTopLookUp mTopLookUp(mtop);

    const unsigned short numTypes   = mtop.ffparams.atnr;
    const t_atom         fillerAtom = {
        0, 0, 0, 0, numTypes, numTypes, ParticleType::Count, -1, 0,
    };

    // In grompp, OpenMP is not initialized and nthreads_get returns 0. We want 1 thread in this case.
    const int gmx_unused nthreads = std::max(gmx_omp_nthreads_get(ModuleMultiThread::Default), 1);
#pragma omp parallel for num_threads(nthreads) schedule(static) firstprivate(mTopLookUp)
    for (int i = 0; i < numTotalAtoms; i++)
    {
        try
        {
            int  g, ag;
            real mA, mB, fac;

            if (index.empty())
            {
                ag = i;
            }
            else
            {
                ag = index[i];
            }
            const bool isValidAtom = isValidGlobalAtom(ag);

            const t_atom& atom = (isValidAtom ? mTopLookUp.getAtomParameters(ag) : fillerAtom);

            if (!mdAtoms->cFREEZE.empty())
            {
                mdAtoms->cFREEZE[i] =
                        (isValidAtom ? getGroupType(groups, SimulationAtomGroupType::Freeze, ag) : 0);
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
                          / opts->tau_t[!mdAtoms->cTC.empty() ? groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag]
                                                              : 0];
                    mA = 0.5 * atom.m * fac;
                    mB = 0.5 * atom.mB * fac;
                }
            }
            else
            {
                mA = atom.m;
                mB = atom.mB;
            }
            if (mdAtoms->nMassPerturbed)
            {
                mdAtoms->massA[i] = mA;
                mdAtoms->massB[i] = mB;
            }
            mdAtoms->massT[i] = mA;

            if (mA == 0.0)
            {
                mdAtoms->invmass[i]           = 0;
                mdAtoms->invMassPerDim[i][XX] = 0;
                mdAtoms->invMassPerDim[i][YY] = 0;
                mdAtoms->invMassPerDim[i][ZZ] = 0;
            }
            else if (!mdAtoms->cFREEZE.empty())
            {
                g = mdAtoms->cFREEZE[i];
                GMX_ASSERT(opts->nFreeze != nullptr, "Must have freeze groups to initialize masses");
                if (opts->nFreeze[g][XX] && opts->nFreeze[g][YY] && opts->nFreeze[g][ZZ])
                {
                    /* Set the mass of completely frozen particles to ALMOST_ZERO
                     * iso 0 to avoid div by zero in lincs or shake.
                     */
                    mdAtoms->invmass[i] = ALMOST_ZERO;
                }
                else
                {
                    /* Note: Partially frozen particles use the normal invmass.
                     * If such particles are constrained, the frozen dimensions
                     * should not be updated with the constrained coordinates.
                     */
                    mdAtoms->invmass[i] = 1.0 / mA;
                }
                for (int d = 0; d < DIM; d++)
                {
                    mdAtoms->invMassPerDim[i][d] = (opts->nFreeze[g][d] ? 0 : 1.0 / mA);
                }
            }
            else
            {
                mdAtoms->invmass[i] = 1.0 / mA;
                for (int d = 0; d < DIM; d++)
                {
                    mdAtoms->invMassPerDim[i][d] = 1.0 / mA;
                }
            }

            mdAtoms->chargeA[i] = atom.q;
            mdAtoms->typeA[i]   = atom.type;
            if (bLJPME)
            {
                real c6              = (isValidAtom
                                                ? mtop.ffparams.iparams[atom.type * (mtop.ffparams.atnr + 1)].lj.c6
                                                : 0.0_real);
                real c12             = (isValidAtom
                                                ? mtop.ffparams.iparams[atom.type * (mtop.ffparams.atnr + 1)].lj.c12
                                                : 0.0_real);
                mdAtoms->sqrt_c6A[i] = std::sqrt(c6);
                if (c6 == 0.0 || c12 == 0)
                {
                    mdAtoms->sigmaA[i] = 1.0;
                }
                else
                {
                    mdAtoms->sigmaA[i] = gmx::sixthroot(c12 / c6);
                }
                mdAtoms->sigma3A[i] = 1 / (mdAtoms->sigmaA[i] * mdAtoms->sigmaA[i] * mdAtoms->sigmaA[i]);
            }
            if (mdAtoms->nPerturbed)
            {
                mdAtoms->bPerturbed[i] = PERTURBED(atom);
                mdAtoms->chargeB[i]    = atom.qB;
                mdAtoms->typeB[i]      = atom.typeB;
                if (bLJPME)
                {
                    real c6              = (isValidAtom ? mtop.ffparams
                                                     .iparams[atom.typeB * (mtop.ffparams.atnr + 1)]
                                                     .lj.c6
                                                        : 0.0_real);
                    real c12             = (isValidAtom ? mtop.ffparams
                                                      .iparams[atom.typeB * (mtop.ffparams.atnr + 1)]
                                                      .lj.c12
                                                        : 0.0_real);
                    mdAtoms->sqrt_c6B[i] = std::sqrt(c6);
                    if (c6 == 0.0 || c12 == 0)
                    {
                        mdAtoms->sigmaB[i] = 1.0;
                    }
                    else
                    {
                        mdAtoms->sigmaB[i] = gmx::sixthroot(c12 / c6);
                    }
                    mdAtoms->sigma3B[i] =
                            1 / (mdAtoms->sigmaB[i] * mdAtoms->sigmaB[i] * mdAtoms->sigmaB[i]);
                }
            }
            mdAtoms->ptype[i] = atom.ptype;

            if (isValidAtom)
            {
                if (!mdAtoms->cTC.empty())
                {
                    mdAtoms->cTC[i] =
                            groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag];
                }
                mdAtoms->cENER[i] = getGroupType(groups, SimulationAtomGroupType::EnergyOutput, ag);
                if (!mdAtoms->cACC.empty())
                {
                    mdAtoms->cACC[i] = groups.groupNumbers[SimulationAtomGroupType::Acceleration][ag];
                }
                if (!mdAtoms->cVCM.empty())
                {
                    mdAtoms->cVCM[i] =
                            groups.groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval][ag];
                }
                if (!mdAtoms->cORF.empty())
                {
                    mdAtoms->cORF[i] =
                            getGroupType(groups, SimulationAtomGroupType::OrientationRestraintsFit, ag);
                }

                if (!mdAtoms->cU1.empty())
                {
                    mdAtoms->cU1[i] = groups.groupNumbers[SimulationAtomGroupType::User1][ag];
                }
                if (!mdAtoms->cU2.empty())
                {
                    mdAtoms->cU2[i] = groups.groupNumbers[SimulationAtomGroupType::User2][ag];
                }
            }
            else
            {
                // As fillers have no mass and interactions, we can add them to group 0 without side-effects
                if (!mdAtoms->cTC.empty())
                {
                    mdAtoms->cTC[i] = 0;
                }
                mdAtoms->cENER[i] = 0;
                if (!mdAtoms->cACC.empty())
                {
                    mdAtoms->cACC[i] = 0;
                }
                if (!mdAtoms->cVCM.empty())
                {
                    mdAtoms->cVCM[i] = 0;
                }
                GMX_ASSERT(mdAtoms->cORF.empty(),
                           "Combination of orientation restraints and fillers is not supported");

                if (!mdAtoms->cU1.empty())
                {
                    mdAtoms->cU1[i] = -1;
                }
                if (!mdAtoms->cU2.empty())
                {
                    mdAtoms->cU2[i] = -1;
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    if (numTotalAtoms > 0)
    {
        /* Pad invmass with 0 so a SIMD MD update does not change v and x */
        for (int i = numTotalAtoms; i < mdAtoms->invmass.paddedSize(); i++)
        {
            mdAtoms->invmass[i] = 0;
        }
    }

    mdAtoms->numHomeAtoms = numHomeAtoms;
    /* We set mass, invmass, invMassPerDim and tmass for lambda=0.
     * For free-energy runs, these should be updated using update_mdatoms().
     */
    mdAtoms->tmass      = mdAtoms->tmassA;
    mdAtoms->massLambda = 0;
}

void update_mdatoms(gmx::MDAtoms* mdAtoms, real massLambda)
{
    if (mdAtoms->nMassPerturbed && massLambda != mdAtoms->massLambda)
    {
        real L1 = 1 - massLambda;

        const int numTotalAtoms = mdAtoms->massT.size();

        /* Update masses of perturbed atoms for the change in mass lambda */
        int gmx_unused nthreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0; i < numTotalAtoms; i++)
        {
            if (mdAtoms->bPerturbed[i])
            {
                mdAtoms->massT[i] = L1 * mdAtoms->massA[i] + massLambda * mdAtoms->massB[i];
                /* Atoms with invmass 0 or ALMOST_ZERO are massless or frozen
                 * and their invmass does not depend on lambda.
                 */
                if (mdAtoms->invmass[i] > 1.1 * ALMOST_ZERO)
                {
                    mdAtoms->invmass[i] = 1.0 / mdAtoms->massT[i];
                    for (int d = 0; d < DIM; d++)
                    {
                        if (mdAtoms->invMassPerDim[i][d] > 1.1 * ALMOST_ZERO)
                        {
                            mdAtoms->invMassPerDim[i][d] = mdAtoms->invmass[i];
                        }
                    }
                }
            }
        }

        /* Update the system mass for the change in mass lambda */
        mdAtoms->tmass = L1 * mdAtoms->tmassA + massLambda * mdAtoms->tmassB;
    }

    mdAtoms->massLambda = massLambda;
}
