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

MDAtoms::MDAtoms(const gmx_mtop_t&          mtop,
                 const t_inputrec&          inputrec,
                 const bool                 rankHasPmeGpuTask,
                 const bool                 useGpuForUpdate,
                 const DeviceStreamManager* deviceStreamManager) :
    mtop_(mtop),
    inputrec_(inputrec),
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

    for (int i = 0; i < mtop.natoms; i++)
    {
        if (getGroupType(mtop.groups, SimulationAtomGroupType::MassCenterVelocityRemoval, i) > 0)
        {
            bVCMgrps = true;
            break;
        }
    }

    // Determine the total system mass and perturbed atom counts
    double totalMassA = 0.0;
    double totalMassB = 0.0;

    gmx_mtop_atomloop_block_t aloop = gmx_mtop_atomloop_block_init(mtop);
    const t_atom*             atom;
    int                       nmol;

    while (gmx_mtop_atomloop_block_next(aloop, &atom, &nmol))
    {
        totalMassA += nmol * atom->m;
        totalMassB += nmol * atom->mB;

        if (inputrec_.efep != FreeEnergyPerturbationType::No && PERTURBED(*atom))
        {
            nPerturbed++;
            if (atom->mB != atom->m)
            {
                nMassPerturbed += nmol;
            }
            if (atom->qB != atom->q)
            {
                nChargePerturbed += nmol;
            }
            if (atom->typeB != atom->type)
            {
                nTypePerturbed += nmol;
            }
        }
    }

    tmassA = totalMassA;
    tmassB = totalMassB;
}

/* Only the GPU update and constraints copy invmass and cTC to the GPU, asynchronously.
 * atoms2md() can reallocate them, because the search step first waits for the coordinates
 * on the host, so those copies are done. */
void MDAtoms::updateAfterPartition(const int           numTotalAtoms,
                                   ArrayRef<const int> index,
                                   const int           numHomeAtomsArg,
                                   const real          massLambda)
{
    gmx_bool         bLJPME;
    const t_grpopts* opts;

    bLJPME = usingLJPme(inputrec_.vdwtype);

    opts = &inputrec_.opts;

    const SimulationGroups& groups = mtop_.groups;

    // MDAtoms needs to allocate space for home atoms and ghost atoms
    // for force, constraint, and virtual-site operations.

    if (nMassPerturbed)
    {
        massA.resize(numTotalAtoms);
        massB.resize(numTotalAtoms);
    }
    massT.resize(numTotalAtoms);
    invmass.resizeWithPadding(numTotalAtoms);
    invMassPerDim.resize(numTotalAtoms);
    chargeA.resizeWithPadding(numTotalAtoms);
    if (nPerturbed > 0)
    {
        chargeB.resizeWithPadding(numTotalAtoms);
    }
    typeA.resize(numTotalAtoms);
    if (nPerturbed)
    {
        typeB.resize(numTotalAtoms);
    }
    if (bLJPME)
    {
        sqrt_c6A.resize(numTotalAtoms);
        sigmaA.resize(numTotalAtoms);
        sigma3A.resize(numTotalAtoms);
        if (nPerturbed)
        {
            sqrt_c6B.resize(numTotalAtoms);
            sigmaB.resize(numTotalAtoms);
            sigma3B.resize(numTotalAtoms);
        }
    }
    ptype.resize(numTotalAtoms);
    if (opts->ngtc > 1)
    {
        cTC.resize(numTotalAtoms);
        /* We always copy cTC with domain decomposition */
    }
    cENER.resize(numTotalAtoms);
    if (inputrec_.useConstantAcceleration)
    {
        cACC.resize(numTotalAtoms);
    }
    if (inputrecFrozenAtoms(&inputrec_))
    {
        cFREEZE.resize(numTotalAtoms);
    }
    if (bVCMgrps)
    {
        cVCM.resize(numTotalAtoms);
    }
    if (nPerturbed)
    {
        bPerturbed.resize(numTotalAtoms);
    }

    // Note that these user groups are empty
    // when there is only one group present.
    // Therefore, when adding code, the user should use something like:
    // gprnrU1 = (cU1.empty() ? 0 : cU1[localatindex])
    if (!mtop_.groups.groupNumbers[SimulationAtomGroupType::User1].empty())
    {
        cU1.resize(numTotalAtoms);
    }
    if (!mtop_.groups.groupNumbers[SimulationAtomGroupType::User2].empty())
    {
        cU2.resize(numTotalAtoms);
    }

    MTopLookUp mTopLookUp(mtop_);

    const unsigned short numTypes   = mtop_.ffparams.atnr;
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

            if (!cFREEZE.empty())
            {
                cFREEZE[i] =
                        (isValidAtom ? getGroupType(groups, SimulationAtomGroupType::Freeze, ag) : 0);
            }
            if (EI_ENERGY_MINIMIZATION(inputrec_.eI))
            {
                /* Displacement is proportional to F, masses used for constraints */
                mA = 1.0;
                mB = 1.0;
            }
            else if (inputrec_.eI == IntegrationAlgorithm::BD)
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
                if (inputrec_.bd_fric > 0)
                {
                    mA = 0.5 * inputrec_.bd_fric * inputrec_.delta_t;
                    mB = 0.5 * inputrec_.bd_fric * inputrec_.delta_t;
                }
                else
                {
                    /* The friction coefficient is mass/tau_t */
                    fac = inputrec_.delta_t
                          / opts->tau_t[!cTC.empty() ? groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag] : 0];
                    mA = 0.5 * atom.m * fac;
                    mB = 0.5 * atom.mB * fac;
                }
            }
            else
            {
                mA = atom.m;
                mB = atom.mB;
            }
            if (nMassPerturbed)
            {
                massA[i] = mA;
                massB[i] = mB;
            }
            // When masses are perturbed, all perturbed-atom
            // mass-related quantities are instead assigned in
            // updateMdAtomsForMassLambdas().
            if (!nMassPerturbed or !bPerturbed[i])
            {
                massT[i] = mA;

                if (mA == 0.0)
                {
                    invmass[i]           = 0;
                    invMassPerDim[i][XX] = 0;
                    invMassPerDim[i][YY] = 0;
                    invMassPerDim[i][ZZ] = 0;
                }
                else if (!cFREEZE.empty())
                {
                    g = cFREEZE[i];
                    GMX_ASSERT(opts->nFreeze != nullptr,
                               "Must have freeze groups to initialize masses");
                    if (opts->nFreeze[g][XX] && opts->nFreeze[g][YY] && opts->nFreeze[g][ZZ])
                    {
                        /* Set the mass of completely frozen particles to ALMOST_ZERO
                         * iso 0 to avoid div by zero in lincs or shake.
                         */
                        invmass[i] = ALMOST_ZERO;
                    }
                    else
                    {
                        /* Note: Partially frozen particles use the normal invmass.
                         * If such particles are constrained, the frozen dimensions
                         * should not be updated with the constrained coordinates.
                         */
                        invmass[i] = 1.0 / mA;
                    }
                    for (int d = 0; d < DIM; d++)
                    {
                        invMassPerDim[i][d] = (opts->nFreeze[g][d] ? 0 : 1.0 / mA);
                    }
                }
                else
                {
                    invmass[i] = 1.0 / mA;
                    for (int d = 0; d < DIM; d++)
                    {
                        invMassPerDim[i][d] = 1.0 / mA;
                    }
                }
            }

            chargeA[i] = atom.q;
            typeA[i]   = atom.type;
            if (bLJPME)
            {
                real c6 =
                        (isValidAtom
                                 ? mtop_.ffparams.iparams[atom.type * (mtop_.ffparams.atnr + 1)].lj.c6
                                 : 0.0_real);
                real c12 =
                        (isValidAtom
                                 ? mtop_.ffparams.iparams[atom.type * (mtop_.ffparams.atnr + 1)].lj.c12
                                 : 0.0_real);
                sqrt_c6A[i] = std::sqrt(c6);
                if (c6 == 0.0 || c12 == 0)
                {
                    sigmaA[i] = 1.0;
                }
                else
                {
                    sigmaA[i] = gmx::sixthroot(c12 / c6);
                }
                sigma3A[i] = 1 / (sigmaA[i] * sigmaA[i] * sigmaA[i]);
            }
            if (nPerturbed)
            {
                bPerturbed[i] = PERTURBED(atom);
                chargeB[i]    = atom.qB;
                typeB[i]      = atom.typeB;
                if (bLJPME)
                {
                    real c6     = (isValidAtom ? mtop_.ffparams
                                                     .iparams[atom.typeB * (mtop_.ffparams.atnr + 1)]
                                                     .lj.c6
                                               : 0.0_real);
                    real c12    = (isValidAtom ? mtop_.ffparams
                                                      .iparams[atom.typeB * (mtop_.ffparams.atnr + 1)]
                                                      .lj.c12
                                               : 0.0_real);
                    sqrt_c6B[i] = std::sqrt(c6);
                    if (c6 == 0.0 || c12 == 0)
                    {
                        sigmaB[i] = 1.0;
                    }
                    else
                    {
                        sigmaB[i] = gmx::sixthroot(c12 / c6);
                    }
                    sigma3B[i] = 1 / (sigmaB[i] * sigmaB[i] * sigmaB[i]);
                }
            }
            ptype[i] = atom.ptype;

            if (isValidAtom)
            {
                if (!cTC.empty())
                {
                    cTC[i] = groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ag];
                }
                cENER[i] = getGroupType(groups, SimulationAtomGroupType::EnergyOutput, ag);
                if (!cACC.empty())
                {
                    cACC[i] = groups.groupNumbers[SimulationAtomGroupType::Acceleration][ag];
                }
                if (!cVCM.empty())
                {
                    cVCM[i] = groups.groupNumbers[SimulationAtomGroupType::MassCenterVelocityRemoval][ag];
                }

                if (!cU1.empty())
                {
                    cU1[i] = groups.groupNumbers[SimulationAtomGroupType::User1][ag];
                }
                if (!cU2.empty())
                {
                    cU2[i] = groups.groupNumbers[SimulationAtomGroupType::User2][ag];
                }
            }
            else
            {
                // As fillers have no mass and interactions, we can add them to group 0 without side-effects
                if (!cTC.empty())
                {
                    cTC[i] = 0;
                }
                cENER[i] = 0;
                if (!cACC.empty())
                {
                    cACC[i] = 0;
                }
                if (!cVCM.empty())
                {
                    cVCM[i] = 0;
                }
                if (!cU1.empty())
                {
                    cU1[i] = -1;
                }
                if (!cU2.empty())
                {
                    cU2[i] = -1;
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    if (numTotalAtoms > 0)
    {
        /* Pad invmass with 0 so a SIMD MD update does not change v and x */
        for (int i = numTotalAtoms; i < invmass.paddedSize(); i++)
        {
            invmass[i] = 0;
        }
    }

    numHomeAtoms = numHomeAtomsArg;
    /* We set mass, invmass, invMassPerDim and tmass for lambda=0.
     * For free-energy runs, these should be updated using updateMdAtomsForMassLambda().
     */
    tmass = tmassA;
    updateForMassLambda(massLambda);
}

void MDAtoms::updateForMassLambda(const real massLambda)
{
    if (nMassPerturbed)
    {
        real L1 = 1 - massLambda;

        const int numTotalAtoms = massT.size();

        /* Update masses of perturbed atoms for the change in mass lambda */
        int gmx_unused nthreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0; i < numTotalAtoms; i++)
        {
            if (bPerturbed[i])
            {
                massT[i] = L1 * massA[i] + massLambda * massB[i];
                /* Atoms with invmass 0 or ALMOST_ZERO are massless or frozen
                 * and their invmass does not depend on lambda.
                 */
                if (invmass[i] > 1.1 * ALMOST_ZERO)
                {
                    invmass[i] = 1.0 / massT[i];
                    for (int d = 0; d < DIM; d++)
                    {
                        if (invMassPerDim[i][d] > 1.1 * ALMOST_ZERO)
                        {
                            invMassPerDim[i][d] = invmass[i];
                        }
                    }
                }
            }
        }

        /* Update the system mass for the change in mass lambda */
        tmass = L1 * tmassA + massLambda * tmassB;
    }
}

} // namespace gmx
