/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * This implements basic nblib tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 */

#include "forcecalculator.h"


//! Sets up and runs the kernel calls
//! TODO Refactor this function to return a handle to dispatchNonbondedKernel
//!      that callers can manipulate directly.
static void setupAndRunInstance(NBKernelSystem          &system,
                                const NBKernelOptions   &options,
                                const bool              &printTimings)
{
    // We set the interaction cut-off to the pairlist cut-off
    interaction_const_t   ic   = setupInteractionConst(options);
    t_nrnb                nrnb = { 0 };
    gmx_enerdata_t        enerd(1, 0);

    gmx::StepWorkload     stepWork;
    stepWork.computeForces = true;
    if (options.computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }

    std::unique_ptr<nonbonded_verlet_t> nbv           = setupNbnxmInstance(options, system);
    const PairlistSet                  &pairlistSet   = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    const gmx::index                    numPairs      = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    gmx_cycles_t                        cycles        = gmx_cycles_read();

    t_forcerec forceRec;
    forceRec.ntype = system.numAtomTypes;
    forceRec.nbfp  = system.nonbondedParameters;
    snew(forceRec.shift_vec, SHIFTS);
    calc_shifts(system.box, forceRec.shift_vec);

    std::vector<gmx::RVec> currentCoords = system.coordinates;
    for (int iter = 0; iter < options.numIterations; iter++)
    {
        // Run the kernel without force clearing
        nbv->dispatchNonbondedKernel(gmx::InteractionLocality::Local,
                                     ic, stepWork, enbvClearFNo, forceRec,
                                     &enerd,
                                     &nrnb);
        // There is one output data structure per thread
        std::vector<nbnxn_atomdata_output_t> nbvAtomsOut = nbv->nbat.get()->out;
        integrateCoordinates(nbvAtomsOut, options, system.box, currentCoords);
    }
    system.coordinates = currentCoords;

    cycles = gmx_cycles_read() - cycles;
    if (printTimings)
    {
        printTimingsOutput(options, system, numPairs, cycles);
    }
}
