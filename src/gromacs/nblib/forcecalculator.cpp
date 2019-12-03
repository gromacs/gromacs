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
 * Implements nblib ForceCalculator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "forcecalculator.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"

#include "integrator.h"


namespace nblib {

ForceCalculator::ForceCalculator(NBKernelSystem          &system,
                                 const NBKernelOptions   &options)
                                 : nbKernelSystem_(system), nbKernelOptions_(options)
{}

//! Print timings outputs
void ForceCalculator::printTimingsOutput(const NBKernelOptions &options,
                                         const NBKernelSystem  &system,
                                         const gmx::index      &numPairs,
                                         gmx_cycles_t           cycles)
{
    const gmx::EnumerationArray<BenchMarkKernels, std::string>  kernelNames   = { "auto", "no", "4xM", "2xMM" };
    const gmx::EnumerationArray<BenchMarkCombRule, std::string> combruleNames = { "geom.", "LB", "none" };

    // Generate an, accurate, estimate of the number of non-zero pair interactions
    const real                          atomDensity          = system.coordinates.size()/det(system.box);
    const real                          numPairsWithinCutoff = atomDensity*4.0/3.0*M_PI*std::pow(options.pairlistCutoff, 3);
    const real                          numUsefulPairs       = system.coordinates.size()*0.5*(numPairsWithinCutoff + 1);
#if GMX_SIMD
    if (options.nbnxmSimd != BenchMarkKernels::SimdNo)
    {
        fprintf(stdout, "SIMD width:           %d\n", GMX_SIMD_REAL_WIDTH);
    }
#endif
    fprintf(stdout, "System size:          %zu atoms\n", system.coordinates.size());
    fprintf(stdout, "Cut-off radius:       %g nm\n", options.pairlistCutoff);
    fprintf(stdout, "Number of threads:    %d\n", options.numThreads);
    fprintf(stdout, "Number of iterations: %d\n", options.numIterations);
    fprintf(stdout, "Compute energies:     %s\n",
            options.computeVirialAndEnergy ? "yes" : "no");
    if (options.coulombType != BenchMarkCoulomb::ReactionField)
    {
        fprintf(stdout, "Ewald excl. corr.:    %s\n",
                options.nbnxmSimd == BenchMarkKernels::SimdNo || options.useTabulatedEwaldCorr ? "table" : "analytical");
    }
    printf("\n");

    fprintf(stdout, "Coulomb LJ   comb. SIMD    Mcycles  Mcycles/it.   %s\n",
            options.cyclesPerPair ? "cycles/pair" : "pairs/cycle");
    fprintf(stdout, "                                                total    useful\n");

    fprintf(stdout, "%-7s %-4s %-5s %-4s ",
            options.coulombType == BenchMarkCoulomb::Pme ? "Ewald" : "RF",
            options.useHalfLJOptimization ? "half" : "all",
            combruleNames[options.ljCombinationRule].c_str(),
            kernelNames[options.nbnxmSimd].c_str());

    const double dCycles = static_cast<double>(cycles);
    if (options.cyclesPerPair)
    {
        fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n",
                cycles*1e-6,
                dCycles/options.numIterations*1e-6,
                dCycles/(options.numIterations*numPairs),
                dCycles/(options.numIterations*numUsefulPairs));
    }
    else
    {
        fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n",
                dCycles*1e-6,
                dCycles/options.numIterations*1e-6,
                options.numIterations*numPairs/dCycles,
                options.numIterations*numUsefulPairs/dCycles);
    }
}


//! Sets up and runs the kernel calls
//! TODO Refactor this function to return a handle to dispatchNonbondedKernel
//!      that callers can manipulate directly.
void ForceCalculator::compute(const bool printTimings)
{
    // We set the interaction cut-off to the pairlist cut-off
    interaction_const_t   ic   = setupInteractionConst(nbKernelOptions_);
    t_nrnb                nrnb = { 0 };
    gmx_enerdata_t        enerd(1, 0);

    gmx::StepWorkload     stepWork;
    stepWork.computeForces = true;
    if (nbKernelOptions_.computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }

    std::unique_ptr<nonbonded_verlet_t> nbv           = setupNbnxmInstance(nbKernelOptions_, nbKernelSystem_);
    const PairlistSet                  &pairlistSet   = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    const gmx::index                    numPairs      = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    gmx_cycles_t                        cycles        = gmx_cycles_read();

    t_forcerec forceRec;
    forceRec.ntype = nbKernelSystem_.numAtoms;
    forceRec.nbfp  = nbKernelSystem_.nonbondedParameters;
    snew(forceRec.shift_vec, SHIFTS);
    calc_shifts(nbKernelSystem_.box, forceRec.shift_vec);

    std::vector<gmx::RVec> currentCoords = nbKernelSystem_.coordinates;
    for (int iter = 0; iter < nbKernelOptions_.numIterations; iter++)
    {
        // Run the kernel without force clearing
        nbv->dispatchNonbondedKernel(gmx::InteractionLocality::Local,
                                     ic, stepWork, enbvClearFNo, forceRec,
                                     &enerd,
                                     &nrnb);
        // There is one output data structure per thread
        std::vector<nbnxn_atomdata_output_t> nbvAtomsOut = nbv->nbat.get()->out;
        integrateCoordinates(nbvAtomsOut, nbKernelOptions_, nbKernelSystem_.box, currentCoords);
    }
    nbKernelSystem_.coordinates = currentCoords;

    cycles = gmx_cycles_read() - cycles;
    if (printTimings)
    {
        printTimingsOutput(nbKernelOptions_, nbKernelSystem_, numPairs, cycles);
    }
}

} // namespace nblib
