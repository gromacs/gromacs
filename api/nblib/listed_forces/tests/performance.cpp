/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <chrono>
#include <iostream>
#include <numeric>

#include "gromacs/mdlib/calcvir.h"
#include "gromacs/topology/topology.h"

#include "nblib/listed_forces/aggregate_transformations.hpp"
#include "nblib/listed_forces/calculator.h"
#include "nblib/listed_forces/conversionscommon.h"
#include "nblib/listed_forces/traits.h"
#include "nblib/tpr.h"

#include "gmxcalculator.h"
#include "linear_chain_input.hpp"
#include "poly_dimethyl_butene_input.hpp"

namespace nblib
{
namespace test
{
namespace
{

std::vector<real> compareForces(gmx::ArrayRef<const gmx::RVec> probe, gmx::ArrayRef<const gmx::RVec> ref)
{
    using T = real;
    if (probe.size() != ref.size())
    {
        throw InputException("compareForces input mismatch\n");
    }

    int            numParticles = probe.size();
    std::vector<T> relErr(numParticles);

    for (size_t i = 0; i < numParticles; ++i)
    {
        T num   = norm(probe[i] - ref[i]);
        T denom = norm(ref[i]);
        T error;
        if (denom == T(0))
        {
            error = (num == T(0)) ? T(0) : INFINITY;
        }
        else
        {
            error = num / denom;
        }

        relErr[i] = error;
    }

    std::sort(begin(relErr), end(relErr));

    return relErr;
}

void aggregationEfficiency(ListedInteractionData interactions)
{
    using AggregateType = TypeListElement_t<0, AggregateTypes>;

    real bondsB4 = pickType<typename AggregateType::TwoCenterAggregateType>(interactions).indices.size();
    real anglesB4 =
            pickType<typename AggregateType::ThreeCenterAggregateType>(interactions).indices.size();
    real pairsB4 = pickType<typename AggregateType::PairAggregateType>(interactions).indices.size();

    auto t0 = std::chrono::high_resolution_clock::now();
    createAggregates(interactions);
    auto t1 = std::chrono::high_resolution_clock::now();
    printf("migration time %f ms\n", std::chrono::duration<double, std::milli>(t1 - t0).count());

    real bondsAfter =
            pickType<typename AggregateType::TwoCenterAggregateType>(interactions).indices.size();
    real anglesAfter =
            pickType<typename AggregateType::ThreeCenterAggregateType>(interactions).indices.size();
    real pairsAfter = pickType<typename AggregateType::PairAggregateType>(interactions).indices.size();

    real bondsRatio  = (bondsB4 > 0) ? 1 - bondsAfter / bondsB4 : 0;
    real anglesRatio = (anglesB4 > 0) ? 1 - anglesAfter / anglesB4 : 0;
    real pairsRatio  = (pairsB4 > 0) ? 1 - pairsAfter / pairsB4 : 0;

    printf("aggregate integration efficiency\nbonds: %f, angles: %f, pairs: %f\n",
           bondsRatio,
           anglesRatio,
           pairsRatio);
}

//! \brief compare whether the forces between a pair of gmx/nblib datasets match up
void compareForces(const ListedInteractionData&   nblibInteractions,
                   const InteractionDefinitions&  idefs,
                   const gmx_ffparams_t&          ffparams,
                   gmx::ArrayRef<const gmx::RVec> x,
                   gmx::ArrayRef<const real>      q,
                   Box                            box,
                   int                            numThreads)
{
    ListedForceCalculator nblibCalculator(nblibInteractions, x.size(), numThreads, box);
    ListedGmxCalculator   gmxCalculator(idefs, ffparams, x.size(), numThreads, box);

    ListedEnergies         nblibEnergies;
    std::vector<gmx::RVec> nblibForces(x.size(), gmx::RVec{ 0, 0, 0 });
    std::vector<gmx::RVec> nblibShiftForces(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 });
    nblibCalculator.compute(x, q, nblibForces, nblibShiftForces, nblibEnergies);

    ListedEnergies         gmxEnergies;
    std::vector<gmx::RVec> gmxForces(x.size(), Vec3{ 0, 0, 0 });
    std::vector<gmx::RVec> gmxShiftForces(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 });
    gmxCalculator.compute(x, q, gmxForces, gmxShiftForces, gmxEnergies);

    auto errors = compareForces(nblibForces, gmxForces);
    std::cout << "max relative force error: " << errors.back() << std::endl;
    auto shiftErrors = compareForces(nblibShiftForces, gmxShiftForces);
    std::cout << "max relative shift force difference: " << shiftErrors.back() << std::endl;

    t_forcerec forcerec;
    forcerec.shift_vec.resize(gmx::c_numShiftVectors);
    calc_shifts(box.legacyMatrix(), forcerec.shift_vec);
    std::vector<gmx::RVec> nblibVirial(3, gmx::RVec{ 0, 0, 0 });
    std::vector<gmx::RVec> gmxVirial(3, gmx::RVec{ 0, 0, 0 });
    calc_vir(gmx::c_numShiftVectors,
             gmx::as_rvec_array(forcerec.shift_vec.data()),
             gmx::as_rvec_array(nblibShiftForces.data()),
             gmx::as_rvec_array(nblibVirial.data()),
             false,
             box.legacyMatrix());
    calc_vir(gmx::c_numShiftVectors,
             gmx::as_rvec_array(forcerec.shift_vec.data()),
             gmx::as_rvec_array(gmxShiftForces.data()),
             gmx::as_rvec_array(gmxVirial.data()),
             false,
             box.legacyMatrix());

    auto virialErrors = compareForces(nblibVirial, gmxVirial);
    std::cout << "max relative virial error: " << virialErrors.back() << std::endl;
}

//! \brief benchmark fixture for bonded calculators
class listedBenchmarkGmx
{
public:
    listedBenchmarkGmx(const InteractionDefinitions&  idefs,
                       const gmx_ffparams_t&          ffparams,
                       gmx::ArrayRef<const gmx::RVec> x,
                       gmx::ArrayRef<const real>      q,
                       Box                            box,
                       int                            reps_,
                       int                            nThr) :
        coordinates_(x.begin(), x.end()),
        charges_(q.begin(), q.end()),
        shiftForces_(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 }),
        numThreads(nThr),
        reps(reps_),
        forceBuffer_(x.size()),
        lfCalculator_(idefs, ffparams, x.size(), numThreads, box),
        calculatorEnergies{ 0 }
    {
        std::fill(calculatorEnergies.begin(), calculatorEnergies.end(), 0);
    }

    void operator()()
    {
        for (int i = 0; i < reps; ++i)
        {
            lfCalculator_.compute(coordinates_, charges_, forceBuffer_, shiftForces_, calculatorEnergies);
        }
    }

    const ListedEnergies& energies() const { return calculatorEnergies; }

private:
    int numThreads;
    int reps;

    std::vector<gmx::RVec>        coordinates_;
    std::vector<real>             charges_;
    std::vector<gmx::RVec>        shiftForces_;
    gmx::ArrayRef<std::nullptr_t> noShifts_;

    nblib::ListedGmxCalculator lfCalculator_;
    std::vector<Vec3>          forceBuffer_;
    ListedEnergies             calculatorEnergies;
};

//! \brief benchmark fixture for bonded calculators
class listedBenchmarkNblib
{
public:
    listedBenchmarkNblib(const ListedInteractionData&   interactions,
                         gmx::ArrayRef<const gmx::RVec> x,
                         gmx::ArrayRef<const real>      q,
                         Box                            box,
                         int                            reps_,
                         int                            nThr) :
        coordinates_(x.begin(), x.end()),
        charges_(q.begin(), q.end()),
        shiftForces_(gmx::c_numShiftVectors, gmx::RVec{ 0, 0, 0 }),
        numThreads(nThr),
        reps(reps_),
        forceBuffer_(x.size()),
        lfCalculator_(interactions, x.size(), numThreads, box),
        calculatorEnergies{ 0 }
    {
        std::fill(calculatorEnergies.begin(), calculatorEnergies.end(), 0);
    }


    void operator()()
    {
        for (int i = 0; i < reps; ++i)
        {
            lfCalculator_.compute(coordinates_, charges_, forceBuffer_, shiftForces_, calculatorEnergies);
        }
    }

    const ListedEnergies& energies() const { return calculatorEnergies; }

private:
    int numThreads;
    int reps;

    std::vector<gmx::RVec>        coordinates_;
    std::vector<real>             charges_;
    std::vector<gmx::RVec>        shiftForces_;
    gmx::ArrayRef<std::nullptr_t> noShifts_;

    nblib::ListedForceCalculator lfCalculator_;
    std::vector<Vec3>            forceBuffer_;
    ListedEnergies               calculatorEnergies;
};


template<class F>
double timeit(F&& functionToTime)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    functionToTime();
    auto t2 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double, std::milli>(t2 - t1).count();
}

void testNblib(const ListedInteractionData&   interactions,
               gmx::ArrayRef<const gmx::RVec> x,
               gmx::ArrayRef<const real>      q,
               Box                            box,
               int                            reps,
               int                            nThreads)
{
    listedBenchmarkNblib runner(interactions, x, q, box, reps, nThreads);

    [[maybe_unused]] double warmup  = timeit(runner);
    double                  elapsed = timeit(runner);

    auto energies = runner.energies();

    printf("nblib time elapsed %f ms:\n", elapsed);
    printf("  HarmonicBond %f\n", energies[FindIndex<HarmonicBondType, AllListedTypes>{}]);
    printf("  HarmonicAngle %f\n", energies[FindIndex<HarmonicAngle, AllListedTypes>{}]);
    printf("  ProperDihedral %f\n", energies[FindIndex<ProperDihedral, AllListedTypes>{}]);
    printf("  ImproperDihedral %f\n", energies[FindIndex<ImproperDihedral, AllListedTypes>{}]);
    printf("  VanDerWaals %f\n", energies[VdwIndex{}]);
    printf("  Coulomb %f\n", energies[CoulombIndex{}]);
}

void testGmx(const InteractionDefinitions&  idefs,
             const gmx_ffparams_t&          ffparams,
             gmx::ArrayRef<const gmx::RVec> x,
             gmx::ArrayRef<const real>      q,
             Box                            box,
             int                            reps,
             int                            nThreads)
{
    listedBenchmarkGmx runner(idefs, ffparams, x, q, box, reps, nThreads);

    [[maybe_unused]] double warmup  = timeit(runner);
    double                  elapsed = timeit(runner);

    auto energies = runner.energies();

    printf("gmx time elapsed %f ms:\n", elapsed);
    printf("  HarmonicBond %f\n", energies[FindIndex<HarmonicBondType, AllListedTypes>{}]);
    printf("  HarmonicAngle %f\n", energies[FindIndex<HarmonicAngle, AllListedTypes>{}]);
    printf("  ProperDihedral %f\n", energies[FindIndex<ProperDihedral, AllListedTypes>{}]);
    printf("  ImproperDihedral %f\n", energies[FindIndex<ImproperDihedral, AllListedTypes>{}]);
    printf("  VanDerWaals %f\n", energies[VdwIndex{}]);
    printf("  Coulomb %f\n", energies[CoulombIndex{}]);
}

} // namespace
} // namespace test
} // namespace nblib

using namespace nblib::test;
using namespace nblib;

int main(int argc, char* argv[])
{
    const std::string filepath   = argv[1];
    int               numthreads = std::stoi(argv[2]);
    nblib::TprReader  tpr(filepath);
    std::cout << tpr.coordinates_.size() << std::endl;

    ListedInteractionData interactions = convertToNblibInteractions(*tpr.interactionDefinitions_);
    auto [idef, ffparams]              = convertToGmxInteractions(interactions);

    ListedInteractionData aggregated = interactions;
    aggregationEfficiency(aggregated);

    auto printStats = [](auto& param) {
        std::cout << param.indices.size() << " " << param.parametersA.size() << " "
                  << param.parametersB.size() << " "
                  << typeid(std::decay_t<decltype(param.parametersA[0])>).name() << std::endl;
    };

    testNblib(interactions, tpr.coordinates_, tpr.charges_, tpr.getBox(), 1, numthreads);
    testGmx(*idef, *ffparams, tpr.coordinates_, tpr.charges_, tpr.getBox(), 1, numthreads);
    compareForces(interactions,
                  *tpr.interactionDefinitions_,
                  *tpr.ffparams_,
                  tpr.coordinates_,
                  tpr.charges_,
                  tpr.getBox(),
                  numthreads);
    std::cout << std::endl;

    int                    numParticles = 5003;
    nblib::LinearChainData linearChain(numParticles);
    std::vector<real>      zeroQ(linearChain.x.size(), 0);
    auto [linear_idef, linear_ffparams] = convertToGmxInteractions(linearChain.interactions);

    testNblib(linearChain.interactions, linearChain.x, zeroQ, linearChain.box, 100, 4);
    testGmx(*linear_idef, *linear_ffparams, linearChain.x, zeroQ, linearChain.box, 100, 4);
    compareForces(linearChain.interactions,
                  *linear_idef,
                  *linear_ffparams,
                  linearChain.x,
                  zeroQ,
                  linearChain.box,
                  4);
    std::cout << std::endl;

    // Note: this test case needs double precision for forces to match precisely
    nblib::PolyDimethylButene polyData(100);
    zeroQ.resize(polyData.x.size(), 0);
    auto [poly_idef, poly_ffparams] = convertToGmxInteractions(polyData.interactions);

    // run benchmarks, correctness only visible in the potential energy
    testNblib(polyData.interactions, polyData.x, zeroQ, polyData.box, 1, 4);
    testGmx(*poly_idef, *poly_ffparams, polyData.x, zeroQ, polyData.box, 1, 4);
    compareForces(polyData.interactions, *poly_idef, *poly_ffparams, polyData.x, zeroQ, polyData.box, 4);
}
