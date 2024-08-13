/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 *
 * \brief This file contains the main function for the non-bonded kernel benchmark
 *
 * \author Berk Hess <hess@kth.se>
 */

#include "gmxpre.h"

#include "nonbonded_bench.h"

#include <memory>
#include <vector>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/nbnxm/benchmark/bench_setup.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class CommandLineModuleSettings;


namespace
{

class NonbondedBenchmark : public ICommandLineOptionsModule
{
public:
    NonbondedBenchmark() {}

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

private:
    int                     sizeFactor_ = 1;
    NbnxmKernelBenchOptions benchmarkOptions_;
};

void NonbondedBenchmark::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    std::vector<const char*> desc = {
        "[THISMODULE] runs benchmarks for one or more so-called Nbnxm",
        "non-bonded pair kernels. The non-bonded pair kernels are",
        "the most compute intensive part of MD simulations",
        "and usually comprise 60 to 90 percent of the runtime.",
        "For this reason they are highly optimized and several different",
        "setups are available to compute the same physical interactions.",
        "In addition, there are different physical treatments of Coulomb",
        "interactions and optimizations for atoms without Lennard-Jones",
        "interactions. There are also different physical treatments of",
        "Lennard-Jones interactions, but only a plain cut-off is supported",
        "in this tool, as that is by far the most common treatment.",
        "And finally, while force output is always necessary, energy output",
        "is only required at certain steps. In total there are",
        "12 relevant combinations of options. The combinations double to 24",
        "when two different SIMD setups are supported. These combinations",
        "can be run with a single invocation using the [TT]-all[tt] option.",
        "The behavior of each kernel is affected by caching behavior,",
        "which is determined by the hardware used together with the system size",
        "and the cut-off radius. The larger the number of atoms per thread,",
        "the more L1 cache is needed to avoid L1 cache misses.",
        "The cut-off radius mainly affects the data reuse: a larger cut-off",
        "results in more data reuse and makes the kernel less sensitive to cache",
        "misses.[PAR]",
        "OpenMP parallelization is used to utilize multiple hardware threads",
        "within a compute node. In these benchmarks there is no interaction",
        "between threads, apart from starting and closing a single OpenMP",
        "parallel region per iteration. Additionally, threads interact",
        "through sharing and evicting data from shared caches.",
        "The number of threads to use is set with the [TT]-nt[tt] option.",
        "Thread affinity is important, especially with SMT and shared",
        "caches. Affinities can be set through the OpenMP library using",
        "the GOMP_CPU_AFFINITY environment variable.[PAR]",
        "The benchmark tool times one or more kernels by running them",
        "repeatedly for a number of iterations set by the [TT]-iter[tt]",
        "option. An initial kernel call is done to avoid additional initial",
        "cache misses. Times are recording in cycles read from efficient,",
        "high accuracy counters in the CPU. Note that these often do not",
        "correspond to actual clock cycles. For each kernel, the tool",
        "reports the total number of cycles, cycles per iteration,",
        "and (total and useful) pair interactions per cycle.",
        "Because a cluster pair list is used instead of an atom pair list,",
        "interactions are also computed for some atom pairs that are beyond",
        "the cut-off distance. These pairs are not useful (except for",
        "additional buffering, but that is not of interest here),",
        "only a side effect of the cluster-pair setup. The SIMD 2xMM kernel",
        "has a higher useful pair ratio then the 4xM kernel due to a smaller",
        "cluster size, but a lower total pair throughput.",
        "It is best to run this, or for that matter any, benchmark",
        "with locked CPU clocks, as thermal throttling can significantly",
        "affect performance. If that is not an option, the [TT]-warmup[TT]",
        "option can be used to run initial, untimed iterations to warm up",
        "the processor.[PAR]",
        "The most relevant regime is between 0.1 to 1 millisecond per",
        "iteration. Thus it is useful to run with system sizes that cover",
        "both ends of this regime.[PAR]",
        "The [TT]-simd[tt] and [TT]-table[tt] options select different",
        "implementations to compute the same physics. The choice of these",
        "options should ideally be optimized for the target hardware.",
        "Historically, we only found tabulated Ewald correction to be useful",
        "on 2-wide SIMD or 4-wide SIMD without FMA support. As all modern",
        "architectures are wider and support FMA, we do not use tables by",
        "default. The only exceptions are kernels without SIMD, which only",
        "support tables.",
        "Options [TT]-coulomb[tt], [TT]-combrule[tt] and [TT]-halflj[tt]",
        "depend on the force field and composition of the simulated system.",
        "The optimization of computing Lennard-Jones interactions for only",
        "half of the atoms in a cluster is useful for water, which does not",
        "use Lennard-Jones on hydrogen atoms in most water models.",
        "In the MD engine, any clusters where at most half of the atoms",
        "have LJ interactions will automatically use this kernel.",
        "And finally, the [TT]-energy[tt] option selects the computation",
        "of energies, which are usually only needed infrequently."
    };

    settings->setHelpText(desc);

    static const EnumerationArray<NbnxmBenchMarkKernels, const char*> c_nbnxmSimdStrings = {
        { "auto", "no", "4xm", "2xmm" }
    };
    static const EnumerationArray<NbnxmBenchMarkCombRule, const char*> c_combRuleStrings = {
        { "geometric", "lb", "none" }
    };
    static const EnumerationArray<NbnxmBenchMarkCoulomb, const char*> c_coulombTypeStrings = {
        { "ewald", "reaction-field" }
    };

    options->addOption(
            IntegerOption("size").store(&sizeFactor_).description("The system size is 3000 atoms times this value"));
    options->addOption(
            IntegerOption("nt").store(&benchmarkOptions_.numThreads).description("The number of OpenMP threads to use"));
    options->addOption(EnumOption<NbnxmBenchMarkKernels>("simd")
                               .store(&benchmarkOptions_.nbnxmSimd)
                               .enumValue(c_nbnxmSimdStrings)
                               .description("SIMD type, auto runs all supported SIMD setups or no "
                                            "SIMD when SIMD is not supported"));
    options->addOption(EnumOption<NbnxmBenchMarkCoulomb>("coulomb")
                               .store(&benchmarkOptions_.coulombType)
                               .enumValue(c_coulombTypeStrings)
                               .description("The functional form for the Coulomb interactions"));
    options->addOption(
            BooleanOption("table")
                    .store(&benchmarkOptions_.useTabulatedEwaldCorr)
                    .description("Use lookup table for Ewald correction instead of analytical"));
    options->addOption(EnumOption<NbnxmBenchMarkCombRule>("combrule")
                               .store(&benchmarkOptions_.ljCombinationRule)
                               .enumValue(c_combRuleStrings)
                               .description("The LJ combination rule"));
    options->addOption(BooleanOption("halflj")
                               .store(&benchmarkOptions_.useHalfLJOptimization)
                               .description("Use optimization for LJ on half of the atoms"));
    options->addOption(BooleanOption("energy")
                               .store(&benchmarkOptions_.computeVirialAndEnergy)
                               .description("Compute energies in addition to forces"));
    options->addOption(
            BooleanOption("all").store(&benchmarkOptions_.doAll).description("Run all 12 combinations of options for coulomb, halflj, combrule"));
    options->addOption(RealOption("cutoff")
                               .store(&benchmarkOptions_.pairlistCutoff)
                               .description("Pair-list and interaction cut-off distance"));
    options->addOption(IntegerOption("iter")
                               .store(&benchmarkOptions_.numIterations)
                               .description("The number of iterations for each kernel"));
    options->addOption(IntegerOption("warmup")
                               .store(&benchmarkOptions_.numWarmupIterations)
                               .description("The number of iterations for initial warmup"));
    options->addOption(BooleanOption("cycles")
                               .store(&benchmarkOptions_.cyclesPerPair)
                               .description("Report cycles/pair instead of pairs/cycle"));
    options->addOption(
            BooleanOption("time").store(&benchmarkOptions_.reportTime).description("Report micro-seconds instead of cycles"));
    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Csv)
                               .outputFile()
                               .store(&benchmarkOptions_.outputFile)
                               .defaultBasename("nonbonded-benchmark")
                               .description("Also output results in csv format"));
}

void NonbondedBenchmark::optionsFinished()
{
    // We compute the Ewald coefficient here to avoid a dependency of the Nbnxm on the Ewald module
    const real ewald_rtol          = 1e-5;
    benchmarkOptions_.ewaldcoeff_q = calc_ewaldcoeff_q(benchmarkOptions_.pairlistCutoff, ewald_rtol);
}

int NonbondedBenchmark::run()
{
    bench(sizeFactor_, benchmarkOptions_);

    return 0;
}

} // namespace

const char NonbondedBenchmarkInfo::name[] = "nonbonded-benchmark";
const char NonbondedBenchmarkInfo::shortDescription[] =
        "Benchmarking tool for the non-bonded pair kernels.";

ICommandLineOptionsModulePointer NonbondedBenchmarkInfo::create()
{
    return ICommandLineOptionsModulePointer(std::make_unique<NonbondedBenchmark>());
}

} // namespace gmx
