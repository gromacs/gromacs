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
/*! \libinternal \file
 *
 * \brief This file declares helper functionality for legacy option handling for mdrun
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 * \inlibraryapi
 */
#ifndef GMX_MDRUN_LEGACYMDRUNOPTIONS_H
#define GMX_MDRUN_LEGACYMDRUNOPTIONS_H

#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/domdec/options.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/real.h"

#include "replicaexchange.h"

struct gmx_output_env_t;

namespace gmx
{
template<typename>
class ArrayRef;

/*! \libinternal
 * \brief This class provides the same command-line option
 * functionality to both CLI and API sessions.
 *
 * This class should not exist, but is necessary now to introduce
 * support for the CLI and API without duplicating code. It should be
 * eliminated following the TODOs below.
 *
 * \warning Instances provide lifetime scope for members that do not have
 *  effective lifetime management or which are frequently accessed unsafely.
 *  The caller is responsible for keeping a LegacyMdrunOptions object alive
 *  for as long as any consumers, direct or transitive.
 *
 * \todo Modules in mdrun should acquire proper option handling so
 *       that all of these declarations and defaults are local to the
 *       modules.
 *
 * \todo Contextual aspects, such as working directory
 *       and environment variable handling are more properly
 *       the role of SimulationContext, and should be moved there.
 */
class LegacyMdrunOptions
{
public:
    //! Ongoing collection of mdrun options
    MdrunOptions mdrunOptions;
    //! Options for the domain decomposition.
    DomdecOptions domdecOptions;
    //! Parallelism-related user options.
    gmx_hw_opt_t hw_opt;
    //! Command-line override for the duration of a neighbor list with the Verlet scheme.
    int nstlist_cmdline = 0;
    //! Parameters for replica-exchange simulations.
    ReplicaExchangeParameters replExParams;

    //! Filename options to fill from command-line argument values.
    std::vector<t_filenm> filenames = { { { efTPR, nullptr, nullptr, ffREAD },
                                          { efTRN, "-o", nullptr, ffWRITE },
                                          { efCOMPRESSED, "-x", nullptr, ffOPTWR },
                                          { efCPT, "-cpi", nullptr, ffOPTRD | ffALLOW_MISSING },
                                          { efCPT, "-cpo", nullptr, ffOPTWR },
                                          { efSTO, "-c", "confout", ffWRITE },
                                          { efEDR, "-e", "ener", ffWRITE },
                                          { efLOG, "-g", "md", ffWRITE },
                                          { efXVG, "-dhdl", "dhdl", ffOPTWR },
                                          { efXVG, "-field", "field", ffOPTWR },
                                          { efXVG, "-table", "table", ffOPTRD },
                                          { efXVG, "-tablep", "tablep", ffOPTRD },
                                          { efXVG, "-tableb", "table", ffOPTRDMULT },
                                          { efTRX, "-rerun", "rerun", ffOPTRD },
                                          { efXVG, "-tpi", "tpi", ffOPTWR },
                                          { efXVG, "-tpid", "tpidist", ffOPTWR },
                                          { efEDI, "-ei", "sam", ffOPTRD },
                                          { efXVG, "-eo", "edsam", ffOPTWR },
                                          { efXVG, "-px", "pullx", ffOPTWR },
                                          { efXVG, "-pf", "pullf", ffOPTWR },
                                          { efXVG, "-ro", "rotation", ffOPTWR },
                                          { efLOG, "-ra", "rotangles", ffOPTWR },
                                          { efLOG, "-rs", "rotslabs", ffOPTWR },
                                          { efLOG, "-rt", "rottorque", ffOPTWR },
                                          { efMTX, "-mtx", "nm", ffOPTWR },
                                          { efRND, "-multidir", nullptr, ffOPTRDMULT },
                                          { efXVG, "-awh", "awhinit", ffOPTRD },
                                          { efDAT, "-membed", "membed", ffOPTRD },
                                          { efTOP, "-mp", "membed", ffOPTRD },
                                          { efNDX, "-mn", "membed", ffOPTRD },
                                          { efXVG, "-if", "imdforces", ffOPTWR },
                                          { efXVG, "-swap", "swapions", ffOPTWR } } };

    //! Print a warning if any force is larger than this (in kJ/mol nm).
    real pforce = -1;

    //! The value of the -append option
    bool appendOption = true;

    /*! \brief Output context for writing text files
     *
     * \todo Clarify initialization, ownership, and lifetime. */
    gmx_output_env_t* oenv = nullptr;

    /*! \brief Command line options, defaults, docs and storage for them to fill. */
    /*! \{ */
    rvec        realddxyz                                                           = { 0, 0, 0 };
    const char* ddrank_opt_choices[static_cast<int>(DdRankOrder::Count) + 1]        = { nullptr,
                                                                                 "interleave",
                                                                                 "pp_pme",
                                                                                 "cartesian",
                                                                                 nullptr };
    const char* dddlb_opt_choices[static_cast<int>(DlbOption::Count) + 1]           = { nullptr,
                                                                              "auto",
                                                                              "no",
                                                                              "yes",
                                                                              nullptr };
    const char* thread_aff_opt_choices[static_cast<int>(ThreadAffinity::Count) + 1] = { nullptr,
                                                                                        "auto",
                                                                                        "on",
                                                                                        "off",
                                                                                        nullptr };
    const char* nbpu_opt_choices[5]    = { nullptr, "auto", "cpu", "gpu", nullptr };
    const char* pme_opt_choices[5]     = { nullptr, "auto", "cpu", "gpu", nullptr };
    const char* pme_fft_opt_choices[5] = { nullptr, "auto", "cpu", "gpu", nullptr };
    const char* bonded_opt_choices[5]  = { nullptr, "auto", "cpu", "gpu", nullptr };
    const char* update_opt_choices[5]  = { nullptr, "auto", "cpu", "gpu", nullptr };
    const char* devicesSelectedByUser  = "";
    const char* userGpuTaskAssignment  = "";


    ImdOptions& imdOptions = mdrunOptions.imdOptions;

    t_pargs pa[48] = {

        { "-dd", FALSE, etRVEC, { &realddxyz }, "Domain decomposition grid, 0 is optimize" },
        { "-ddorder", FALSE, etENUM, { ddrank_opt_choices }, "DD rank order" },
        { "-npme",
          FALSE,
          etINT,
          { &domdecOptions.numPmeRanks },
          "Number of separate ranks to be used for PME, -1 is guess" },
        { "-nt",
          FALSE,
          etINT,
          { &hw_opt.nthreads_tot },
          "Total number of threads to start (0 is guess)" },
        { "-ntmpi",
          FALSE,
          etINT,
          { &hw_opt.nthreads_tmpi },
          "Number of thread-MPI ranks to start (0 is guess)" },
        { "-ntomp",
          FALSE,
          etINT,
          { &hw_opt.nthreads_omp },
          "Number of OpenMP threads per MPI rank to start (0 is guess)" },
        { "-ntomp_pme",
          FALSE,
          etINT,
          { &hw_opt.nthreads_omp_pme },
          "Number of OpenMP threads per MPI rank to start (0 is -ntomp)" },
        { "-pin",
          FALSE,
          etENUM,
          { thread_aff_opt_choices },
          "Whether mdrun should try to set thread affinities" },
        { "-pinoffset",
          FALSE,
          etINT,
          { &hw_opt.core_pinning_offset },
          "The lowest logical core number to which mdrun should pin the first thread" },
        { "-pinstride",
          FALSE,
          etINT,
          { &hw_opt.core_pinning_stride },
          "Pinning distance in logical cores for threads, use 0 to minimize the number of threads "
          "per physical core" },
        { "-gpu_id",
          FALSE,
          etSTR,
          { &devicesSelectedByUser },
          "List of unique GPU device IDs available to use" },
        { "-gputasks",
          FALSE,
          etSTR,
          { &userGpuTaskAssignment },
          "List of GPU device IDs, mapping each task on a node to a device. "
          "Tasks include PP and PME (if present)." },
        { "-ddcheck",
          FALSE,
          etBOOL,
          { &domdecOptions.ddBondedChecking },
          "Check for all bonded interactions with DD" },
        { "-ddbondcomm",
          FALSE,
          etBOOL,
          { &domdecOptions.useBondedCommunication },
          "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
        { "-rdd",
          FALSE,
          etREAL,
          { &domdecOptions.minimumCommunicationRange },
          "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial "
          "coordinates" },
        { "-rcon",
          FALSE,
          etREAL,
          { &domdecOptions.constraintCommunicationRange },
          "Maximum distance for P-LINCS (nm), 0 is estimate" },
        { "-dlb", FALSE, etENUM, { dddlb_opt_choices }, "Dynamic load balancing (with DD)" },
        { "-dds",
          FALSE,
          etREAL,
          { &domdecOptions.dlbScaling },
          "Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in "
          "order to "
          "provide a margin in which dynamic load balancing can act while preserving the minimum "
          "cell size." },
        { "-ddcsx",
          FALSE,
          etSTR,
          { &domdecOptions.cellSizeX },
          "HIDDENA string containing a vector of the relative sizes in the x "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsy",
          FALSE,
          etSTR,
          { &domdecOptions.cellSizeY },
          "HIDDENA string containing a vector of the relative sizes in the y "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsz",
          FALSE,
          etSTR,
          { &domdecOptions.cellSizeZ },
          "HIDDENA string containing a vector of the relative sizes in the z "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-nb", FALSE, etENUM, { nbpu_opt_choices }, "Calculate non-bonded interactions on" },
        { "-nstlist",
          FALSE,
          etINT,
          { &nstlist_cmdline },
          "Set nstlist when using a Verlet buffer tolerance (0 is guess)" },
        { "-tunepme",
          FALSE,
          etBOOL,
          { &mdrunOptions.tunePme },
          "Optimize PME load between PP/PME ranks or GPU/CPU" },
        { "-pme", FALSE, etENUM, { pme_opt_choices }, "Perform PME calculations on" },
        { "-pmefft", FALSE, etENUM, { pme_fft_opt_choices }, "Perform PME FFT calculations on" },
        { "-bonded", FALSE, etENUM, { bonded_opt_choices }, "Perform bonded calculations on" },
        { "-update", FALSE, etENUM, { update_opt_choices }, "Perform update and constraints on" },
        { "-v", FALSE, etBOOL, { &mdrunOptions.verbose }, "Be loud and noisy" },
        { "-pforce", FALSE, etREAL, { &pforce }, "Print all forces larger than this (kJ/mol nm)" },
        { "-reprod",
          FALSE,
          etBOOL,
          { &mdrunOptions.reproducible },
          "Avoid optimizations that affect binary reproducibility; "
          "this can significantly reduce performance" },
        { "-cpt",
          FALSE,
          etREAL,
          { &mdrunOptions.checkpointOptions.period },
          "Checkpoint interval (minutes)" },
        { "-cpnum",
          FALSE,
          etBOOL,
          { &mdrunOptions.checkpointOptions.keepAndNumberCheckpointFiles },
          "Keep and number checkpoint files" },
        { "-append",
          FALSE,
          etBOOL,
          { &appendOption },
          "Append to previous output files when continuing from checkpoint instead of adding the "
          "simulation part number to all file names" },
        { "-nsteps",
          FALSE,
          etINT64,
          { &mdrunOptions.numStepsCommandline },
          "Run this number of steps (-1 means infinite, -2 means use mdp option, smaller is "
          "invalid)" },
        { "-maxh",
          FALSE,
          etREAL,
          { &mdrunOptions.maximumHoursToRun },
          "Terminate after 0.99 times this time (hours)" },
        { "-replex",
          FALSE,
          etINT,
          { &replExParams.exchangeInterval },
          "Attempt replica exchange periodically with this period (steps)" },
        { "-nex",
          FALSE,
          etINT,
          { &replExParams.numExchanges },
          "Number of random exchanges to carry out each exchange interval (N^3 is one suggestion). "
          " -nex zero or not specified gives neighbor replica exchange." },
        { "-reseed",
          FALSE,
          etINT,
          { &replExParams.randomSeed },
          "Seed for replica exchange, -1 is generate a seed" },
        { "-imdport", FALSE, etINT, { &imdOptions.port }, "HIDDENIMD listening port" },
        { "-imdwait",
          FALSE,
          etBOOL,
          { &imdOptions.wait },
          "HIDDENPause the simulation while no IMD client is connected" },
        { "-imdterm",
          FALSE,
          etBOOL,
          { &imdOptions.terminatable },
          "HIDDENAllow termination of the simulation from IMD client" },
        { "-imdpull",
          FALSE,
          etBOOL,
          { &imdOptions.pull },
          "HIDDENAllow pulling in the simulation from IMD client" },
        { "-rerunvsite",
          FALSE,
          etBOOL,
          { &mdrunOptions.rerunConstructVsites },
          "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
        { "-confout",
          FALSE,
          etBOOL,
          { &mdrunOptions.writeConfout },
          "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last "
          "step" },
        { "-stepout",
          FALSE,
          etINT,
          { &mdrunOptions.verboseStepPrintInterval },
          "HIDDENFrequency of writing the remaining wall clock time for the run" },
        { "-resetstep",
          FALSE,
          etINT,
          { &mdrunOptions.timingOptions.resetStep },
          "HIDDENReset cycle counters after these many time steps" },
        { "-resethway",
          FALSE,
          etBOOL,
          { &mdrunOptions.timingOptions.resetHalfway },
          "HIDDENReset the cycle counters after half the number of steps or halfway "
          "[TT]-maxh[tt]" }
    };
    /*! \} */

    //! Parses the command-line input and prepares to start mdrun.
    int updateFromCommandLine(int argc, char** argv, ArrayRef<const char*> desc);

    ~LegacyMdrunOptions();
};

} // end namespace gmx

#endif
