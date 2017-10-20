/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief Declares the routine running the inetgrators.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_RUNNER_H
#define GMX_MDLIB_RUNNER_H

#include <cstdio>

#include <array>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "repl_ex.h"

struct gmx_output_env_t;
struct ReplicaExchangeParameters;
struct t_commrec;

namespace gmx
{

/*! \libinternal \brief Runner object for supporting setup and execution of mdrun.
 *
 * This class has responsibility for the lifetime of data structures
 * that exist for the life of the simulation, e.g. for logging and
 * communication.
 *
 * \todo Most of the attributes should be declared by specific modules
 * as command-line options. Accordingly, they do not conform to the
 * naming scheme, because that would make for a lot of noise in the
 * diff, only to have it change again when the options move to their
 * modules.
 *
 * \todo Preparing logging and MPI contexts could probably be a
 * higher-level responsibility, so that an Mdrunner would get made
 * without needing to re-initialize these components (as currently
 * happens always for the master rank, and differently for the spawned
 * ranks with thread-MPI).
 *
 * Interaction diagram for mdrun preparation phase:
 *
 * \msc
 *     Mdrunner,
 *     Options,
 *     MpiRuntime,
 *     MultiSim,
 *     MpiEnv,
 *     Hardware,
 *     TaskStrategy,
 *     Inputrec,
 *     Heuristics,
 *     TmpiSetup,
 *     DutySetup,
 *     DomDec,
 *     TaskAssignment,
 *     NbModule,
 *     PmeModule,
 *     Integrator;
 *
 *     Options box Options [ label="commandline and environment variables" ],
 *     MpiRuntime box MpiRuntime [ label="different for Tmpi" ],
 *     MultiSim box MultiSim [ label="trivial for Tmpi" ],
 *     Hardware box Hardware [ label="detected on each node" ],
 *     TaskStrategy box TaskStrategy [ label="are NB or PME on GPUs" ],
 *     TmpiSetup box TmpiSetup [ label="currently\nget_\nnthreads_\nmpi()" ],
 *     DutySetup box DutySetup [ label="rank NB and/or PME duty" ],
 *     DomDec box DomDec [ label="Also MPMD PME" ],
 *     TaskAssignment box TaskAssignment [ label="responsibility for its node" ];
 *
 *     MpiRuntime => MpiEnv [ label="sets up" ];
 *     Mdrunner => Options [ label="parse\nuser input" ];
 *     Options => MultiSim [ label="sets up" ];
 *     MpiEnv => MultiSim [ label="sets up" ];
 *     Mdrunner => Hardware [ label="detect hardware" ];
 *     Mdrunner => Inputrec [ label="read from .tpr" ];
 *     Options => TaskStrategy [ label="gives manual control of" ];
 *     MultiSim => TaskStrategy [ label="influences" ];
 *     Hardware => TaskStrategy [ label="determines default behaviour of" ];
 *     Options => TmpiSetup [ label="gives manual control of" ];
 *     Hardware => TmpiSetup [ label="determines default behaviour of" ];
 *     TaskStrategy => TmpiSetup [ label="determines default behaviour of" ];
 *     Heuristics => TmpiSetup [ label="determines\ndefault\nbehaviour of" ];
 *     Inputrec => TmpiSetup [ label="influences" ];
 *     TmpiSetup => MpiEnv [ label="sets up" ];
 *     --- [ label="thread-MPI now behaves the same as real MPI" ];
 *     Options => DomDec [ label="gives manual control of" ];
 *     MpiEnv => DomDec [ label="context for" ];
 *     DomDec => MpiEnv [ label="continues to fill" ];
 *     Options => DutySetup [ label="gives manual control of" ];
 *     DomDec => DutySetup [ label="handles default duty assignment" ];
 *     Options => TaskAssignment [ label="gives manual control of" ];
 *     Hardware => TaskAssignment [ label="determines default behaviour of" ];
 *     TaskStrategy => TaskAssignment [ label="determines" ];
 *     DutySetup => TaskAssignment [ label="determines" ];
 *     Options => NbModule [ label="gives manual control of" ];
 *     Hardware => NbModule [ label="chooses default buffer size" ];
 *     TaskAssignment => NbModule [ label="gives device\nIDs to" ];
 *     TaskAssignment => PmeModule [ label="gives device\nIDs to" ];
 *     DomDec => Integrator [ label="used by" ];
 *     NbModule => Integrator [ label="used by" ];
 *     PmeModule => Integrator [ label="used by" ];
 *
 * \endmsc
 */
/*   
 *
 *
 *     runner,
 *     module [ URL="\ref gmx::TrajectoryAnalysisModule" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     runner box module [ label="caller owns runner and module objects" ];
 *     module => data [ label="create (in constructor)" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     runner => runner [ label="parse user input" ];
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     runner => runner [ label="initialize topology\nand selections" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     module => data [ label="initialize" ];
 *     runner => runner [ label="read frame 0" ];
 *     runner => module [ label="initAfterFirstFrame()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAfterFirstFrame()" ];
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => module [ label="analyzeFrame(0)",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     module => data [ label="add data",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     module => data [ label="finishFrame()",
 *                      URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner => data [ label="finishFrameSerial()",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => module [ label="analyzeFrame(1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     ...;
 *     --- [ label="loop over frames ends" ];
 *     runner => module [ label="finishAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::finishAnalysis()" ];
 *     module => data [ label="post-process data" ];
 *     runner => module [ label="writeOutput()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::writeOutput()" ];
 * \endmsc
 */
class Mdrunner
{
    private:
        //! Parallelism-related user options.
        gmx_hw_opt_t             hw_opt;
        //! Filenames and properties from command-line argument values.
        std::array<t_filenm, 34> filenames =
        {{{ efTPR, nullptr,     nullptr,     ffREAD },
          { efTRN, "-o",        nullptr,     ffWRITE },
          { efCOMPRESSED, "-x", nullptr,     ffOPTWR },
          { efCPT, "-cpi",      nullptr,     ffOPTRD | ffALLOW_MISSING },
          { efCPT, "-cpo",      nullptr,     ffOPTWR },
          { efSTO, "-c",        "confout",   ffWRITE },
          { efEDR, "-e",        "ener",      ffWRITE },
          { efLOG, "-g",        "md",        ffWRITE },
          { efXVG, "-dhdl",     "dhdl",      ffOPTWR },
          { efXVG, "-field",    "field",     ffOPTWR },
          { efXVG, "-table",    "table",     ffOPTRD },
          { efXVG, "-tablep",   "tablep",    ffOPTRD },
          { efXVG, "-tableb",   "table",     ffOPTRDMULT },
          { efTRX, "-rerun",    "rerun",     ffOPTRD },
          { efXVG, "-tpi",      "tpi",       ffOPTWR },
          { efXVG, "-tpid",     "tpidist",   ffOPTWR },
          { efEDI, "-ei",       "sam",       ffOPTRD },
          { efXVG, "-eo",       "edsam",     ffOPTWR },
          { efXVG, "-devout",   "deviatie",  ffOPTWR },
          { efXVG, "-runav",    "runaver",   ffOPTWR },
          { efXVG, "-px",       "pullx",     ffOPTWR },
          { efXVG, "-pf",       "pullf",     ffOPTWR },
          { efXVG, "-ro",       "rotation",  ffOPTWR },
          { efLOG, "-ra",       "rotangles", ffOPTWR },
          { efLOG, "-rs",       "rotslabs",  ffOPTWR },
          { efLOG, "-rt",       "rottorque", ffOPTWR },
          { efMTX, "-mtx",      "nm",        ffOPTWR },
          { efRND, "-multidir", nullptr,     ffOPTRDMULT},
          { efDAT, "-membed",   "membed",    ffOPTRD },
          { efTOP, "-mp",       "membed",    ffOPTRD },
          { efNDX, "-mn",       "membed",    ffOPTRD },
          { efXVG, "-if",       "imdforces", ffOPTWR },
          { efXVG, "-swap",     "swapions",  ffOPTWR }}};
        /*! \brief Filename arguments.
         *
         * Provided for compatibility with old C-style code accessing
         * command-line arguments that are file names. */
        t_filenm *fnm = filenames.data();
        /*! \brief Number of filename argument values.
         *
         * Provided for compatibility with old C-style code accessing
         * command-line arguments that are file names. */
        int nfile = filenames.size();
        //! Output context for writing text files
        gmx_output_env_t                *oenv = nullptr;
        //! Ongoing collection of mdrun options
        MdrunOptions                     mdrunOptions;
        //! Options for the domain decomposition.
        DomdecOptions                    domdecOptions;
        //! Target short-range interations for "cpu", "gpu", or "auto". Default is "auto".
        const char                      *nbpu_opt = nullptr;
        //! Command-line override for the duration of a neighbor list with the Verlet scheme.
        int                              nstlist_cmdline = 0;
        //! Number of simulations in multi-simulation set.
        int                              nmultisim = 0;
        //! Parameters for replica-exchange simulations.
        ReplicaExchangeParameters        replExParams;
        //! Print a warning if any force is larger than this (in kJ/mol nm).
        real                             pforce = -1;
        //! Handle to file used for logging.
        FILE                            *fplog;
        //! Handle to communication data structure.
        t_commrec                       *cr;

    public:
        /*! \brief Defaulted constructor.
         *
         * Note that when member variables are not present in the constructor
         * member initialization list (which is true for the default constructor),
         * then they are initialized with any default member initializer specified
         * when they were declared, or default initialized. */
        Mdrunner() = default;
        //! Start running mdrun by calling its C-style main function.
        int mainFunction(int argc, char *argv[]);
        /*! \brief Driver routine, that calls the different simulation methods. */
        int mdrunner();
        //! Called when thread-MPI spawns threads.
        t_commrec *spawnThreads(int numThreadsToLaunch);
        /*! \brief Re-initializes the object after threads spawn.
         *
         * \todo Can this be refactored so that the Mdrunner on a spawned thread is
         * constructed ready to use? */
        void reinitializeOnSpawnedThread();
};

}      // namespace gmx

#endif // GMX_MDLIB_RUNNER_H
