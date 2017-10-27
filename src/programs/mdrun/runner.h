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
 */
class Mdrunner
{
    private:
        //! Parallelism-related user options.
        gmx_hw_opt_t             hw_opt;
        //! Filenames and properties from command-line argument values.
#ifdef BUILD_WITH_FDA
        std::array<t_filenm, 42 > filenames =
#else
        std::array<t_filenm, 34 > filenames =
#endif
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
#ifdef BUILD_WITH_FDA
          { efNDX, "-pfn",      nullptr,     ffOPTRD },
          { efPFI, "-pfi",      "fda",       ffOPTRD },
          { efPFA, "-pfa",      "fda",       ffOPTWR },
          { efPFR, "-pfr",      "fda",       ffOPTWR },
          { efPSA, "-psa",      "fda",       ffOPTWR },
          { efPSR, "-psr",      "fda",       ffOPTWR },
          { efVSA, "-vsa",      "fda",       ffOPTWR },
          { efVMA, "-vma",      "fda",       ffOPTWR },
#endif
          { efXVG, "-ro",       "rotation",  ffOPTWR },
          { efLOG, "-ra",       "rotangles", ffOPTWR },
          { efLOG, "-rs",       "rotslabs",  ffOPTWR },
          { efLOG, "-rt",       "rottorque", ffOPTWR },
          { efMTX, "-mtx",      "nm",        ffOPTWR },
          { efRND, "-multidir", nullptr,     ffOPTRDMULT},
          { efXVG, "-awh",      "awhinit",   ffOPTRD },
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
        //! Target long-range interactions for "cpu", "gpu", or "auto". Default is "auto".
        const char                      *pme_opt = nullptr;
        //! Target long-range interactions FFT/solve stages for "cpu", "gpu", or "auto". Default is "auto".
        const char                      *pme_fft_opt = nullptr;
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
