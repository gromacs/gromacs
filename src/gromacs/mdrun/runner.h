/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017,2018, by the GROMACS development team, led by
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
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_RUNNER_H
#define GMX_MDRUN_RUNNER_H

#include <cstdio>

#include <array>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
struct ReplicaExchangeParameters;
struct t_commrec;

namespace gmx
{

/*!
 * \brief Create the default set of MD filename options.
 *
 * \return Ownership of a new filename option container.
 */
std::unique_ptr<std::array<t_filenm, 34>> makeDefaultMdFilenames();

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
 * \ingroup module_mdrun
 */
class Mdrunner
{
    private:
        //! Parallelism-related user options.
        gmx_hw_opt_t             hw_opt;
        //! Filenames and properties from command-line argument values.

        std::unique_ptr<std::array<t_filenm, 34>> filenames{nullptr};

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
        //! Parameters for replica-exchange simulations.
        ReplicaExchangeParameters        replExParams;
        //! Print a warning if any force is larger than this (in kJ/mol nm).
        real                             pforce = -1;
        //! Handle to file used for logging.
        FILE                            *fplog;
        //! Handle to communication data structure.
        t_commrec                       *cr;
        //! Handle to multi-simulation handler.
        gmx_multisim_t                  *ms;

    public:
        /*! \brief Builder class to manage object creation.
         *
         * This class is a member of gmx::Mdrunner to allow access to private
         * gmx::Mdrunner members.
         *
         * It is non-trivial to establish an initialized gmx::Mdrunner invariant,
         * so objects can be obtained by clients using a Builder, a move, or a
         * clone() operation. Clients cannot default initialize or copy
         * gmx::Mdrunner.
         */
        class BuilderImplementation;

        /*! \brief Constructor.
         *
         * Note that when member variables are not present in the constructor
         * member initialization list (which is true for the default constructor),
         * then they are initialized with any default member initializer specified
         * when they were declared, or default initialized. */
        Mdrunner() = default;

        ~Mdrunner();

        // Copy requires special attention. Use clone methods.
        Mdrunner(const Mdrunner &)            = delete;
        Mdrunner &operator=(const Mdrunner &) = delete;

        // Allow move
//        Mdrunner(Mdrunner &&) noexcept;
//        Mdrunner &operator=(Mdrunner &&) noexcept;

        /*! \brief Driver routine, that calls the different simulation methods. */
        /*!
         * Currently, thread-MPI does not spawn threads until during mdrunner() and parallelism
         * is not initialized until some time during this call...
         */
        int mdrunner();
        //! Called when thread-MPI spawns threads.
        t_commrec *spawnThreads(int numThreadsToLaunch) const;
        /*! \brief Re-initializes the object after threads spawn.
         *
         * \todo Can this be refactored so that the Mdrunner on a spawned thread is
         * constructed ready to use? */
        // Replace with cloneOnSpawnedThread to get a new ready-to-use Mdrunner on the new thread.
//        void reinitializeOnSpawnedThread();

        /*! \brief Initializes a new Mdrunner from the master.
         *
         * Run in a new thread from a const pointer to the master.
         * \returns New Mdrunner instance suitable for running in additional threads.
         */
        std::unique_ptr<Mdrunner> cloneOnSpawnedThread() const;

};

/*! \libinternal
 * \brief Build a gmx::Mdrunner.
 *
 * Client code (such as the `mdrun` CLI program) uses this builder to get an initialized Mdrunner.
 *
 * A builder allows the library to ensure that client code cannot obtain an
 * uninitialized or partially initialized runner by refusing to build() if the
 * client has not provided sufficient or self-consistent direction. Director
 * code can be implemented for different user interfaces, encapsulating any
 * run-time functionality that does not belong in the library MD code, such
 * as command-line option processing or interfacing to external libraries.
 *
 * \ingroup module_mdrun
 *
 * \internal
 *
 * The initial Builder implementation is neither extensible at run time nor
 * at compile time. Future implementations should evolve to compose the runner,
 * rather than just consolidating the parameters for initialization, but there
 * is not yet a firm design for how flexibly module code will be coupled to
 * the builder and how much of the client interface will be in this Builder
 * versus Builders provided by the various modules.
 *
 * As the modules are more clearly encapsulated, these `set` methods should be
 * replaced with `add` methods to attach objects for different modules,
 * and each module can provide its own builder and user interface helpers to
 * the Director code. The runner and client code will also have to be updated
 * as appropriate default behavior is clarified for (a) default behavior of
 * client when user does not provide input, (b) default behavior of builder
 * when client does not provide input, and (c) default behavior of runner
 * when builder does not provide input.
 */
class MdrunnerBuilder final
{
    public:
        /*!
         * \brief A default constructor is available to client code in the GROMACS build environment.
         */
        MdrunnerBuilder();

        /*!
         * \brief Get ownership of an initialized gmx::Mdrunner.
         *
         * After build() is called, the Builder object should not be used
         * again. More clearly defined behavior requires updates to data
         * ownership of input arguments and Mdrunner members.
         *
         * Some input parameters have not been wrapped in `set` or `add` methods
         * yet because their use should probably be modernized and/or encapsulated
         * into appropriate functional modules (rather than being passed directly).
         * These are run-time inputs that can trigger automatically-generated
         * run-time parameters. Further clarification is needed to establish
         * whether the final parameter determined by the runner is the result
         * of direct user input, a client default, a runner default, or a value
         * that was calculated automatically, either at the request of the user
         * or as a default behavior.
         *
         * \param nbpu_opt Target short-range interations for "cpu", "gpu", or "auto".
         * \param pme_opt Target long-range interactions for "cpu", "gpu", or "auto".
         * \param pme_fft_opt Target long-range interactions FFT/solve stages for "cpu", "gpu", or "auto".

         * \return ownership of a new Mdrunner.
         */
        std::unique_ptr<Mdrunner> build(const char* nbpu_opt,
                                        const char* pme_opt,
                                        const char* pme_fft_opt);

        /*!
         * \brief Provide access to the commrec to use.
         *
         * \param communicator non-owning handle to comm pointer.
         * \return
         */
        MdrunnerBuilder &setCommunications(t_commrec** communicator);

        /*!
         * \brief Provide access to the multisim communicator to use.
         *
         * \param multisim non-owning handle to multisim comm pointer.
         * \return
         */
        MdrunnerBuilder &addMultiSim(gmx_multisim_t** multisim);

        /*!
         * \brief Provide access to the output environment resources to use.
         *
         * \param outputEnvironment non-owning handle to output context.
         * \param logFile non-owning handle to log filehandle.
         * \return
         */
        MdrunnerBuilder &setOutputContext(gmx_output_env_t** outputEnvironment, FILE** logFile);

        /*!
         * \brief Set Mdrun options not owned by some other module.
         *
         * \param options structure to copy
         * \param forceWarningThreshold Print a warning if any force is larger than this (in kJ/mol nm)
         * \return
         */
        MdrunnerBuilder &setExtraMdrunOptions(const MdrunOptions &options,
                                              real                forceWarningThreshold);

        /*!
         * \brief Set the domain decomposition module.
         *
         * \param options options with which to construct domain decomposition.
         * \return
         */
        MdrunnerBuilder &setDomdec(const DomdecOptions &options);

        /*!
         * \brief Set parallelism resource management.
         *
         * \param options parallelism options to copy.
         * \return
         */
        MdrunnerBuilder &setHardwareOptions(const gmx_hw_opt_t &options);

        /*!
         * \brief Set Verlet list manager.
         *
         * \param rebuildInterval override for the duration of a neighbor list with the Verlet scheme.
         * \return
         */
        MdrunnerBuilder &setVerletList(int rebuildInterval);

        /*!
         * \brief Set replica exchange manager.
         *
         * \param params parameters with which to set up replica exchange.
         * \return
         */
        MdrunnerBuilder &setReplicaExchange(const ReplicaExchangeParameters &params);

        /*!
         * \brief Set I/O files.
         *
         * Borrow a container of t_filenm. Note that std::array<t_filenm> is
         * not copyable, moveable, or assignable. Each element must be uniquely
         * created in one place, but some members of the t_filenm elements are
         * updated in various scopes.
         *
         * It is the responsibility of the calling code to make sure that the
         * pointer remains valid for the lifetime of the MdrunnerBuilder and
         * its product, the Runner.
         *
         * \param filenames Borrowed pointer to container of t_filename.
         */
        MdrunnerBuilder &setFilenames(std::unique_ptr<std::array<t_filenm, 34>> filenames);

        ~MdrunnerBuilder();


    private:
        std::unique_ptr<Mdrunner::BuilderImplementation> impl_;
};

}      // namespace gmx

#endif // GMX_MDRUN_RUNNER_H
