/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Declares the checkpoint handler class.
 *
 * This class sets the signal to checkpoint based on the elapsed simulation time,
 * and handles the signal if it is received. When handling the signal (via
 * decideIfCheckpointingThisStep()), it is deciding whether a checkpoint should be saved
 * at the current step. This can be due to a received signal, or if the current simulation
 * step is the last. This information can be queried via the isCheckpointingStep() function.
 *
 * The setting and handling is implemented in private functions. They are only called
 * if a respective boolean is true. For the trivial case of no checkpointing (or no checkpoint
 * signal setting on any other rank than master), the translation unit of the calling
 * function is therefore never left. In the future, this will be achieved by adding
 * (or not adding) handlers / setters to the task graph.
 *
 * TODO: Update this to add that this is no actually reading and writing checkpoints
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_CHECKPOINTHANDLER_H
#define GMX_MDLIB_CHECKPOINTHANDLER_H

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/icheckpointclient.h"
#include "gromacs/mdrunutility/accumulateglobals.h"

struct DomdecOptions;
struct gmx_file_position_t;
struct gmx_walltime_accounting;
struct ObservablesHistory;
struct t_commrec;
struct t_fileio;
struct t_inputrec;
class t_state;
struct t_trxframe;

namespace gmx
{
class CheckpointHelper;

/*! \brief Checkpoint signals
 *
 * Signals set and read by CheckpointHandler. Possible signals include
 *   * nothing to signal
 *   * do checkpoint (at next NS step)
 */
enum class CheckpointSignal
{
    noSignal, doCheckpoint
};

/*! \libinternal
 * \brief Class handling the checkpoint signal
 *
 * Master rank sets the checkpointing signal periodically
 * All ranks receive checkpointing signal and set the respective flag
 */
class CheckpointHandler final : public IAccumulateGlobalsClient
{
    public:
        /*! \brief CheckpointHandler constructor
         *
         * Needs a non-null pointer to the builder to register for global communication, and
         * (const) references to data it needs to determine whether a signal needs to be
         * set or handled.
         */
        CheckpointHandler(
            compat::not_null<AccumulateGlobalsBuilder*> accumulateGlobalsBuilder,
            bool                                        neverUpdateNeighborList,
            bool                                        isMaster,
            bool                                        writeFinalCheckpoint,
            real                                        checkpointingPeriod,
            std::map<CheckpointKeyword, compat::not_null<ICheckpointClient*> > clients,
            const std::string &writeCheckpointFilename,
            int simulationPart);

        /*! \brief Decides whether a checkpointing signal needs to be set
         *
         * Checkpointing signal is set based on the elapsed run time and the checkpointing
         * interval.
         */
        void setSignal(gmx_walltime_accounting *walltime_accounting) const
        {
            if (rankCanSetSignal_)
            {
                setSignalImpl(walltime_accounting);
            }
        }

        /*! \brief Checks whether a signal has been received
         *
         */
        void handleSignal()
        {
            if (checkpointingIsActive_ &&
                static_cast<CheckpointSignal>(lround(signalView_[0])) == CheckpointSignal::doCheckpoint)
            {
                signalView_[0]      = static_cast<double>(CheckpointSignal::noSignal);
                hasUnhandledSignal_ = true;
            }
        }

        /*! \brief Decides whether a checkpoint shall be written at this step
         *
         * Checkpointing is done if this is not the initial step, and
         *   * a signal has been set and the current step is a neighborlist creation
         *     step, or
         *   * the current step is the last step and a the simulation is writing
         *     configurations.
         */
        void decideIfCheckpointingThisStep(bool bNS, bool bFirstStep, bool bLastStep)
        {
            if (checkpointingIsActive_)
            {
                decideIfCheckpointingThisStepImpl(bNS, bFirstStep, bLastStep);
            }
        }

        //! Query decision in decideIfCheckpointingThisStep()
        bool isCheckpointingStep() const
        {
            return checkpointThisStep_;
        }

        /* From IAccumulateGlobalsClient */
        //! Return the number of values needed to pass signal.
        int getNumGlobalsRequired() const override;
        //! Store where to write and read signal
        void setViewForGlobals(AccumulateGlobals *accumulateGlobals,
                               ArrayRef<double>   view) override;
        //! Called (in debug mode) after MPI reduction is complete.
        void notifyAfterCommunication() override;

        //! Read checkpoint
        void readCheckpoint(
            std::string readCheckpointFilename, FILE **pfplog,
            const t_commrec *cr, const DomdecOptions &domdecOptions, t_inputrec *ir,
            t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory,
            bool bAppendOutputFiles, bool bForceAppend, bool reproducibilityRequested);

        //! Write checkpoint
        void writeCheckpoint(
            int64_t step, double time, const t_commrec *cr,
            bool bExpanded, int elamstats, int eIntegrator,
            t_state *state, ObservablesHistory *observablesHistory,
            bool bNumberAndKeep, FILE *fplog);

        /* Read everything that can be stored in t_trxframe from a checkpoint file */
        static void readCheckpointTrxframe(
            t_fileio   *fp,
            t_trxframe *fr);

        /* Print the complete contents of checkpoint file fn to out */
        static void listCheckpoint(
            std::string readCheckpointFilename,
            FILE       *out);

        /* ! \brief Read simulation step and part from a checkpoint file
         *
         * Used by tune_pme to handle tuning with a checkpoint file as part of the input.
         *
         * \param[in]  readCheckpointFilename  Name of checkpoint file
         * \param[out] simulationPart        The part of the simulation that wrote the checkpoint
         * \param[out] step                   The final step number of the simulation that wrote the checkpoint
         *
         * The output variables will both contain 0 if filename is NULL, the file
         * does not exist, or is not readable. */
        static void readCheckpointPartAndStep(
            std::string  readCheckpointFilename,
            int         *simulationPart,
            int64_t     *step);

        /* ! \brief Read simulation part and output filenames from a checkpoint file
         *
         * Used by mdrun to handle restarts
         *
         * \param[in]  readCheckpointFilename  Name of checkpoint file
         * \param[out] simulationPart         The part of the simulation that wrote the checkpoint
         * \param[out] outputfiles             Container of output file names from the previous run.
         */
        static void readCheckpointSimulationPartAndFilenames(
            std::string                       readCheckpointFilename,
            int                              *simulationPart,
            std::vector<gmx_file_position_t> *outputfiles);

    private:
        void setSignalImpl(gmx_walltime_accounting *walltime_accounting) const;

        void decideIfCheckpointingThisStepImpl(bool bNS, bool bFirstStep, bool bLastStep);

        void readCheckpointImpl(
            compat::not_null<CheckpointHelper*> checkpointHelper, FILE **pfplog,
            const t_commrec *cr, const DomdecOptions &domdecOptions, t_inputrec *ir,
            t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory,
            bool bAppendOutputFiles, bool bForceAppend, bool reproducibilityRequested);

        ArrayRef<double>  signalView_;

        bool              hasUnhandledSignal_;
        bool              checkpointThisStep_;
        int               numberOfNextCheckpoint_;

        const bool        rankCanSetSignal_;
        const bool        checkpointingIsActive_;
        const bool        writeFinalCheckpoint_;
        const bool        neverUpdateNeighborlist_;
        const real        checkpointingPeriod_;

        std::map<CheckpointKeyword, compat::not_null<ICheckpointClient*> > clientMap_;
        std::vector<int>     intData_;
        std::vector<int64_t> int64Data_;
        std::vector<real>    realData_;
        std::vector<double>  doubleData_;
        int                  numClients_;

        const std::string    writeCheckpointFilename_;
        int                  simulationPart_;
};

/*! \libinternal
 * \brief
 *
 *
 */
class CheckpointHandlerBuilder final
{
    public:
        /*! \brief Register checkpoint client
         *
         */
        void registerCheckpointClient(compat::not_null<ICheckpointClient*> client);

        /* \brief Create CheckpointHandler
         *
         */
        std::unique_ptr<CheckpointHandler> getCheckpointHandler(
            compat::not_null<AccumulateGlobalsBuilder*> accumulateGlobalsBuilder,
            bool                                        neverUpdateNeighborList,
            bool                                        isMaster,
            bool                                        writeFinalCheckpoint,
            real                                        checkpointingPeriod,
            const std::string                          &writeCheckpointFilename,
            int                                         simulationPart);
    private:
        /*! \brief Registered checkpoint clients
         *
         */
        std::map<CheckpointKeyword, compat::not_null<ICheckpointClient*> > checkpointClients_;
};
}      // namespace gmx

#endif // GMX_MDLIB_CHECKPOINTHANDLER_H
