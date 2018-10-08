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
/*! \internal \file
 * \brief
 * Defines the checkpoint handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "checkpointhandler.h"

#include "buildinfo.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/checkpointhelper.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

namespace gmx
{

// TODO: Put that in the right spot - for now here to make our life easier
/*! \brief Name of checkpointing keywords
 *
 */
static std::string getCheckpointKeywordName(CheckpointKeyword keyword)
{
    static const std::map<CheckpointKeyword, std::string> checkpointKeywordNames = {
        {CheckpointKeyword::state, "t_state"},
        {CheckpointKeyword::ekinstate, "ekinstate_t"},
        {CheckpointKeyword::swaphist, "swaphistory_t"},
        {CheckpointKeyword::enerhist, "energyhistory_t"},
        {CheckpointKeyword::dfhist, "df_history_t"},
        {CheckpointKeyword::edsamhist, "edsamhistory_t"},
        {CheckpointKeyword::correlationGridHistory, "CorrelationGridHistory"},
        {CheckpointKeyword::awhBiasHistory, "AwhBiasHistory"},
        {CheckpointKeyword::awhHistory, "AwhHistory"}
    };

    return checkpointKeywordNames.at(keyword);
};

CheckpointHandler::CheckpointHandler(
        compat::not_null<AccumulateGlobalsBuilder*> accumulateGlobalsBuilder,
        bool                                        neverUpdateNeighborList,
        bool                                        isMaster,
        bool                                        writeFinalCheckpoint,
        real                                        checkpointingPeriod,
        std::map<CheckpointKeyword, compat::not_null<ICheckpointClient*> > clients,
        const std::string &writeCheckpointFilename,
        int simulationPart) :
    hasUnhandledSignal_(false),
    checkpointThisStep_(false),
    numberOfNextCheckpoint_(1),
    rankCanSetSignal_(checkpointingPeriod >= 0 && isMaster),
    checkpointingIsActive_(checkpointingPeriod >= 0),
    writeFinalCheckpoint_(writeFinalCheckpoint),
    neverUpdateNeighborlist_(neverUpdateNeighborList),
    checkpointingPeriod_(checkpointingPeriod),
    clientMap_(clients),
    numClients_(static_cast<int>(clients.size())),
    writeCheckpointFilename_(writeCheckpointFilename),
    simulationPart_(simulationPart)
{
    if (checkpointingIsActive_)
    {
        // In multi-sim settings, checkpointing signal needs to be shared across simulations
        accumulateGlobalsBuilder->registerClient(
                compat::not_null<IAccumulateGlobalsClient*>(this), true);
    }
}

void CheckpointHandler::setSignalImpl(
        gmx_walltime_accounting_t walltime_accounting) const
{
    if (!hasUnhandledSignal_ &&
        (checkpointingPeriod_ == 0 ||
         walltime_accounting_get_time_since_start(walltime_accounting) >=
         numberOfNextCheckpoint_ * checkpointingPeriod_ * 60.0))
    {
        signalView_[0] = static_cast<double>(CheckpointSignal::doCheckpoint);
    }
}

void CheckpointHandler::decideIfCheckpointingThisStepImpl(
        bool bNS, bool bFirstStep, bool bLastStep)
{
    if (((hasUnhandledSignal_ && (bNS || neverUpdateNeighborlist_)) ||
         (bLastStep && writeFinalCheckpoint_)) &&
        !bFirstStep)
    {
        checkpointThisStep_ = true;
        hasUnhandledSignal_ = false;
        numberOfNextCheckpoint_++;
    }
}

int CheckpointHandler::getNumGlobalsRequired() const
{
    return 1;
}

void CheckpointHandler::setViewForGlobals(
        AccumulateGlobals *accumulateGlobals gmx_unused,
        ArrayRef<double>   view)
{
    // No need for AccumulateGlobals pointer

    signalView_ = view;
    // set signal empty
    signalView_[0] = static_cast<double>(CheckpointSignal::noSignal);
}

void CheckpointHandler::notifyAfterCommunication()
{
    // Not sure if we'll need that for anything...
}

void CheckpointHandler::writeCheckpoint(
        int64_t step, double time, const t_commrec *cr,
        bool bExpanded, int elamstats, int eIntegrator,
        t_state *state, ObservablesHistory *observablesHistory,
        bool bNumberAndKeep, FILE *fplog)
{
    char timebuf[STRLEN];
    gmx_format_current_time(timebuf, STRLEN);

    if (fplog)
    {
        fprintf(fplog, "Writing checkpoint, step %s at %s\n\n",
                std::to_string(step).c_str(), timebuf);
    }

    // Create helper & write header
    const ivec       ones             = {1, 1, 1};
    CheckpointHelper checkpointHelper = CheckpointHelper(
                gmx_version(), GMX_DOUBLE, getProgramContext().fullBinaryPath(),
                timebuf, BUILD_TIME, BUILD_USER, BUILD_HOST,
                DOMAINDECOMP(cr) ? cr->dd->nnodes : cr->nnodes,
                DOMAINDECOMP(cr) ? cr->npmenodes : 0,
                DOMAINDECOMP(cr) ? cr->dd->nc : ones,
                simulationPart_, step, time,
                gmx_fio_get_output_file_positions(), writeCheckpointFilename_);

    auto filePointer = checkpointHelper.openWriteCheckpointFile(step);

    checkpointHelper.doHeader(compat::make_unique<XDRWrapper>(filePointer), false);

    if (checkpointHelper.doFiles(compat::make_unique<XDRWrapper>(filePointer), false))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    // Add client info to header
    auto xdrWrapper = compat::make_unique<XDRWrapper>(filePointer);
    xdrWrapper->doCptIntErr(&numClients_);
    for (const auto &entry : clientMap_)
    {
        const auto &client  = entry.second;
        auto        keyword = static_cast<int>(client->getKeyword());
        auto        version = client->getVersion();
        xdrWrapper->doCptIntErr(&keyword);
        xdrWrapper->doCptIntErr(&version);
        client->notifyWrite();
    }

    // Write int / real / double data
    for (const auto &entry : clientMap_)
    {
        const auto &client  = entry.second;
        if (xdrWrapper->doArrayRef(client->getIntView()) ||
            xdrWrapper->doArrayRef(client->getInt64View()) ||
            xdrWrapper->doArrayRef(client->getRealView()) ||
            xdrWrapper->doArrayRef(client->getDoubleView()))
        {
            gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
        }
    }

    /* START LEGACY CODE
     *
     * The code above would do all the work - but for now, the clients are not implemented. The code
     * below uses the legacy functions to achieve the same.
     *
     */

    CheckpointHelper::writeLegacy(
            compat::make_unique<XDRWrapper>(filePointer),
            bExpanded, elamstats, eIntegrator,
            state, observablesHistory);

    // END LEGACY CODE

    // Write footer
    checkpointHelper.doFooter(compat::make_unique<XDRWrapper>(filePointer), false);
    checkpointHelper.closeWriteCheckpointFile(filePointer, bNumberAndKeep);

#if GMX_FAHCORE
    /*code for alternate checkpointing scheme.  moved from top of loop over
       steps */
    fcRequestCheckPoint();
    if (fcCheckPointParallel( cr->nodeid, NULL, 0) == 0)
    {
        gmx_fatal( 3, __FILE__, __LINE__, "Checkpoint error on step %d\n", step );
    }
#endif  /* end GMX_FAHCORE block */
}

void CheckpointHandler::readCheckpoint(
        std::string readCheckpointFilename, FILE **pfplog,
        const t_commrec *cr, const DomdecOptions &domdecOptions, t_inputrec *ir,
        t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory,
        bool bAppendOutputFiles, bool bForceAppend, bool reproducibilityRequested)
{
    bool isLegacy = false;
    if (SIMMASTER(cr))
    {
        isLegacy = CheckpointHelper::getCheckpointVersion(readCheckpointFilename) == CheckpointVersion::legacy;
    }
    gmx_bcast(sizeof(isLegacy), &isLegacy, cr);
    if (isLegacy)
    {
        legacy::load_checkpoint(
                readCheckpointFilename.c_str(), pfplog, cr, domdecOptions.numCells,
                ir, state, bReadEkin, observablesHistory,
                bAppendOutputFiles, bForceAppend, reproducibilityRequested);
        return;
    }

    auto checkpointHelper = CheckpointHelper(readCheckpointFilename);
    if (SIMMASTER(cr))
    {
        readCheckpointImpl(
                compat::not_null<CheckpointHelper*>(&checkpointHelper), pfplog,
                cr, domdecOptions, ir, state, bReadEkin, observablesHistory,
                bAppendOutputFiles, bForceAppend, reproducibilityRequested);

        // TODO: Copied from legacy checkpoint - clean that up!
        ir->fepvals->init_fep_state = state->fep_state;  /* there should be a better way to do this than setting it here.
                                                        Investigate for 5.0. */
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(checkpointHelper.step_), &checkpointHelper.step_, cr);
        gmx_bcast(sizeof(*bReadEkin), bReadEkin, cr);
    }
    ir->bContinuation    = TRUE;
    if (ir->nsteps >= 0)
    {
        ir->nsteps          += ir->init_step - checkpointHelper.step_;
    }
    ir->init_step        = checkpointHelper.step_;
    ir->simulation_part  = checkpointHelper.simulationPart_ + 1;

    this->simulationPart_ = checkpointHelper.simulationPart_ + 1;
}


void CheckpointHandler::readCheckpointImpl(
        compat::not_null<CheckpointHelper*> checkpointHelper, FILE **pfplog,
        const t_commrec *cr, const DomdecOptions &domdecOptions, t_inputrec *ir,
        t_state *state, gmx_bool *bReadEkin, ObservablesHistory *observablesHistory,
        bool bAppendOutputFiles, bool bForceAppend, bool reproducibilityRequested)
{
    FILE *fplog = *pfplog;

    // Create helper & read header
    auto filePointer = checkpointHelper->openReadCheckpointFile();
    checkpointHelper->doHeader(compat::make_unique<XDRWrapper>(filePointer), false);

    if (bAppendOutputFiles && checkpointHelper->doublePrecision_ != GMX_DOUBLE)
    {
        gmx_fatal(FARGS, "Output file appending requested, but the code and "
                  "checkpoint file precision (single/double) don't match");
    }

    /* This will not be written if we do appending, since fplog is still NULL then */
    // TODO: More elegant way to check this?
    if (fplog)
    {
        checkpointHelper->printHeader(fplog);
    }

    if (MASTER(cr))
    {
        checkpointHelper->checkMatch(
                fplog, reproducibilityRequested, cr, domdecOptions);
    }

    if (checkpointHelper->doFiles(compat::make_unique<XDRWrapper>(filePointer), true))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    // Get client info from header
    auto xdrWrapper = compat::make_unique<XDRWrapper>(filePointer);
    int  fileNumClients;
    xdrWrapper->doCptIntErr(&fileNumClients);
    if (fileNumClients != numClients_)
    {
        gmx_fatal(FARGS, "Mismatch in number of checkpointing clients. File has %d, while current "
                  "simulation has %d.", fileNumClients, numClients_);
    }
    for (int i = 0; i < fileNumClients; ++i)
    {
        int intKeyword, version;
        xdrWrapper->doCptIntErr(&intKeyword);
        xdrWrapper->doCptIntErr(&version);
        auto keyword = static_cast<CheckpointKeyword>(intKeyword);
        if (!clientMap_.count(keyword))
        {
            gmx_fatal(FARGS, "Checkpointing client mismatch. Client %s found in file is not registered"
                      "in current simulation.", getCheckpointKeywordName(keyword).c_str());
        }
        if (clientMap_.at(keyword)->getVersion() != version)
        {
            gmx_fatal(FARGS, "Version mismatch for checkpointing client %s. File has %d, while current "
                      "simulation has %d.", getCheckpointKeywordName(keyword).c_str(), version, numClients_);
        }
    }

    // Read int / real / double data
    for (const auto &entry : clientMap_)
    {
        const auto &client  = entry.second;
        if (xdrWrapper->doArrayRef(client->getIntView()) ||
            xdrWrapper->doArrayRef(client->getInt64View()) ||
            xdrWrapper->doArrayRef(client->getRealView()) ||
            xdrWrapper->doArrayRef(client->getDoubleView()))
        {
            gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
        }
    }

    for (const auto &entry : clientMap_)
    {
        const auto &client = entry.second;
        client->notifyRead();
    }

    /* START LEGACY CODE
     *
     * The code above would do all the work - but for now, the clients are not implemented. The code
     * below uses the legacy functions to achieve the same.
     *
     */

    CheckpointHelper::readLegacy(
            compat::make_unique<XDRWrapper>(filePointer), writeCheckpointFilename_, ir,
            state, bReadEkin, observablesHistory);

    // END LEGACY CODE

    checkpointHelper->doFooter(compat::make_unique<XDRWrapper>(filePointer), true);

    checkpointHelper->closeReadCheckpointFile(filePointer);

    if (bAppendOutputFiles)
    {
        checkpointHelper->appendOutputFiles(fplog, pfplog, bForceAppend);
    }
}

void CheckpointHandler::readCheckpointTrxframe(
        t_fileio   *fp,
        t_trxframe *fr)
{
    bool isLegacy =
        CheckpointHelper::getCheckpointVersion(gmx_fio_getname(fp)) == CheckpointVersion::legacy;

    if (isLegacy)
    {
        legacy::read_checkpoint_trxframe(fp, fr);
    }
    else
    {
        // TODO: Implement this
        gmx_fatal(FARGS, "Reading checkpoint file to trx not implemented for new "
                  "checkpoint format.");
    }
}

void CheckpointHandler::listCheckpoint(
        std::string readCheckpointFilename,
        FILE       *out)
{
    bool isLegacy =
        CheckpointHelper::getCheckpointVersion(readCheckpointFilename) == CheckpointVersion::legacy;

    if (isLegacy)
    {
        legacy::list_checkpoint(readCheckpointFilename.c_str(), out);
    }
    else
    {
        // TODO: Implement this
        gmx_fatal(FARGS, "Listing checkpoint file not implemented for new "
                  "checkpoint format.");
    }
}

void CheckpointHandler::readCheckpointPartAndStep(
        std::string readCheckpointFilename,
        int        *simulationPart,
        int64_t    *step)
{
    bool isLegacy =
        CheckpointHelper::getCheckpointVersion(readCheckpointFilename) == CheckpointVersion::legacy;

    if (isLegacy)
    {
        legacy::read_checkpoint_part_and_step(
                readCheckpointFilename.c_str(), simulationPart, step);
    }
    else
    {
        CheckpointHelper checkpointHelper = CheckpointHelper(readCheckpointFilename);
        auto             fp               = checkpointHelper.openReadCheckpointFile();
        checkpointHelper.doHeader(compat::make_unique<XDRWrapper>(fp), true);
        checkpointHelper.closeReadCheckpointFile(fp);

        *simulationPart = checkpointHelper.simulationPart_;
        *step           = checkpointHelper.step_;
    }
}

void CheckpointHandler::readCheckpointSimulationPartAndFilenames(
        std::string                       readCheckpointFilename,
        int                              *simulationPart,
        std::vector<gmx_file_position_t> *outputfiles)
{
    bool             isLegacy =
        CheckpointHelper::getCheckpointVersion(readCheckpointFilename) == CheckpointVersion::legacy;
    CheckpointHelper checkpointHelper = CheckpointHelper(readCheckpointFilename);
    auto             fp               = checkpointHelper.openReadCheckpointFile();
    if (isLegacy)
    {
        legacy::read_checkpoint_simulation_part_and_filenames(
                fp, simulationPart, outputfiles);
    }
    else
    {
        checkpointHelper.doHeader(compat::make_unique<XDRWrapper>(fp), true);
        checkpointHelper.doFiles(compat::make_unique<XDRWrapper>(fp), true);
        *simulationPart = checkpointHelper.simulationPart_;
        // TODO: Check this actually works...
        *outputfiles = std::move(checkpointHelper.outputfiles_);
    }
    checkpointHelper.closeReadCheckpointFile(fp);
}

void CheckpointHandlerBuilder::registerCheckpointClient(gmx::compat::not_null<ICheckpointClient *> client)
{
    auto module = client->getKeyword();
    GMX_ASSERT(static_cast<CheckpointKeyword>(0) <= module && module < CheckpointKeyword::count, "Invalid module");
    GMX_RELEASE_ASSERT(
            checkpointClients_.insert(
                    std::pair<CheckpointKeyword, compat::not_null<ICheckpointClient*> >(module, client)).second,
            "Multiple clients writing to the same keyword.");
}

std::unique_ptr<CheckpointHandler> CheckpointHandlerBuilder::getCheckpointHandler(
        compat::not_null<AccumulateGlobalsBuilder*> accumulateGlobalsBuilder,
        bool                                        neverUpdateNeighborList,
        bool                                        isMaster,
        bool                                        writeFinalCheckpoint,
        real                                        checkpointingPeriod,
        const std::string                          &writeCheckpointFilename,
        int                                         simulationPart)
{
    return compat::make_unique<CheckpointHandler>(
            accumulateGlobalsBuilder, neverUpdateNeighborList,
            isMaster, writeFinalCheckpoint, checkpointingPeriod,
            checkpointClients_, writeCheckpointFilename, simulationPart);
}

} // namespace gmx
