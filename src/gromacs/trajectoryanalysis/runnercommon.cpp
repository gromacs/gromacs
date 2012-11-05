/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::TrajectoryAnalysisRunnerCommon.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "oenv.h"
#include "rmpbc.h"
#include "smalloc.h"
#include "statutil.h"
#include "tpxio.h"
#include "vec.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/runnercommon.h"
#include "gromacs/trajectoryanalysis/scheduler.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "analysissettings-impl.h"

namespace gmx
{

class TrajectoryAnalysisRunnerCommon::Impl
{
    public:
        // TODO: Change runner to const after fixing problem with counting frames.
        Impl(TrajectoryAnalysisRunnerCommon &runner,
            TrajectoryAnalysisSettings *settings);
        ~Impl();

        void finishTrajectory();

        RoundRobinScheduler scheduler_;
        TrajectoryAnalysisSettings &settings_;
        TopologyInformation     topInfo_;

        bool                    bHelp_;
        bool                    bShowHidden_;
        bool                    bQuiet_;
        //! Name of the trajectory file (empty if not provided).
        std::string             trjfile_;
        //! Name of the topology file (empty if no topology provided).
        std::string             topfile_;
        //! Name of the index file (empty if no index file provided).
        std::string             ndxfile_;
        double                  startTime_;
        double                  endTime_;
        double                  deltaTime_;

        gmx_ana_indexgrps_t    *grps_;
        bool                    bTrajOpen_;
        //! The current frame, or \p NULL if no frame loaded yet.
        t_trxframe             *fr;
        gmx_rmpbc_t             gpbc_;
        //! Used to store the status variable from read_first_frame().
        t_trxstatus            *status_;
        output_env_t            oenv_;
        //! Next frame that will be read in by readNextFrame()
        size_t                  nextFrameNumber_;
};

TrajectoryAnalysisRunnerCommon::Impl::Impl(TrajectoryAnalysisRunnerCommon
        &runner, TrajectoryAnalysisSettings *settings)
    : settings_(*settings),
      bHelp_(false), bShowHidden_(false), bQuiet_(false),
      startTime_(0.0), endTime_(0.0), deltaTime_(0.0),
      grps_(NULL),
      bTrajOpen_(false), fr(NULL), gpbc_(NULL), status_(NULL), oenv_(NULL),
      // For this class, first frame is read before calling readNextFrame(), so
      // next frame number needs to be manually set to 1.
      nextFrameNumber_(1)
{
}


TrajectoryAnalysisRunnerCommon::Impl::~Impl()
{
    if (grps_ != NULL)
    {
        gmx_ana_indexgrps_free(grps_);
    }
    finishTrajectory();
    if (fr)
    {
        // There doesn't seem to be a function for freeing frame data
        sfree(fr->x);
        sfree(fr->v);
        sfree(fr->f);
        sfree(fr);
    }
    if (oenv_ != NULL)
    {
        output_env_done(oenv_);
    }
}


void
TrajectoryAnalysisRunnerCommon::Impl::finishTrajectory()
{
    if (bTrajOpen_)
    {
        close_trx(status_);
        bTrajOpen_ = false;
    }
    if (gpbc_ != NULL)
    {
        gmx_rmpbc_done(gpbc_);
        gpbc_ = NULL;
    }
}

/*********************************************************************
 * TrajectoryAnalysisRunnerCommon
 */

TrajectoryAnalysisRunnerCommon::TrajectoryAnalysisRunnerCommon(
        TrajectoryAnalysisSettings *settings)
    : impl_(new Impl(*this, settings))
{
}


TrajectoryAnalysisRunnerCommon::~TrajectoryAnalysisRunnerCommon()
{
}


void
TrajectoryAnalysisRunnerCommon::initOptions(Options *options)
{
    TrajectoryAnalysisSettings &settings = impl_->settings_;

    // Add options for help.
    options->addOption(BooleanOption("h").store(&impl_->bHelp_)
                           .description("Print help and quit"));
    options->addOption(BooleanOption("hidden").store(&impl_->bShowHidden_)
                           .hidden()
                           .description("Show hidden options"));
    options->addOption(BooleanOption("quiet").store(&impl_->bQuiet_)
                           .hidden()
                           .description("Hide options in normal run"));

    // Add common file name arguments.
    options->addOption(FileNameOption("f")
                           .filetype(eftTrajectory).inputFile()
                           .store(&impl_->trjfile_)
                           .defaultBasename("traj")
                           .description("Input trajectory or single configuration"));
    options->addOption(FileNameOption("s")
                           .filetype(eftTopology).inputFile()
                           .store(&impl_->topfile_)
                           .defaultBasename("topol")
                           .description("Input structure"));
    options->addOption(FileNameOption("n")
                           .filetype(eftIndex).inputFile()
                           .store(&impl_->ndxfile_)
                           .defaultBasename("index")
                           .description("Extra index groups"));
    options->addOption(SelectionFileOption("sf"));

    // Add options for trajectory time control.
    options->addOption(DoubleOption("b").store(&impl_->startTime_).timeValue()
                           .description("First frame (%t) to read from trajectory"));
    options->addOption(DoubleOption("e").store(&impl_->endTime_).timeValue()
                           .description("Last frame (%t) to read from trajectory"));
    options->addOption(DoubleOption("dt").store(&impl_->deltaTime_).timeValue()
                           .description("Only use frame if t MOD dt == first time (%t)"));

    // Add time unit option.
    settings.impl_->timeUnitManager.addTimeUnitOption(options, "tu");

    // Add plot options.
    settings.impl_->plotSettings.addOptions(options);

    // Add common options for trajectory processing.
    if (!settings.hasFlag(TrajectoryAnalysisSettings::efNoUserRmPBC))
    {
        options->addOption(BooleanOption("rmpbc").store(&settings.impl_->bRmPBC)
                               .description("Make molecules whole for each frame"));
    }
    if (!settings.hasFlag(TrajectoryAnalysisSettings::efNoUserPBC))
    {
        options->addOption(BooleanOption("pbc").store(&settings.impl_->bPBC)
                               .description("Use periodic boundary conditions for distance calculation"));
    }
}


void
TrajectoryAnalysisRunnerCommon::scaleTimeOptions(Options *options)
{
    impl_->settings_.impl_->timeUnitManager.scaleTimeOptions(options);
}


bool
TrajectoryAnalysisRunnerCommon::optionsFinished(Options *options)
{
    if (impl_->bHelp_)
    {
        return false;
    }

    impl_->settings_.impl_->plotSettings.setTimeUnit(
            impl_->settings_.impl_->timeUnitManager.timeUnit());

    if (impl_->trjfile_.empty() && impl_->topfile_.empty())
    {
        GMX_THROW(InconsistentInputError("No trajectory or topology provided, nothing to do!"));
    }

    if (options->isSet("b"))
        setTimeValue(TBEGIN, impl_->startTime_);
    if (options->isSet("e"))
        setTimeValue(TEND, impl_->endTime_);
    if (options->isSet("dt"))
        setTimeValue(TDELTA, impl_->deltaTime_);

    return true;
}


void
TrajectoryAnalysisRunnerCommon::initIndexGroups(SelectionCollection *selections)
{
    if (impl_->ndxfile_.empty())
    {
        // TODO: Initialize default selections
        selections->setIndexGroups(NULL);
    }
    else
    {
        gmx_ana_indexgrps_init(&impl_->grps_, NULL, impl_->ndxfile_.c_str());
        selections->setIndexGroups(impl_->grps_);
    }
}


void
TrajectoryAnalysisRunnerCommon::doneIndexGroups(SelectionCollection *selections)
{
    if (impl_->grps_ != NULL)
    {
        selections->setIndexGroups(NULL);
        gmx_ana_indexgrps_free(impl_->grps_);
        impl_->grps_ = NULL;
    }
}


void
TrajectoryAnalysisRunnerCommon::initTopology(SelectionCollection *selections)
{
    const TrajectoryAnalysisSettings &settings = impl_->settings_;
    bool bRequireTop
        = settings.hasFlag(TrajectoryAnalysisSettings::efRequireTop)
          || selections->requiresTopology();
    if (bRequireTop && impl_->topfile_.empty())
    {
        GMX_THROW(InconsistentInputError("No topology provided, but one is required for analysis"));
    }

    // Load the topology if requested.
    if (!impl_->topfile_.empty())
    {
        char  title[STRLEN];

        snew(impl_->topInfo_.top_, 1);
        impl_->topInfo_.bTop_ = read_tps_conf(impl_->topfile_.c_str(), title,
                impl_->topInfo_.top_, &impl_->topInfo_.ePBC_,
                &impl_->topInfo_.xtop_, NULL, impl_->topInfo_.boxtop_, TRUE);
        if (hasTrajectory()
            && !settings.hasFlag(TrajectoryAnalysisSettings::efUseTopX))
        {
            sfree(impl_->topInfo_.xtop_);
            impl_->topInfo_.xtop_ = NULL;
        }
    }

    // Read the first frame if we don't know the maximum number of atoms
    // otherwise.
    int  natoms = -1;
    if (!impl_->topInfo_.hasTopology())
    {
        initFirstFrame();
        natoms = impl_->fr->natoms;
    }
    selections->setTopology(impl_->topInfo_.topology(), natoms);

    /*
    if (impl_->bSelDump)
    {
        gmx_ana_poscalc_coll_print_tree(stderr, impl_->pcc);
        fprintf(stderr, "\n");
    }
    */
}


void
TrajectoryAnalysisRunnerCommon::initFirstFrame()
{
    // Return if we have already initialized the trajectory.
    if (impl_->fr)
    {
        return;
    }

    time_unit_t time_unit
        = static_cast<time_unit_t>(impl_->settings_.timeUnit() + 1);
    output_env_init(&impl_->oenv_, 0, NULL, time_unit, FALSE, exvgNONE, 0, 0);

    int frflags = impl_->settings_.frflags();
    frflags |= TRX_NEED_X;

    snew(impl_->fr, 1);

    const TopologyInformation &top = impl_->topInfo_;
    if (hasTrajectory())
    {
        if (!read_first_frame(impl_->oenv_, &impl_->status_,
                              impl_->trjfile_.c_str(), impl_->fr, frflags))
        {
            GMX_THROW(FileIOError("Could not read coordinates from trajectory"));
        }
        impl_->bTrajOpen_ = true;

        if (top.hasTopology() && impl_->fr->natoms > top.topology()->atoms.nr)
        {
            GMX_THROW(InconsistentInputError(formatString(
                      "Trajectory (%d atoms) does not match topology (%d atoms)",
                      impl_->fr->natoms, top.topology()->atoms.nr)));
        }
        // Check index groups if they have been initialized based on the topology.
        /*
        if (top)
        {
            for (int i = 0; i < impl_->sel->nr(); ++i)
            {
                gmx_ana_index_check(impl_->sel->sel(i)->indexGroup(),
                                    impl_->fr->natoms);
            }
        }
        */
    }
    else
    {
        // Prepare a frame from topology information.
        // TODO: Initialize more of the fields.
        if (frflags & (TRX_NEED_V))
        {
            GMX_THROW(NotImplementedError("Velocity reading from a topology not implemented"));
        }
        if (frflags & (TRX_NEED_F))
        {
            GMX_THROW(InvalidInputError("Forces cannot be read from a topology"));
        }
        impl_->fr->flags  = frflags;
        impl_->fr->natoms = top.topology()->atoms.nr;
        impl_->fr->bX     = TRUE;
        snew(impl_->fr->x, impl_->fr->natoms);
        memcpy(impl_->fr->x, top.xtop_,
               sizeof(*impl_->fr->x) * impl_->fr->natoms);
        impl_->fr->bBox   = TRUE;
        copy_mat(const_cast<rvec *>(top.boxtop_), impl_->fr->box);
    }

    impl_->scheduler_.init();
    readNextFrame();

    set_trxframe_ePBC(impl_->fr, top.ePBC());
    if (top.hasTopology() && impl_->settings_.hasRmPBC())
    {
        impl_->gpbc_ = gmx_rmpbc_init(&top.topology()->idef, top.ePBC(),
                                      impl_->fr->natoms, impl_->fr->box);
    }
}

// TODO: Find a more efficient and correct way to count frames.
// Problems caused by this poor implementation:
// 1) Inefficient since it rewinds and reads all frames.
// 2) Function should be const.
// 3) Block scheduler reference to runner should be const.
// 4) Impl constructor in this file should also input a const runner.
// 5) Rewinds input file to front and manually sets nextFrameNumber.
// 6) Calls read_next_frame directly to bypass scheduler.

size_t
TrajectoryAnalysisRunnerCommon::getTotalNumberOfFrames()
{
    if (!hasTrajectory()) {
        return 0;
    }

    rewind_trj(impl_->status_);
    size_t numFramesRead=0;
    // TODO: Avoid reading into official buffers
    while (read_next_frame(impl_->oenv_, impl_->status_, impl_->fr))
    {
        numFramesRead++;
    }

    rewind_trj(impl_->status_);
    impl_->nextFrameNumber_ = 0;
    return numFramesRead;
}


size_t
TrajectoryAnalysisRunnerCommon::getCurrentFrameNumber()
{
    return impl_->nextFrameNumber_ - 1;
}


bool
TrajectoryAnalysisRunnerCommon::readNextFrame()
{
    bool bContinue = hasTrajectory() && !impl_->scheduler_.done();
    if (bContinue)
    {
        size_t targetFrameNumber = impl_->scheduler_.selectNextFrameNumber();
        if (impl_->nextFrameNumber_ > targetFrameNumber)
        {
            rewind_trj (impl_->status_);
            impl_->nextFrameNumber_ = 0;
        }
        for (; impl_->nextFrameNumber_ <= targetFrameNumber;
             ++impl_->nextFrameNumber_)
        {
            // TODO: Avoid reading into official buffers when skipping.
            bContinue = read_next_frame(impl_->oenv_, impl_->status_, impl_->fr);
            if (!bContinue) {
                break;
            }
        }
    }

    if (bContinue)
    {
        return true;
    }
    else
    {
        impl_->finishTrajectory();
        return false;
    }
}


void
TrajectoryAnalysisRunnerCommon::initFrame()
{
    if (impl_->gpbc_ != NULL)
    {
        gmx_rmpbc_trxfr(impl_->gpbc_, impl_->fr);
    }
}


TrajectoryAnalysisRunnerCommon::HelpFlags
TrajectoryAnalysisRunnerCommon::helpFlags() const
{
    HelpFlags flags = 0;

    if (!impl_->bQuiet_)
    {
        flags |= efHelpShowOptions;
        if (impl_->bHelp_)
        {
            flags |= efHelpShowDescriptions;
        }
        if (impl_->bShowHidden_)
        {
            flags |= efHelpShowHidden;
        }
    }
    return flags;
}

bool
TrajectoryAnalysisRunnerCommon::hasTrajectory() const
{
    return !impl_->trjfile_.empty();
}


const TopologyInformation &
TrajectoryAnalysisRunnerCommon::topologyInformation() const
{
    return impl_->topInfo_;
}


t_trxframe &
TrajectoryAnalysisRunnerCommon::frame() const
{
    GMX_RELEASE_ASSERT(impl_->fr != NULL, "Frame not available when accessed");
    return *impl_->fr;
}

} // namespace gmx
