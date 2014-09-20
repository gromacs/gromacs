/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::TrajectoryAnalysisRunnerCommon.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "runnercommon.h"

#include <string.h>

#include "gromacs/fileio/timecontrol.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "analysissettings-impl.h"

namespace gmx
{

class TrajectoryAnalysisRunnerCommon::Impl
{
    public:
        Impl(TrajectoryAnalysisSettings *settings);
        ~Impl();

        void finishTrajectory();

        TrajectoryAnalysisSettings &settings_;
        TopologyInformation         topInfo_;

        //! Name of the trajectory file (empty if not provided).
        std::string                 trjfile_;
        //! Name of the topology file (empty if no topology provided).
        std::string                 topfile_;
        //! Name of the index file (empty if no index file provided).
        std::string                 ndxfile_;
        double                      startTime_;
        double                      endTime_;
        double                      deltaTime_;

        gmx_ana_indexgrps_t        *grps_;
        bool                        bTrajOpen_;
        //! The current frame, or \p NULL if no frame loaded yet.
        t_trxframe                 *fr;
        gmx_rmpbc_t                 gpbc_;
        //! Used to store the status variable from read_first_frame().
        t_trxstatus                *status_;
        output_env_t                oenv_;
};


TrajectoryAnalysisRunnerCommon::Impl::Impl(TrajectoryAnalysisSettings *settings)
    : settings_(*settings),
      startTime_(0.0), endTime_(0.0), deltaTime_(0.0),
      grps_(NULL),
      bTrajOpen_(false), fr(NULL), gpbc_(NULL), status_(NULL), oenv_(NULL)
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
    : impl_(new Impl(settings))
{
}


TrajectoryAnalysisRunnerCommon::~TrajectoryAnalysisRunnerCommon()
{
}


void
TrajectoryAnalysisRunnerCommon::initOptions(Options *options)
{
    TrajectoryAnalysisSettings &settings = impl_->settings_;

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

    options->addOption(SelectionFileOption("sf"));
}


void
TrajectoryAnalysisRunnerCommon::scaleTimeOptions(Options *options)
{
    impl_->settings_.impl_->timeUnitManager.scaleTimeOptions(options);
}


void
TrajectoryAnalysisRunnerCommon::optionsFinished(Options *options)
{
    impl_->settings_.impl_->plotSettings.setTimeUnit(
            impl_->settings_.impl_->timeUnitManager.timeUnit());

    if (impl_->trjfile_.empty() && impl_->topfile_.empty())
    {
        GMX_THROW(InconsistentInputError("No trajectory or topology provided, nothing to do!"));
    }

    if (options->isSet("b"))
    {
        setTimeValue(TBEGIN, impl_->startTime_);
    }
    if (options->isSet("e"))
    {
        setTimeValue(TEND, impl_->endTime_);
    }
    if (options->isSet("dt"))
    {
        setTimeValue(TDELTA, impl_->deltaTime_);
    }
}


void
TrajectoryAnalysisRunnerCommon::initIndexGroups(SelectionCollection *selections,
                                                bool                 bUseDefaults)
{
    if (impl_->ndxfile_.empty())
    {
        if (!bUseDefaults)
        {
            selections->setIndexGroups(NULL);
            return;
        }
        initTopology(selections);
    }
    const char *const ndxfile
        = (!impl_->ndxfile_.empty() ? impl_->ndxfile_.c_str() : NULL);
    gmx_ana_indexgrps_init(&impl_->grps_, impl_->topInfo_.topology(), ndxfile);
    selections->setIndexGroups(impl_->grps_);
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
    // Return immediately if the topology has already been loaded.
    if (impl_->topInfo_.hasTopology())
    {
        return;
    }

    const TrajectoryAnalysisSettings &settings = impl_->settings_;
    const bool bRequireTop
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
    output_env_init(&impl_->oenv_, getProgramContext(), time_unit, FALSE, exvgNONE, 0);

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

    set_trxframe_ePBC(impl_->fr, top.ePBC());
    if (top.hasTopology() && impl_->settings_.hasRmPBC())
    {
        impl_->gpbc_ = gmx_rmpbc_init(&top.topology()->idef, top.ePBC(),
                                      impl_->fr->natoms);
    }
}


bool
TrajectoryAnalysisRunnerCommon::readNextFrame()
{
    bool bContinue = false;
    if (hasTrajectory())
    {
        bContinue = read_next_frame(impl_->oenv_, impl_->status_, impl_->fr);
    }
    if (!bContinue)
    {
        impl_->finishTrajectory();
    }
    return bContinue;
}


void
TrajectoryAnalysisRunnerCommon::initFrame()
{
    if (impl_->gpbc_ != NULL)
    {
        gmx_rmpbc_trxfr(impl_->gpbc_, impl_->fr);
    }
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
