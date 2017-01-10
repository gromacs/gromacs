/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <algorithm>
#include <string>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "analysissettings-impl.h"

namespace gmx
{

class TrajectoryAnalysisRunnerCommon::Impl : public ITopologyProvider
{
    public:
        explicit Impl(TrajectoryAnalysisSettings *settings);
        ~Impl();

        bool hasTrajectory() const { return !trjfile_.empty(); }

        void initTopology(bool required);
        void initFirstFrame();
        void initFrameIndexGroup();
        void finishTrajectory();

        // From ITopologyProvider
        virtual gmx_mtop_t *getTopology(bool required)
        {
            initTopology(required);
            return topInfo_.mtop_;
        }
        virtual int getAtomCount()
        {
            if (!topInfo_.hasTopology())
            {
                if (trajectoryGroup_.isValid())
                {
                    GMX_THROW(InconsistentInputError("-fgroup is only supported when -s is also specified"));
                }
                // Read the first frame if we don't know the maximum number of
                // atoms otherwise.
                initFirstFrame();
                return fr->natoms;
            }
            return -1;
        }

        TrajectoryAnalysisSettings &settings_;
        TopologyInformation         topInfo_;

        //! Name of the trajectory file (empty if not provided).
        std::string                 trjfile_;
        //! Name of the topology file (empty if no topology provided).
        std::string                 topfile_;
        Selection                   trajectoryGroup_;
        double                      startTime_;
        double                      endTime_;
        double                      deltaTime_;
        bool                        bStartTimeSet_;
        bool                        bEndTimeSet_;
        bool                        bDeltaTimeSet_;

        bool                        bTrajOpen_;
        //! The current frame, or \p NULL if no frame loaded yet.
        t_trxframe                 *fr;
        gmx_rmpbc_t                 gpbc_;
        //! Used to store the status variable from read_first_frame().
        t_trxstatus                *status_;
        gmx_output_env_t           *oenv_;
};


TrajectoryAnalysisRunnerCommon::Impl::Impl(TrajectoryAnalysisSettings *settings)
    : settings_(*settings),
      startTime_(0.0), endTime_(0.0), deltaTime_(0.0),
      bStartTimeSet_(false), bEndTimeSet_(false), bDeltaTimeSet_(false),
      bTrajOpen_(false), fr(nullptr), gpbc_(nullptr), status_(nullptr), oenv_(nullptr)
{
}


TrajectoryAnalysisRunnerCommon::Impl::~Impl()
{
    finishTrajectory();
    if (fr != nullptr)
    {
        // There doesn't seem to be a function for freeing frame data
        sfree(fr->x);
        sfree(fr->v);
        sfree(fr->f);
        sfree(fr->index);
        sfree(fr);
    }
    if (oenv_ != nullptr)
    {
        output_env_done(oenv_);
    }
}

void
TrajectoryAnalysisRunnerCommon::Impl::initTopology(bool required)
{
    // Return immediately if the topology has already been loaded.
    if (topInfo_.hasTopology())
    {
        return;
    }

    if (required && topfile_.empty())
    {
        GMX_THROW(InconsistentInputError("No topology provided, but one is required for analysis"));
    }

    // Load the topology if requested.
    if (!topfile_.empty())
    {
        snew(topInfo_.mtop_, 1);
        readConfAndTopology(topfile_.c_str(), &topInfo_.bTop_, topInfo_.mtop_,
                            &topInfo_.ePBC_, &topInfo_.xtop_, nullptr,
                            topInfo_.boxtop_);
        // TODO: Only load this here if the tool actually needs it; selections
        // take care of themselves.
        for (int i = 0; i < topInfo_.mtop_->nmoltype; ++i)
        {
            gmx_moltype_t &moltype = topInfo_.mtop_->moltype[i];
            if (!moltype.atoms.haveMass)
            {
                // Try to read masses from database, be silent about missing masses
                atomsSetMassesBasedOnNames(&moltype.atoms, FALSE);
            }
        }
        if (hasTrajectory()
            && !settings_.hasFlag(TrajectoryAnalysisSettings::efUseTopX))
        {
            sfree(topInfo_.xtop_);
            topInfo_.xtop_ = nullptr;
        }
    }
}

void
TrajectoryAnalysisRunnerCommon::Impl::initFirstFrame()
{
    // Return if we have already initialized the trajectory.
    if (fr != nullptr)
    {
        return;
    }
    time_unit_t time_unit
        = static_cast<time_unit_t>(settings_.timeUnit() + 1);
    output_env_init(&oenv_, getProgramContext(), time_unit, FALSE, exvgNONE, 0);

    int frflags = settings_.frflags();
    frflags |= TRX_NEED_X;

    snew(fr, 1);

    if (hasTrajectory())
    {
        if (!read_first_frame(oenv_, &status_, trjfile_.c_str(), fr, frflags))
        {
            GMX_THROW(FileIOError("Could not read coordinates from trajectory"));
        }
        bTrajOpen_ = true;

        if (topInfo_.hasTopology())
        {
            const int topologyAtomCount = topInfo_.topology()->atoms.nr;
            if (fr->natoms > topologyAtomCount)
            {
                const std::string message
                    = formatString("Trajectory (%d atoms) does not match topology (%d atoms)",
                                   fr->natoms, topologyAtomCount);
                GMX_THROW(InconsistentInputError(message));
            }
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
        fr->natoms = topInfo_.topology()->atoms.nr;
        fr->bX     = TRUE;
        snew(fr->x, fr->natoms);
        memcpy(fr->x, topInfo_.xtop_,
               sizeof(*fr->x) * fr->natoms);
        fr->bBox   = TRUE;
        copy_mat(topInfo_.boxtop_, fr->box);
    }

    set_trxframe_ePBC(fr, topInfo_.ePBC());
    if (topInfo_.hasTopology() && settings_.hasRmPBC())
    {
        gpbc_ = gmx_rmpbc_init(&topInfo_.topology()->idef, topInfo_.ePBC(),
                               fr->natoms);
    }
}

void
TrajectoryAnalysisRunnerCommon::Impl::initFrameIndexGroup()
{
    if (!trajectoryGroup_.isValid())
    {
        return;
    }
    GMX_RELEASE_ASSERT(bTrajOpen_,
                       "Trajectory index only makes sense with a real trajectory");
    if (trajectoryGroup_.atomCount() != fr->natoms)
    {
        const std::string message = formatString(
                    "Selection specified with -fgroup has %d atoms, but "
                    "the trajectory (-f) has %d atoms.",
                    trajectoryGroup_.atomCount(), fr->natoms);
        GMX_THROW(InconsistentInputError(message));
    }
    fr->bIndex = TRUE;
    snew(fr->index, trajectoryGroup_.atomCount());
    std::copy(trajectoryGroup_.atomIndices().begin(),
              trajectoryGroup_.atomIndices().end(),
              fr->index);
}

void
TrajectoryAnalysisRunnerCommon::Impl::finishTrajectory()
{
    if (bTrajOpen_)
    {
        close_trx(status_);
        bTrajOpen_ = false;
    }
    if (gpbc_ != nullptr)
    {
        gmx_rmpbc_done(gpbc_);
        gpbc_ = nullptr;
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


ITopologyProvider *
TrajectoryAnalysisRunnerCommon::topologyProvider()
{
    return impl_.get();
}


void
TrajectoryAnalysisRunnerCommon::initOptions(IOptionsContainer *options,
                                            TimeUnitBehavior  *timeUnitBehavior)
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

    // Add options for trajectory time control.
    options->addOption(DoubleOption("b")
                           .store(&impl_->startTime_).storeIsSet(&impl_->bStartTimeSet_)
                           .timeValue()
                           .description("First frame (%t) to read from trajectory"));
    options->addOption(DoubleOption("e")
                           .store(&impl_->endTime_).storeIsSet(&impl_->bEndTimeSet_)
                           .timeValue()
                           .description("Last frame (%t) to read from trajectory"));
    options->addOption(DoubleOption("dt")
                           .store(&impl_->deltaTime_).storeIsSet(&impl_->bDeltaTimeSet_)
                           .timeValue()
                           .description("Only use frame if t MOD dt == first time (%t)"));

    // Add time unit option.
    timeUnitBehavior->setTimeUnitFromEnvironment();
    timeUnitBehavior->addTimeUnitOption(options, "tu");
    timeUnitBehavior->setTimeUnitStore(&impl_->settings_.impl_->timeUnit);

    options->addOption(SelectionOption("fgroup")
                           .store(&impl_->trajectoryGroup_)
                           .onlySortedAtoms().onlyStatic()
                           .description("Atoms stored in the trajectory file "
                                        "(if not set, assume first N atoms)"));

    // Add plot options.
    settings.impl_->plotSettings.initOptions(options);

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
TrajectoryAnalysisRunnerCommon::optionsFinished()
{
    if (impl_->trjfile_.empty() && impl_->topfile_.empty())
    {
        GMX_THROW(InconsistentInputError("No trajectory or topology provided, nothing to do!"));
    }

    if (impl_->trajectoryGroup_.isValid() && impl_->trjfile_.empty())
    {
        GMX_THROW(InconsistentInputError("-fgroup only makes sense together with a trajectory (-f)"));
    }

    impl_->settings_.impl_->plotSettings.setTimeUnit(impl_->settings_.timeUnit());

    if (impl_->bStartTimeSet_)
    {
        setTimeValue(TBEGIN, impl_->startTime_);
    }
    if (impl_->bEndTimeSet_)
    {
        setTimeValue(TEND, impl_->endTime_);
    }
    if (impl_->bDeltaTimeSet_)
    {
        setTimeValue(TDELTA, impl_->deltaTime_);
    }
}


void
TrajectoryAnalysisRunnerCommon::initTopology()
{
    const bool topologyRequired =
        impl_->settings_.hasFlag(TrajectoryAnalysisSettings::efRequireTop);
    impl_->initTopology(topologyRequired);
}


void
TrajectoryAnalysisRunnerCommon::initFirstFrame()
{
    impl_->initFirstFrame();
}


void
TrajectoryAnalysisRunnerCommon::initFrameIndexGroup()
{
    impl_->initFrameIndexGroup();
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
    if (impl_->gpbc_ != nullptr)
    {
        gmx_rmpbc_trxfr(impl_->gpbc_, impl_->fr);
    }
}


bool
TrajectoryAnalysisRunnerCommon::hasTrajectory() const
{
    return impl_->hasTrajectory();
}


const TopologyInformation &
TrajectoryAnalysisRunnerCommon::topologyInformation() const
{
    return impl_->topInfo_;
}


t_trxframe &
TrajectoryAnalysisRunnerCommon::frame() const
{
    GMX_RELEASE_ASSERT(impl_->fr != nullptr, "Frame not available when accessed");
    return *impl_->fr;
}

} // namespace gmx
