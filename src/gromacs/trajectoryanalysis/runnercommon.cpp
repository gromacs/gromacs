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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"
#include "gromacs/utility/gmxassert.h"

#include "analysissettings-impl.h"

namespace gmx
{

class TrajectoryAnalysisRunnerCommon::Impl
{
    public:
        Impl(TrajectoryAnalysisSettings *settings);
        ~Impl();

        void finishTrajectory();

        TrajectoryAnalysisSettings &_settings;
        Options                 _options;
        TopologyInformation     _topInfo;

        bool                    _bHelp;
        bool                    _bShowHidden;
        bool                    _bQuiet;
        //! Name of the trajectory file (empty if not provided).
        std::string             _trjfile;
        //! Name of the topology file (empty if no topology provided).
        std::string             _topfile;
        //! Name of the index file (empty if no index file provided).
        std::string             _ndxfile;
        double                  _startTime;
        double                  _endTime;
        double                  _deltaTime;

        gmx_ana_indexgrps_t    *_grps;
        bool                    _bTrajOpen;
        //! The current frame, or \p NULL if no frame loaded yet.
        t_trxframe          *fr;
        gmx_rmpbc_t             _gpbc;
        //! Used to store the status variable from read_first_frame().
        t_trxstatus            *_status;
        output_env_t            _oenv;
};


TrajectoryAnalysisRunnerCommon::Impl::Impl(TrajectoryAnalysisSettings *settings)
    : _settings(*settings), _options("common", "Common analysis control"),
      _bHelp(false), _bShowHidden(false), _bQuiet(false),
      _startTime(0.0), _endTime(0.0), _deltaTime(0.0),
      _grps(NULL),
      _bTrajOpen(false), fr(NULL), _gpbc(NULL), _status(NULL), _oenv(NULL)
{
}


TrajectoryAnalysisRunnerCommon::Impl::~Impl()
{
    if (_grps != NULL)
    {
        gmx_ana_indexgrps_free(_grps);
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
    if (_oenv != NULL)
    {
        output_env_done(_oenv);
    }
}


void
TrajectoryAnalysisRunnerCommon::Impl::finishTrajectory()
{
    if (_bTrajOpen)
    {
        close_trx(_status);
        _bTrajOpen = false;
    }
    if (_gpbc != NULL)
    {
        gmx_rmpbc_done(_gpbc);
        _gpbc = NULL;
    }
}

/*********************************************************************
 * TrajectoryAnalysisRunnerCommon
 */

TrajectoryAnalysisRunnerCommon::TrajectoryAnalysisRunnerCommon(
        TrajectoryAnalysisSettings *settings)
    : _impl(new Impl(settings))
{
}


TrajectoryAnalysisRunnerCommon::~TrajectoryAnalysisRunnerCommon()
{
}


Options &
TrajectoryAnalysisRunnerCommon::initOptions()
{
    TrajectoryAnalysisSettings &settings = _impl->_settings;
    Options &options = _impl->_options;

    // Add options for help.
    options.addOption(BooleanOption("h").store(&_impl->_bHelp)
                          .description("Print help and quit"));
    options.addOption(BooleanOption("hidden").store(&_impl->_bShowHidden)
                          .hidden()
                          .description("Show hidden options"));
    options.addOption(BooleanOption("quiet").store(&_impl->_bQuiet)
                          .hidden()
                          .description("Hide options in normal run"));

    // Add common file name arguments.
    options.addOption(FileNameOption("f")
                          .filetype(eftTrajectory).inputFile()
                          .store(&_impl->_trjfile)
                          .description("Input trajectory"));
    options.addOption(FileNameOption("s")
                          .filetype(eftTopology).inputFile()
                          .store(&_impl->_topfile)
                          .description("Input topology"));
    options.addOption(FileNameOption("n")
                          .filetype(eftIndex).inputFile()
                          .store(&_impl->_ndxfile)
                          .description("Extra index groups"));
    options.addOption(SelectionFileOption("sf"));

    // Add options for trajectory time control.
    options.addOption(DoubleOption("b").store(&_impl->_startTime).timeValue()
                          .description("First frame (%t) to read from trajectory"));
    options.addOption(DoubleOption("e").store(&_impl->_endTime).timeValue()
                          .description("Last frame (%t) to read from trajectory"));
    options.addOption(DoubleOption("dt").store(&_impl->_deltaTime).timeValue()
                          .description("Only use frame if t MOD dt == first time (%t)"));

    // Add time unit option.
    settings._impl->timeUnitManager.addTimeUnitOption(&options, "tu");

    // Add plot options.
    settings._impl->plotSettings.addOptions(&options);

    // Add common options for trajectory processing.
    if (!settings.hasFlag(TrajectoryAnalysisSettings::efNoUserRmPBC))
    {
        options.addOption(BooleanOption("rmpbc").store(&settings._impl->bRmPBC)
                              .description("Make molecules whole for each frame"));
    }
    if (!settings.hasFlag(TrajectoryAnalysisSettings::efNoUserPBC))
    {
        options.addOption(BooleanOption("pbc").store(&settings._impl->bPBC)
                              .description("Use periodic boundary conditions for distance calculation"));
    }

    return _impl->_options;
}


void
TrajectoryAnalysisRunnerCommon::scaleTimeOptions(Options *options)
{
    _impl->_settings._impl->timeUnitManager.scaleTimeOptions(options);
}


bool
TrajectoryAnalysisRunnerCommon::initOptionsDone()
{
    if (_impl->_bHelp)
    {
        return false;
    }

    _impl->_settings._impl->plotSettings.setTimeUnit(
            _impl->_settings._impl->timeUnitManager.timeUnit());

    if (_impl->_trjfile.empty() && _impl->_topfile.empty())
    {
        GMX_THROW(InconsistentInputError("No trajectory or topology provided, nothing to do!"));
    }

    if (_impl->_options.isSet("b"))
        setTimeValue(TBEGIN, _impl->_startTime);
    if (_impl->_options.isSet("e"))
        setTimeValue(TEND, _impl->_endTime);
    if (_impl->_options.isSet("dt"))
        setTimeValue(TDELTA, _impl->_deltaTime);

    return true;
}


void
TrajectoryAnalysisRunnerCommon::initIndexGroups(SelectionCollection *selections)
{
    if (_impl->_ndxfile.empty())
    {
        // TODO: Initialize default selections
        selections->setIndexGroups(NULL);
    }
    else
    {
        gmx_ana_indexgrps_init(&_impl->_grps, NULL, _impl->_ndxfile.c_str());
        selections->setIndexGroups(_impl->_grps);
    }
}


void
TrajectoryAnalysisRunnerCommon::doneIndexGroups(SelectionCollection *selections)
{
    if (_impl->_grps != NULL)
    {
        selections->setIndexGroups(NULL);
        gmx_ana_indexgrps_free(_impl->_grps);
        _impl->_grps = NULL;
    }
}


void
TrajectoryAnalysisRunnerCommon::initTopology(SelectionCollection *selections)
{
    const TrajectoryAnalysisSettings &settings = _impl->_settings;
    bool bRequireTop
        = settings.hasFlag(TrajectoryAnalysisSettings::efRequireTop)
          || selections->requiresTopology();
    if (bRequireTop && _impl->_topfile.empty())
    {
        GMX_THROW(InconsistentInputError("No topology provided, but one is required for analysis"));
    }

    // Load the topology if requested.
    if (!_impl->_topfile.empty())
    {
        char  title[STRLEN];

        snew(_impl->_topInfo._top, 1);
        _impl->_topInfo._bTop = read_tps_conf(_impl->_topfile.c_str(), title,
                _impl->_topInfo._top, &_impl->_topInfo._ePBC,
                &_impl->_topInfo._xtop, NULL, _impl->_topInfo._boxtop, TRUE);
        if (hasTrajectory()
            && !settings.hasFlag(TrajectoryAnalysisSettings::efUseTopX))
        {
            sfree(_impl->_topInfo._xtop);
            _impl->_topInfo._xtop = NULL;
        }
    }

    // Read the first frame if we don't know the maximum number of atoms
    // otherwise.
    int  natoms = -1;
    if (!_impl->_topInfo.hasTopology())
    {
        initFirstFrame();
        natoms = _impl->fr->natoms;
    }
    selections->setTopology(_impl->_topInfo.topology(), natoms);

    /*
    if (_impl->bSelDump)
    {
        gmx_ana_poscalc_coll_print_tree(stderr, _impl->pcc);
        fprintf(stderr, "\n");
    }
    */
}


void
TrajectoryAnalysisRunnerCommon::initFirstFrame()
{
    // Return if we have already initialized the trajectory.
    if (_impl->fr)
    {
        return;
    }
    snew(_impl->_oenv, 1);
    output_env_init_default(_impl->_oenv);
    _impl->_oenv->time_unit
        = static_cast<time_unit_t>(_impl->_settings.timeUnit() + 1);

    int frflags = _impl->_settings.frflags();
    frflags |= TRX_NEED_X;

    snew(_impl->fr, 1);

    const TopologyInformation &top = _impl->_topInfo;
    if (hasTrajectory())
    {
        if (!read_first_frame(_impl->_oenv, &_impl->_status,
                              _impl->_trjfile.c_str(), _impl->fr, frflags))
        {
            GMX_THROW(FileIOError("Could not read coordinates from trajectory"));
        }
        _impl->_bTrajOpen = true;

        if (top.hasTopology() && _impl->fr->natoms > top.topology()->atoms.nr)
        {
            GMX_THROW(InconsistentInputError(formatString(
                      "Trajectory (%d atoms) does not match topology (%d atoms)",
                      _impl->fr->natoms, top.topology()->atoms.nr)));
        }
        // Check index groups if they have been initialized based on the topology.
        /*
        if (top)
        {
            for (int i = 0; i < _impl->sel->nr(); ++i)
            {
                gmx_ana_index_check(_impl->sel->sel(i)->indexGroup(),
                                    _impl->fr->natoms);
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
        _impl->fr->flags  = frflags;
        _impl->fr->natoms = top.topology()->atoms.nr;
        _impl->fr->bX     = TRUE;
        snew(_impl->fr->x, _impl->fr->natoms);
        memcpy(_impl->fr->x, top._xtop,
               sizeof(*_impl->fr->x) * _impl->fr->natoms);
        _impl->fr->bBox   = TRUE;
        copy_mat(const_cast<rvec *>(top._boxtop), _impl->fr->box);
    }

    set_trxframe_ePBC(_impl->fr, top.ePBC());
    if (top.hasTopology() && _impl->_settings.hasRmPBC())
    {
        _impl->_gpbc = gmx_rmpbc_init(&top.topology()->idef, top.ePBC(),
                                      _impl->fr->natoms, _impl->fr->box);
    }
}


bool
TrajectoryAnalysisRunnerCommon::readNextFrame()
{
    bool bContinue = false;
    if (hasTrajectory())
    {
        bContinue = read_next_frame(_impl->_oenv, _impl->_status, _impl->fr);
    }
    if (!bContinue)
    {
        _impl->finishTrajectory();
    }
    return bContinue;
}


void
TrajectoryAnalysisRunnerCommon::initFrame()
{
    if (_impl->_gpbc != NULL)
    {
        gmx_rmpbc_trxfr(_impl->_gpbc, _impl->fr);
    }
}


TrajectoryAnalysisRunnerCommon::HelpFlags
TrajectoryAnalysisRunnerCommon::helpFlags() const
{
    HelpFlags flags = 0;

    if (!_impl->_bQuiet)
    {
        flags |= efHelpShowOptions;
        if (_impl->_bHelp)
        {
            flags |= efHelpShowDescriptions;
        }
        if (_impl->_bShowHidden)
        {
            flags |= efHelpShowHidden;
        }
    }
    return flags;
}

bool
TrajectoryAnalysisRunnerCommon::hasTrajectory() const
{
    return !_impl->_trjfile.empty();
}


const TopologyInformation &
TrajectoryAnalysisRunnerCommon::topologyInformation() const
{
    return _impl->_topInfo;
}


t_trxframe &
TrajectoryAnalysisRunnerCommon::frame() const
{
    GMX_RELEASE_ASSERT(_impl->fr != NULL, "Frame not available when accessed");
    return *_impl->fr;
}

} // namespace gmx
