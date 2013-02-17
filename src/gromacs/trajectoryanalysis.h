/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \defgroup module_trajectoryanalysis Framework for Trajectory Analysis
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for implementing trajectory analysis modules.
 *
 * This module implements a framework for implementing flexible trajectory
 * analysis routines.  It provides a base class for implementing analysis as
 * reusable modules that can be used from different contexts and can also
 * support per-frame parallelization.  It integrally uses functionality from the
 * following modules:
 *  - \ref module_options
 *  - \ref module_analysisdata
 *  - \ref module_selection
 *
 * The main interface of this module is the gmx::TrajectoryAnalysisModule class.
 * Analysis modules should derive from this class, and override the necessary
 * virtual methods to provide the actual initialization and analysis routines.
 * Classes gmx::TrajectoryAnalysisSettings and gmx::TopologyInformation (in
 * addition to classes declared in the above-mentioned modules) are used to pass
 * information to and from these methods.  gmx::TrajectoryAnalysisModuleData can
 * be used in advanced scenarios where the tool requires custom thread-local
 * data for parallel analysis.
 *
 * In addition to the framework for defining analysis modules, this module also
 * provides gmx::TrajectoryAnalysisCommandLineRunner, which implements a
 * command-line program that runs a certain analysis module.
 *
 * Internally, the module also defines a set of trajectory analysis modules that
 * can be instantiated using createTrajectoryAnalysisModule().
 *
 * For an example of how to implement an analysis tool using the framework, see
 * \ref template.cpp.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for trajectory analysis framework
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_H
#define GMX_TRAJECTORYANALYSIS_H

#include "analysisdata.h"
#include "options.h"
#include "selection.h"

#include "selection/nbsearch.h"
#include "trajectoryanalysis/analysismodule.h"
#include "trajectoryanalysis/analysissettings.h"
#include "trajectoryanalysis/cmdlinerunner.h"
#include "utility/exceptions.h"
#include "utility/programinfo.h"
#include "utility/stringutil.h"

#endif
