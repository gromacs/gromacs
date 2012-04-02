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

#include "fatalerror/exceptions.h"
#include "trajectoryanalysis/analysismodule.h"
#include "trajectoryanalysis/analysissettings.h"
#include "trajectoryanalysis/cmdlinerunner.h"
#include "trajectoryanalysis/nbsearch.h"

#endif
