/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013, by the GROMACS development team, led by
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
/*! \defgroup module_trajectoryanalysis Framework for Trajectory Analysis (trajectoryanalysis)
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
 * \if msc
 * The sequence charts below provides an overview of how the trajectory
 * analysis modules typically interact with other components.
 * The first chart shows the overall flow in the case of parallel (threaded)
 * execution and the interaction with the \ref module_analysisdata module.
 * Although parallelization has not yet been implemented, it has influenced the
 * design.
 * \msc
 *     runner [ label="CommandLineRunner", URL="\ref gmx::TrajectoryAnalysisCommandLineRunner" ],
 *     module [ label="analysis\nmodule" ],
 *     thread1 [ label="analysis\nthread 1" ],
 *     thread2 [ label="analysis\nthread 2" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     module box thread2 [ label="single TrajectoryAnalysisModule object",
 *                          URL="\ref gmx::TrajectoryAnalysisModule" ];
 *     module => data [ label="create (in constructor)" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     runner => runner [ label="parse user input" ];
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     module => data [ label="initialize" ];
 *     runner << module;
 *     runner => runner [ label="read frame 0" ];
 *     runner => module [ label="initAfterFirstFrame()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAfterFirstFrame()" ];
 *     runner => thread1 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread1 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread1;
 *     runner => thread2 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread2 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread2;
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => thread1 [ label="analyzeFrame(0)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => thread2 [ label="analyzeFrame(1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="finishFrame(1)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner << thread2 [ label="analyzeFrame(1)" ];
 *     runner => runner [ label="read and initialize frame 2" ];
 *     runner => thread2 [ label="analyzeFrame(2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="finishFrame(0)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     data => data [ label="process frame 0" ];
 *     data => data [ label="process frame 1" ];
 *     runner << thread1 [ label="analyzeFrame(0)" ];
 *     ...;
 *     --- [ label="loop over frames ends" ];
 *     runner => thread1 [ label="finishFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread1 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread1 -> module [ label="accumulate results" ];
 *     runner << thread1;
 *     runner => thread2 [ label="finishFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread2 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread2 -> module [ label="accumulate results" ];
 *     runner << thread2;
 *     runner => module [ label="finishAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::finishAnalysis()" ];
 *     runner => module [ label="writeOutput()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::writeOutput()" ];
 * \endmsc
 *
 * The second chart below shows the interaction with selections and options
 * with focus on selection options.  See \ref module_options for a general
 * overview of how options work.
 * \msc
 *     runner [ label="CommandLineRunner", URL="\ref gmx::TrajectoryAnalysisCommandLineRunner" ],
 *     options [ label="Options", URL="\ref module_options" ],
 *     selection [ label="selections", URL="\ref module_selection" ],
 *     module [ label="analysis module", URL="\ref gmx::TrajectoryAnalysisModule" ];
 *
 *     runner box selection [ label="all these objects are owned by the framework" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     module => options [ label="addOption(SelectionOption)",
 *                         URL="\ref gmx::SelectionOption" ];
 *     ...;
 *     runner << module;
 *     runner => options [ label="parse user input" ];
 *     options => selection [ label="parse selections" ];
 *     options -> module [ label="initialize Selection variables",
 *                         URL="\ref gmx::Selection" ];
 *     runner << options;
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     module => selection [ label="adjust SelectionOptions",
 *                           URL="\ref gmx::SelectionOptionInfo" ];
 *     runner << module;
 *     runner => selection [ label="prompt missing selections" ];
 *     runner => selection [ label="compile selections" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     module -> selection [ label="access compiled selections",
 *                           URL="\ref gmx::Selection" ];
 *     runner << module;
 *     ...;
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="read and initialize frame 0" ];
 *     runner => selection [ label="evaluate selections for frame 0" ];
 *     runner => module [ label="analyzeFrame(0)",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     module -> selection [ label="access evaluated selections",
 *                           URL="\ref gmx::Selection" ];
 *     runner << module;
 *     ...;
 *     --- [ label="loop over frames ends" ];
 *     ...;
 *     runner => module [ label="finishAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::finishAnalysis()" ];
 *     module -> selection [ label="access compiled selections",
 *                           URL="\ref gmx::Selection" ];
 *     runner << module;
 * \endmsc
 * \endif
 *
 * In addition to the framework for defining analysis modules, this module also
 * provides gmx::TrajectoryAnalysisCommandLineRunner, which implements a
 * command-line program that runs a certain analysis module.
 *
 * Internally, the module also defines a set of trajectory analysis modules that
 * can currently be accessed only through gmx::registerTrajectoryAnalysisModules.
 *
 * For an example of how to implement an analysis tool using the framework, see
 * \ref template.cpp.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for trajectory analysis framework
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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
