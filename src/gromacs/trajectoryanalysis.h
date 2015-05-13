/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
 * The sequence charts below provides an overview of how the trajectory
 * analysis modules typically interact with other components.
 * The first chart provides an overview of the call sequence of the most
 * important methods in gmx::TrajectoryAnalysisModule.
 * There is a runner, which is responsible for doing the work that is shared
 * between all trajectory analysis (such as reading the trajectory and
 * processing selections).  The runner then calls different methods in the
 * analysis module at appropriate points to perform the module-specific tasks.
 * The analysis module is responsible for creating and managing
 * gmx::AnalysisData objects, and the chart shows the most important
 * interactions with this module as well.  However, the runner takes
 * responsibility of calling gmx::AnalysisData::finishFrameSerial().
 * Interactions with options (for command-line option processing) and
 * selections is not shown for brevity: see \ref module_options for an overview
 * of how options work, and the second chart for a more detailed view of how
 * selections are accessed from an analysis module.
 * \msc
 *     runner,
 *     module [ URL="\ref gmx::TrajectoryAnalysisModule" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     runner box module [ label="caller owns runner and module objects" ];
 *     module => data [ label="create (in constructor)" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     runner => runner [ label="parse user input" ];
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     runner => runner [ label="initialize topology\nand selections" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     module => data [ label="initialize" ];
 *     runner => runner [ label="read frame 0" ];
 *     runner => module [ label="initAfterFirstFrame()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAfterFirstFrame()" ];
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => module [ label="analyzeFrame(0)",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     module => data [ label="add data",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     module => data [ label="finishFrame()",
 *                      URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner => data [ label="finishFrameSerial()",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => module [ label="analyzeFrame(1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     ...;
 *     --- [ label="loop over frames ends" ];
 *     runner => module [ label="finishAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::finishAnalysis()" ];
 *     module => data [ label="post-process data" ];
 *     runner => module [ label="writeOutput()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::writeOutput()" ];
 * \endmsc
 *
 * The second chart below shows the interaction with selections and options
 * with focus on selection options.  The gmx::TrajectoryAnalysisModule object
 * creates one or more gmx::Selection variables, and uses gmx::SelectionOption
 * to indicate them as the destination for selections.  This happens in
 * gmx::TrajectoryAnalysisModule::initOptions().  After the options have been
 * parsed (includes parsing any options present on the command-line or read
 * from files, but not those provided interactively),
 * gmx::TrajectoryAnalysisModule::optionsFinished() can adjust the selections
 * using gmx::SelectionOptionInfo.  This is done like this to allow the
 * analysis module to influence the interactive prompt of selections based on
 * what command-line options were given.  After optionsFinished() returns, the
 * interactive selection prompt is presented if necessary.  After this point,
 * all access to selections from the analysis module is through the
 * gmx::Selection variables: the runner is responsible for calling methods in
 * the selection library, and these methods update the content referenced by
 * the gmx::Selection variables.  See documentation of
 * gmx::TrajectoryAnalysisModule for details of what the selections contain at
 * each point.
 * \msc
 *     runner,
 *     options [ label="Options", URL="\ref module_options" ],
 *     selection [ label="selections", URL="\ref module_selection" ],
 *     module [ label="module", URL="\ref gmx::TrajectoryAnalysisModule" ];
 *
 *     runner box selection [ label="all these objects are owned by the framework" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     module => options [ label="addOption(SelectionOption)",
 *                         URL="\ref gmx::SelectionOption" ];
 *     module => options [ label="addOption() (other options)",
 *                         URL="\ref gmx::Options::addOption()" ];
 *     ...;
 *     runner << module;
 *     runner => options [ label="parse command-line parameters" ];
 *     options => selection [ label="parse selections" ];
 *     selection -> module [ label="initialize Selection variables",
 *                           URL="\ref gmx::Selection" ];
 *     runner << options;
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     module => selection [ label="adjust SelectionOptions",
 *                           URL="\ref gmx::SelectionOptionInfo" ];
 *     runner << module;
 *     runner => selection [ label="prompt missing selections" ];
 *     selection -> module [ label="initialize Selection variables",
 *                         URL="\ref gmx::Selection" ];
 *     runner => selection [ label="compile selections" ];
 *     selection -> module [ label="change content referenced\nby Selection variables" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     ...;
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="read and initialize frame 0" ];
 *     runner => selection [ label="evaluate selections for frame 0" ];
 *     selection -> module [ label="change content referenced\nby Selection variables" ];
 *     ...;
 * \endmsc
 *
 * The final chart shows the flow within the frame loop in the case of parallel
 * (threaded) execution and the interaction with the \ref module_analysisdata
 * module in this case.  Although parallelization has not yet been implemented,
 * it has influenced the design and needs to be understood if one wants to
 * write modules that can take advantage of the parallelization once it gets
 * implemented.  The parallelization takes part over frames: analyzing a single
 * frame is one unit of work.  When the frame loop is started,
 * gmx::TrajectoryAnalysisModule::startFrames() is called for each thread, and
 * initializes an object that contains thread-local data needed during the
 * analysis.  This includes selection information, gmx::AnalysisDataHandle
 * objects, and possibly other module-specific variables.  Then, the runner
 * reads the frames in sequence and passes the work into the different threads,
 * together with the appropriate thread-local data object.
 * The gmx::TrajectoryAnalysisModule::analyzeFrame() calls are only allowed to modify
 * the thread-local data object; everything else is read-only.  For any output,
 * they pass the information to gmx::AnalysisData, which together with the
 * runner takes care of ordering the data from different frames such that it
 * gets processed in the right order.
 * When all frames are analyzed, gmx::TrajectoryAnalysisModule::finishFrames()
 * is called for each thread-local data object to destroy them and to
 * accumulate possible results from them into the main
 * gmx::TrajectoryAnalysisModule object.
 * Note that in the diagram, some part of the work attributed for the runner
 * (e.g., evaluating selections) will actually be carried out in the analysis
 * threads before gmx::TrajectoryAnalysisModule::analyzeFrame() gets called.
 * \msc
 *     runner,
 *     module [ label="module object" ],
 *     thread1 [ label="analysis\nthread 1" ],
 *     thread2 [ label="analysis\nthread 2" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     module box thread2 [ label="single TrajectoryAnalysisModule object",
 *                          URL="\ref gmx::TrajectoryAnalysisModule" ];
 *     ...;
 *     --- [ label="loop over frames starts" ];
 *     runner => thread1 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread1 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread1 [ label="pdata1" ];
 *     runner => thread2 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread2 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread2 [ label="pdata2" ];
 *     |||;
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => thread1 [ label="analyzeFrame(0, pdata1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => thread2 [ label="analyzeFrame(1, pdata2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="finishFrame(1)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner << thread2 [ label="analyzeFrame() (frame 1)" ];
 *     runner => runner [ label="read and initialize frame 2" ];
 *     runner => thread2 [ label="analyzeFrame(2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="finishFrame(0)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner << thread1 [ label="analyzeFrame() (frame 0)" ];
 *     runner => data [ label="finishFrameSerial() (frame 0)",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     runner => data [ label="finishFrameSerial() (frame 1)",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     ...;
 *     runner => thread1 [ label="finishFrames(pdata1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread1 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread1 -> module [ label="accumulate results" ];
 *     runner << thread1;
 *     runner => thread2 [ label="finishFrames(pdata2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread2 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread2 -> module [ label="accumulate results" ];
 *     runner << thread2;
 *     --- [ label="loop over frames ends" ];
 *     ...;
 * \endmsc
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

#include "gromacs/analysisdata.h"
#include "gromacs/options.h"
#include "gromacs/selection.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#endif
