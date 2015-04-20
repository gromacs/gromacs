/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_analysisdata Parallelizable Handling of Output Data (analysisdata)
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for handling and processing output data from
 * analysis.
 *
 * <H3>Overview</H3>
 *
 * This module provides functionality to do common processing for tabular data
 * in analysis tools.  In addition to providing this common functionality, one
 * major driver for this module is to make it simple to write analysis tools
 * that process frames in parallel: the functionality in this module takes care
 * of necessary synchronization and communication such that output from the
 * frames is collected and output in the correct order.
 * See \ref page_analysisdata for an overview of the high-level functionality
 * and the terminology used.
 *
 * This module consists of two main parts.  The first is formed by the
 * gmx::AbstractAnalysisData class and classes that derive from it:
 * gmx::AnalysisData and gmx::AnalysisArrayData.  These classes are used to
 * process and store raw data as produced by the analysis tool.  They also
 * provide an interface to attach data modules that implement
 * gmx::AnalysisDataModuleInterface.
 *
 * Modules that implement gmx::AnalysisDataModuleInterface form the second part
 * of the module, and they provide functionality to do processing on the data.
 * These modules can also derive from gmx::AbstractAnalysisData, allowing other
 * modules to be attached to them to form a processing chain that best suits
 * the analysis tool.  Typically, such a processing chain ends in a plotting
 * module that writes the data into a file, but the final module can also
 * provide direct access to the processed data, allowing the analysis tool to
 * do custom postprocessing outside the module framework.
 *
 * <H3>Using Data Objects and Modules</H3>
 *
 * To use the functionality in this module, you typically declare one or more
 * AnalysisData objects and set its properties.  You then create some module
 * objects and set their properties (see the list of classes that implement
 * gmx::AnalysisDataModuleInterface) and attach them to the data objects or to
 * one another using gmx::AbstractAnalysisData::addModule().  Then you add the
 * actual data values to the gmx::AnalysisData object, which automatically
 * passes it on to the modules.
 * After all data is added, you may optionally access some results directly
 * from the module objects or from the gmx::AnalysisData object itself.
 * However, in many cases it is sufficient to initially add a plotting module
 * to the processing chain, which will then automatically write the results
 * into a file.
 *
 * For simple processing needs with a small amount of data, an
 * gmx::AnalysisArrayData class is also provided, which keeps all the data in an
 * in-memory array and allows you to manipulate the data as you wish before you
 * pass the data to the attached modules.
 *
 * <H3>Data Modules</H3>
 *
 * Modules that derive from gmx::AnalysisDataModuleInterface can operate in two
 * modes:
 *  - In _serial_ mode, the frames are presented to the module always in the
 *    order of increasing indices, even if they become ready in a different
 *    order in the attached data.
 *  - In _parallel_ mode, the frames are presented in the order that they
 *    become available in the input data, which may not be sequential.
 *    This mode allows the input data to optimize its behavior if it does not
 *    need to store and sort the frames.
 *
 * The figure below shows the sequence of callbacks that the module receives.
 * Arrows show a dependency between callbacks: the event at the start of the
 * arrow always occurs before the event at the end.  The events in the box are
 * repeated for each frame.  Dashed lines within this box show dependencies
 * between these frames:
 *  - In serial mode, all the events are called in a deterministic order, with
 *    each frame completely processed before the next starts.
 *  - In parallel mode, multiple frames can be in progress simultaneously, and
 *    the events for different frames can occur even concurrently on different
 *    threads.  However, frameFinishSerial() events will always occur in
 *    deterministic, sequential order for the frames.  Also, the number of
 *    concurrent frames is limited by the parallelization factor passed to
 *    parallelDataStarted(): only M frames after the last frame for which
 *    frameFinishSerial() has been called can be in progress
 *
 * \dot
 *     digraph datamodule_events {
 *         rankdir = LR
 *         node [ shape=box ]
 *
 *         start  [ label="dataStarted()",
 *                  URL="\ref gmx::AnalysisDataModuleInterface::dataStarted()" ]
 *         pstart [ label="parallelDataStarted()",
 *                  URL="\ref gmx::AnalysisDataModuleInterface::parallelDataStarted()" ]
 *         subgraph cluster_frame {
 *             label = "for each frame"
 *             framestart   [ label="frameStarted()",
 *                            URL="\ref gmx::AnalysisDataModuleInterface::frameStarted()" ]
 *             pointsadd    [ label="pointsAdded()",
 *                            URL="\ref gmx::AnalysisDataModuleInterface::pointsAdded()" ]
 *             framefinish  [ label="frameFinished()",
 *                            URL="\ref gmx::AnalysisDataModuleInterface::frameFinished()" ]
 *             serialfinish [ label="frameFinishedSerial()",
 *                            URL="\ref gmx::AnalysisDataModuleInterface::frameFinishedSerial()" ]
 *         }
 *         finish [ label="dataFinished()",
 *                  URL="\ref gmx::AnalysisDataModuleInterface::dataFinished()" ]
 *
 *         start -> framestart
 *         pstart -> framestart
 *         framestart -> pointsadd
 *         pointsadd -> pointsadd [ label="0..*", dir=back ]
 *         pointsadd -> framefinish
 *         framefinish -> serialfinish
 *         serialfinish -> finish
 *
 *         framestart:se -> serialfinish:sw [ dir=back, style=dashed, weight=0,
 *                                            label="serial: frame n+1\nparallel: frame n+M" ]
 *         serialfinish -> serialfinish [ dir=back, style=dashed,
 *                                        label="frame n+1" ]
 *     }
 * \enddot
 *
 * If the input data supports parallel mode, it calls parallelDataStarted().
 * If the module returns `true` from this method, then it will process the
 * frames in the parallel mode.  If the module returns `false`, it will get the
 * frames in serial order.
 * If the input data does not support parallel mode, it calls dataStarted(),
 * and the module will always get the frames in order.
 *
 * The sequence of when the module methods are called with respect to when data
 * is added to the data object depends on the type of the module and the type
 * of the data.  However, generally the modules do not need to know the details
 * of how this happens, as long as they work with the above state diagram.
 *
 * For parallel processing, the gmx::AnalysisData object itself only provides
 * the infrastructure to support all of the above, including the reordering of
 * the frames for serial processing.  However, the caller is still responsible
 * of the actual thread synchronization, and must call
 * gmx::AnalysisData::finishFrameSerial() for each frame from a suitable
 * context where the serial processing for that frame can be done.  When using
 * the data objects as part of the trajectory analysis framework
 * (\ref page_analysisframework), these calls are handled by the framework.
 *
 * \if libapi
 * <H3>Writing New Data and Module Objects</H3>
 *
 * New data modules can be implemented to perform custom operations that are
 * not supported by the modules provided in this module.  This is done by
 * creating a new class that implements gmx::AnalysisDataModuleInterface.
 * If the new module computes values that can be used as input for other
 * modules, the new class should also derive from gmx::AbstractAnalysisData, and
 * preferably use gmx::AnalysisDataStorage internally to implement storage of
 * values.  See the documentation of the mentioned classes for more details on
 * how to implement custom modules.
 * When implementing a new module, it should be considered whether it can be of
 * more general use, and if so, it should be added to this module.
 *
 * It is also possible to implement new data source objects by deriving a class
 * from gmx::AbstractAnalysisData.  This should not normally be necessary, since
 * this module provides general data source objects for most typical uses.
 * If the classes in this module are not suitable for some specific use, it
 * should be considered whether a new generic class could be added (or an
 * existing extended) instead of implementing a local custom solution.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for analysis data handling.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_H
#define GMX_ANALYSISDATA_H

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/arraydata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/displacement.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/lifetime.h"
#include "gromacs/analysisdata/modules/plot.h"

#endif
