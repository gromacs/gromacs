/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2013, by the GROMACS development team, led by
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
 *
 * This module consists of two main parts.  The first is formed by the
 * gmx::AbstractAnalysisData class and classes that derive from it:
 * gmx::AnalysisData and gmx::AnalysisArrayData.  These classes are used to
 * process and store raw data as produced by the analysis tool.  They also
 * provide an interface to attach data modules that implement
 * gmx::AnalysisDataModuleInterface.
 * Modules that implement this interface form the second part of the module,
 * and they provide functionality to do processing operations on the data.
 * These modules can also derive from gmx::AbstractAnalysisData, allowing other
 * modules to be attached to them to form a processing chain that best suits
 * the analysis tool.  Typically, such a processing chain ends in a plotting
 * module that writes the data into a file, but the final module can also
 * provide direct access to the processed data, allowing the analysis tool to
 * do custom postprocessing outside the module framework.
 *
 * The sequence chart below shows an overview of how analysis data objects
 * and modules interact (currently, multipoint data is an exception to this
 * diagram).
 * \msc
 *     caller,
 *     data [ label="AnalysisData", URL="\ref gmx::AnalysisData" ],
 *     module1 [ label="data module", URL="\ref gmx::AnalysisDataModuleInterface" ];
 *
 *     caller box module1 [ label="caller creates and initializes all objects" ];
 *     caller => data [ label="addModule(module1)",
 *                      URL="\ref gmx::AbstractAnalysisData::addModule() " ];
 *     caller => data [ label="startData()",
 *                      URL="\ref gmx::AnalysisData::startData()" ];
 *     data => module1 [ label="dataStarted()",
 *                       URL="\ref gmx::AnalysisDataModuleInterface::dataStarted()" ];
 *     caller << data [ label="handle",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => data [ label="handle->startFrame(0)",
 *                      URL="\ref gmx::AnalysisDataHandle::startFrame()" ];
 *     caller => data [ label="add data for frame 0",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => data [ label="handle->finishFrame() (frame 0)",
 *                      URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     data => module1 [ label="frameStarted(0)",
 *                       URL="\ref gmx::AnalysisDataModuleInterface::frameStarted()" ];
 *     data => module1 [ label="pointsAdded()",
 *                       URL="\ref gmx::AnalysisDataModuleInterface::pointsAdded()" ];
 *     data => module1 [ label="frameFinished(0)",
 *                       URL="\ref gmx::AnalysisDataModuleInterface::frameFinished()" ];
 *     caller => data [ label="handle->startFrame(1)",
 *                      URL="\ref gmx::AnalysisDataHandle::startFrame()" ];
 *     caller => data [ label="add data for frame 1",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => data [ label="handle2->finishFrame() (frame 1)",
 *                      URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     data => module1 [ label="process frame 1" ];
 *     ... [ label="add more frames" ];
 *     caller => data [ label="handle->finishData()",
 *                      URL="\ref gmx::AnalysisDataHandle::finishData()" ];
 *     data => module1 [ label="dataFinished()",
 *                       URL="\ref gmx::AnalysisDataModuleInterface::dataFinished()" ];
 * \endmsc
 *
 * The second sequence chart below shows the interaction in the case of two
 * gmx::AnalysisDataHandle objects, which can be used to insert data
 * concurrently for multiple frames.  The gmx::AnalysisData object and the
 * handles take care to notify the module such that it always receives the
 * frames in the correct order.
 * \msc
 *     caller,
 *     handle1 [ label="handle1", URL="\ref gmx::AnalysisDataHandle" ],
 *     handle2 [ label="handle2", URL="\ref gmx::AnalysisDataHandle" ],
 *     module1 [ label="data module", URL="\ref gmx::AnalysisDataModuleInterface" ];
 *
 *     caller box handle2 [ label="caller has created both handles using startData()" ];
 *     caller => handle1 [ label="startFrame(0)",
 *                         URL="\ref gmx::AnalysisDataHandle::startFrame()" ];
 *     caller => handle2 [ label="startFrame(1)",
 *                         URL="\ref gmx::AnalysisDataHandle::startFrame()" ];
 *     caller => handle1 [ label="add data for frame 0",
 *                         URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => handle2 [ label="add data for frame 1",
 *                         URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => handle2 [ label="finishFrame() (frame 1)",
 *                         URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     caller => handle2 [ label="startFrame(2)",
 *                         URL="\ref gmx::AnalysisDataHandle::startFrame()" ];
 *     caller => handle1 [ label="finishFrame() (frame 0)",
 *                         URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     handle1 => module1 [ label="process frame 0" ];
 *     handle1 => module1 [ label="process frame 1" ];
 *     caller => handle2 [ label="add data for frame 2",
 *                         URL="\ref gmx::AnalysisDataHandle" ];
 *     caller => handle2 [ label="finishFrame() (frame 2)",
 *                         URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     handle2 => module1 [ label="process frame 2" ];
 *     ...;
 *     caller => handle1 [ label="finishData()",
 *                         URL="\ref gmx::AnalysisDataHandle::finishData()" ];
 *     caller => handle2 [ label="finishData()",
 *                         URL="\ref gmx::AnalysisDataHandle::finishData()" ];
 *     handle2 => module1 [ label="dataFinished()",
 *                          URL="\ref gmx::AnalysisDataModuleInterface::dataFinished()" ];
 * \endmsc
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
 * After all data is added, you may optionally access some
 * results directly from the module objects.  However, in many cases it is
 * sufficient to initially add a plotting module to the processing chain, which
 * will then automatically write the results into a file.
 *
 * For simple processing needs with a small amount of data, an
 * gmx::AnalysisArrayData class is also provided, which keeps all the data in an
 * in-memory array and allows you to manipulate the data as you wish before you
 * pass the data to the attached modules.
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

#include "analysisdata/analysisdata.h"
#include "analysisdata/arraydata.h"
#include "analysisdata/dataframe.h"
#include "analysisdata/modules/average.h"
#include "analysisdata/modules/displacement.h"
#include "analysisdata/modules/histogram.h"
#include "analysisdata/modules/lifetime.h"
#include "analysisdata/modules/plot.h"

#endif
