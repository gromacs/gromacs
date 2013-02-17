/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012, by the GROMACS development team, led by
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
/*! \defgroup module_analysisdata Parallelizable Handling of Output Data
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
 * AbstractAnalysisData class and classes that derive from it:
 * AnalysisData and AnalysisArrayData.  These classes are used to process and
 * store raw data as produced by the analysis tool.  They also provide an
 * interface to attach data modules that implement AnalysisDataModuleInterface.
 * Modules that implement this interface form the second part of the module,
 * and they provide functionality to do processing operations on the data.
 * These modules can also derive from AbstractAnalysisData, allowing other
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
 * AnalysisDataModuleInterface) and attach them to the data objects or to one
 * another using AbstractAnalysisData::addModule().  Then you add the actual
 * data values to the AnalysisData object, which automatically passes it on to
 * the modules.  After all data is added, you may optionally access some
 * results directly from the module objects.  However, in many cases it is
 * sufficient to initially add a plotting module to the processing chain, which
 * will then automatically write the results into a file.
 *
 * For simple processing needs with a small amount of data, an
 * AnalysisArrayData class is also provided, which keeps all the data in an
 * in-memory array and allows you to manipulate the data as you wish before you
 * pass the data to the attached modules.
 *
 * \if libapi
 * <H3>Writing New Data and Module Objects</H3>
 *
 * New data modules can be implemented to perform custom operations that are
 * not supported by the modules provided in this module.  This is done by
 * creating a new class that implements AnalysisDataModuleInterface.
 * If the new module computes values that can be used as input for other
 * modules, the new class should also derive from AbstractAnalysisData, and
 * preferably use AnalysisDataStorage internally to implement storage of
 * values.  See the documentation of the mentioned classes for more details on
 * how to implement custom modules.
 * When implementing a new module, it should be considered whether it can be of
 * more general use, and if so, it should be added to this module.
 *
 * It is also possible to implement new data source objects by deriving a class
 * from AbstractAnalysisData.  This should not normally be necessary, since
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
#include "analysisdata/modules/plot.h"

#endif
