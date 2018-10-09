/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \defgroup module_coordinateio Handling of writing new coordinate files
 * \ingroup group_analysismodules
 * \brief
 * Provides basic functions to to open coordinate files and manipulate structures
 * read from them during data analysis.
 *
 * The module implements the basics for handling of writing coordinate trajectory files
 * and changing metadata in the underlying datastructure. It provides a container for storing
 * modules that change metadata, as well as a manager to write output files that uses
 * those methods. It can be used from within \ref module_trajectoryanalysis, and uses
 * methods from:
 * - \ref module_options
 * - \ref module_selection
 *
 * <H3>Overview</H3>
 * The module provides the infrastructure to perform operations on coordinate data files
 * and structures during data analysis. It implements ways to change the information
 * in coordinate data structures as well as checking that both input data and output
 * method are matching for a given coordinate file writing operation.
 *
 * The module has as its main parts the outputadapter modules implemented using the
 * ICoordinateData interface to change information in a local (deep) copy of t_trxframes
 * stored in the outputmanager,
 * as well as communicating the requirements to write the included information to disk.
 *
 * <H3>Outputadapter</H3>
 * Each OutputAdapter module implements the same IOutputAdapter interface and
 * has to set its requirements for final
 * processing as a flag from the enum in requirementflags. During processing, they implement a custom
 * version of the processFrame directive that modifies one part of the information in
 * the (previously deep copied) t_trxframe  local to the outputmanager
 * before the data gets written to disk.
 *
 *
 * The interaction between the OutputManager and the OutputAdapter modules derived from
 * IOutputAdapter is shown in the diagram below.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   userinput,
   builder [ label="OutputManagerBuilder" ],
   outputmanager [ label="OutputManager" ],
   outputadapters [ label="OutputAdapters" ],
   container [ label="OutputAdapterStorage" ];

   builder box builder [ label="Builds new OutputManager according to requirements" ];
   userinput => builder [ label="Requests new coordinate output" ];
   userinput => outputadapters [ label="Specifies required modules" ];
   outputadapters => userinput [ label="Return or give error for wrong preconditions" ];
   outputadapters => container [ label="Stores outputadapters" ];
   builder => outputmanager [ label="Constructs new manager according to specifications" ];
   builder => container [ label="Requests outputadapters from storage" ];
   container => builder [ label="Gives ownership of stored outputadapters" ];
   builder => outputadapters [ label="Tries to registers adapters" ];
   outputadapters => outputmanager [ label="Communicates module requirements" ];
   outputmanager => outputadapters [ label="Registers or throws error" ];
   outputmanager => builder [ label="Returns finished outputmanager" ];
   builder box builder [ label="outputmanager created, can start to work on input data" ];

 * \endmsc
 *
 * Once the OutputManager object and its registered modules are created, they can be used to
 * iterate over input data to write new coordinate frames.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   analysisloop,
   outputmanager [ label="OutputManager" ],
   builder [ label="OutputManagerBuilder" ],
   outputadapters [ label="OutputAdapters" ] ,
   filewriting;

   builder box builder [ label="Builds new OutputManager" ];

   analysistool => builder [ label="Requests new outputmanager" ];
   builder      => analysistool [ label="Returns constructed outputmanager" ];
   --- [ label="Setup of outputmanager complete, analysis phase begins" ];
    analysistool   => analysisloop [ label="Starts iteration over frames" ];
    analysisloop   => outputmanager [ label="Provides coordinates" ];
    outputmanager  => outputadapters [ label="Provide coordinate frames for changing" ];
    outputadapters => outputmanager [ label="Return after changing data" ];
    outputmanager  => filewriting [ label="Send new coordinates for writing" ];
    filewriting    => outputmanager [ label="Continue after writing to disk" ];
    outputmanager  => analysisloop [ label="Returns after writing" ];
    analysisloop box analysisloop [ label="Iterates over frames" ];
    --- [ label="Analysis complete, object is destructed and files are closed" ];

 *  \endmsc
 *
 *
 *
 * \if libapi
 * <H3>Preparing new Modules</H3>
 *
 * If additional methods are needed to perform changes to the t_trxframe metadata,
 * new modules can be written that again implement the ICoordinateData interface. The new method should
 * follow the approach of the other modules that are present in changing a minimal set of
 * t_trxframe data.
 *
 * The new module should then be registered on the outputmanager in the same way as the existing
 * modules, either through a wrapper that passes the module on to the outputmanagerbuilder, or through
 * a new custom method.
 * \endif
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for coordinate file output.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_H
#define GMX_COORDINATEIO_H

#include "gromacs/coordinateio/builder.h"
#include "gromacs/coordinateio/enums.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/coordinateio/outputmanager.h"

#endif
