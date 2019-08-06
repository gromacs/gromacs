/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Provides basic functions to handle writing of new coordinate files.
 *
 * The methods included in the coordinateio module implement the basics
 * for manipulating and writing coordinate trajectory files
 * and changing metadata in the underlying datastructures. Included are a container for storing
 * modules that change trajectory data, as well as a manager to write output files that uses
 * those methods. It can be used from within \ref module_trajectoryanalysis, and uses
 * methods from:
 * - \ref module_options
 * - \ref module_selection
 *
 * <H3>Overview</H3>
 * The methods in coordinateio provide the infrastructure to perform operations on coordinate data files
 * and structures during data analysis. It implements ways to change the information
 * in coordinate data structures as well as checking that both input data and output
 * method are matching for a given coordinate file writing operation. For this
 * components verify first that all the requirements can be satisfied. Then
 * components are build that will change the coordinate information accordingly.
 *
 * The main parts are the outputadapters implemented using the
 * IOutputAdapter interface to change information in a local (deep) copy of t_trxframes
 * stored in the coordinatefile.
 *
 * <H3>Outputadapter</H3>
 * Each OutputAdapter module implements the same IOutputAdapter interface and
 * has to set its requirements for final
 * processing as a flag from the enum in requirementflags. During processing, they implement a custom
 * version of the processFrame directive that modifies the stored trajectory data before writing
 * a new file to disk.
 *
 *
 * The interaction between the CoordinateFile and the OutputAdapter modules derived from
 * IOutputAdapter is shown in the diagram below.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   builder [ label="CoordinateFileBuilder" ],
   coordinatefile [ label="CoordinateFile" ],
   container [ label="OutputAdapterStorage" ],
   outputadapters [ label="OutputAdapters" ];

   analysistool => builder [ label="Requests new coordinate output" ];
   analysistool => builder [ label="Specifies required OutputAdapters" ];
   builder => outputadapters [ label="Tries to construct new outputadapters" ];
   outputadapters => builder [ label="Return or give error for wrong preconditions" ];
   outputadapters => container [ label="Outputadapters are stored" ];
   container => builder [ label="Gives error if storage conditions are violated" ];
   builder => coordinatefile [ label="Constructs new manager according to specifications" ];
   builder => container [ label="Requests storage object with registered outputadapters" ];
   container => builder [ label="Gives ownership of stored outputadapters" ];
   builder box builder [ label="Tests preconditions of storage object and new coordinatefile" ];
   builder => analysistool [ label="Raise error if preconditions don't match" ];
   builder => coordinatefile [ label="Add storage object to new coordinatefile" ];
   coordinatefile => analysistool [ label="Returns finished coordinatefile" ];
   builder box builder [ label="coordinatefile created, can start to work on input data" ];

 * \endmsc
 *
 * Once the CoordinateFile object and its registered modules are created, they can be used to
 * iterate over input data to write new coordinate frames.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   analysisloop,
   coordinatefile [ label="CoordinateFile" ],
   outputadapters [ label="OutputAdapters" ] ,
   filewriting;

   --- [ label="Setup of coordinatefile complete, analysis phase begins" ];
    analysistool   => analysisloop [ label="Starts iteration over frames" ];
    analysisloop   => coordinatefile [ label="Provides coordinates" ];
    coordinatefile  => outputadapters [ label="Provide coordinate frames for changing" ];
    outputadapters => coordinatefile [ label="Return after changing data" ];
    coordinatefile  => filewriting [ label="Send new coordinates for writing" ];
    filewriting    => coordinatefile [ label="Continue after writing to disk" ];
    coordinatefile  => analysisloop [ label="Returns after writing" ];
    analysisloop box analysisloop [ label="Iterates over frames" ];
    --- [ label="Analysis complete, object is destructed and files are closed" ];

 *  \endmsc
 *
 *
 *
 * \if libapi
 * <H3>Preparing new OutputAdapters</H3>
 *
 * If additional methods are needed to perform changes to the t_trxframe metadata,
 * new OutputAdapters can be written that again implement the IOutputAdapter interface.
 * The new method should follow the approach of the other modules that are present
 * in changing a minimal set of t_trxframe data.
 * \endif
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
/*! \file
 * \libinternal
 * \brief
 * Public API convenience header for coordinate file output.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_H
#define GMX_COORDINATEIO_H

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/coordinateio/outputadapters.h"
#include "gromacs/coordinateio/outputadapters/outputselector.h"
#include "gromacs/coordinateio/outputadapters/setatoms.h"
#include "gromacs/coordinateio/outputadapters/setbox.h"
#include "gromacs/coordinateio/outputadapters/setforces.h"
#include "gromacs/coordinateio/outputadapters/setprecision.h"
#include "gromacs/coordinateio/outputadapters/setstarttime.h"
#include "gromacs/coordinateio/outputadapters/settimestep.h"
#include "gromacs/coordinateio/outputadapters/setvelocities.h"

#include "coordinateio/enums.h"

#endif
