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
/*! \defgroup module_coordinateio Handling of opening and closing coordinate files.
 * \brief
 * Provides basic functions to open and manipulate coordinate files during data analysis.
 *
 * <H3>Overview</H3>
 * The module provides the infrastructure to perform operations on coordinate data files
 * and structures during data analysis. It implements ways to change the information
 * in coordinate data structures as well as checking that both input data and output
 * method are matching for a given coordinate file writing operation.
 *
 * The module has as its main parts the outputadapter modules implemented using the
 * ICoordinateData interface to change information in t_trxframes,
 * as well as communicating the requirements to write the included information to disk.
 *
 * <H3>Outputadapter Modules</H3>
 * Each module derives from the common interfaces and has to set its requirements for final
 * processing as a flag from the enum in requirementflags. During processing, they implement a custom
 * version of the processFrame directive that modifies one part of the information in
 * a t_trxframe before the data gets written to disk.
 *
 *
 * The interaction between the OutputManager and the frameadapter modules derived from
 * ICoordinateOutput is shown in the diagram below.
 *
 * \msc
   wordwraparcs=true;

   analysistool,
   analysisloop,
   outputmanager,
   userinput,
   builder,
   modules,
   filewriting;

   builder box builder [ label="builds new OutputManager" ];

   analysistool => builder [ label="constructs" ];
   userinput => modules [ label="requires construction" ];
   modules box modules [ label="register requirements for output file" ];
   modules => builder [ label="registered on" ];
   builder => outputmanager [ label="creates new" ];
   modules => outputmanager [ label="check requirements" ];
   builder box builder [ label="new object created, not used again afterwards" ];
   --- [ label="Setup of outputmanager complete, analysis phase begins" ];
    analysistool => analysisloop [ label="starts iteration over frames" ];
    analysisloop => outputmanager [ label="provides coordinates" ];
    outputmanager => modules [ label="provide coordinate frames for changing" ];
    modules => outputmanager [ label="return after changing data" ];
    outputmanager => filewriting [ label="Send new coordinates for writing" ];
    outputmanager => analysisloop [ label="returns after writing" ];
    analysisloop box analysisloop [ label="iterates over frames" ];

 *  \endmsc
 *
 *
 *
 * \if libapi
 * <H3>Preparing new Modules</H3>
 *
 * To perform additional changes to a t_trxframe, new modules can be written that again
 * implement the ICoordinateData interface. The new method should
 * follow the approach of the other modules that are present in changing a minimal set of
 * t_trxframe data. The loading of the module can then be added to the CoordinateOutputUserOptions
 * class that acts as a wrapper that can be called from methods that want to use the
 * outputadapters.
 * \endif
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for file input/output.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_H
#define GMX_COORDINATEIO_H

#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/coordinateio/outputmanager.h"
#include "gromacs/coordinateio/requirementflags.h"
#include "gromacs/coordinateio/modules/outputselector.h"
#include "gromacs/coordinateio/modules/setatoms.h"
#include "gromacs/coordinateio/modules/setbox.h"
#include "gromacs/coordinateio/modules/setforces.h"
#include "gromacs/coordinateio/modules/setprecision.h"
#include "gromacs/coordinateio/modules/settime.h"
#include "gromacs/coordinateio/modules/setvelocities.h"

#endif
