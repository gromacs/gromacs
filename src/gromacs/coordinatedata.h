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
/*! \defgroup module_coordinatedata Handling of Coordinate Frame Data (coordinatedata)
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for handling and processing
 * coordinate frame data for analysis.
 *
 * <H3>Overview</H3>
 * This module provides the framework for frameconverter modules that have the
 * task to operate on coordinate frame data structures before writing them to disk
 * or passing them on to other analysis tools.
 *
 * The module provides the IFrameConverter interface to write tools that change the
 * coordinates of t_trxframe datastructures according to the input provided by the user
 * when constructing the corresponding object.
 *
 * An additional method is provided to chain different modules together and allow
 * the processing of several different changes one after the other.
 *
 * <H3>Frameconverter Modules</H3>
 * The common modules derive from the interface mentioned above and expose functions to
 * set a coordinate frame for modification, process the data inside the frame object, and
 * return the frame to be worked on by other tools.
 *
 * \if libapi
 * <H3>Writing new frameconverters</H3>
 * A new frameconverter should implement the IFrameConverter interface and should follow
 * the approach shown in the existing examples to change the coordinates in a frame
 * in a way that can be passed on if needed to more objects afterwards.
 * \endif
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for coordinate data handling.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinatedata
 */
#ifndef GMX_COORDINATEDATA_H
#define GMX_COORDINATEDATA_H

#include "gromacs/coordinatedata/frameconverters.h"

#endif
