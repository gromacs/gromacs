/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_selection Parsing and Evaluation of Analysis Selections (selection)
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for initializing and evaluating selections.
 *
 * The core of the selection engine is accessed through
 * gmx::SelectionCollection, which manages a set of selections.
 * Documentation for that class explains the general selection mechanisms.
 *
 * For each selection that is parsed using a gmx::SelectionCollection, a
 * gmx::Selection handle is returned and can be used to access information
 * about that selection.  gmx::SelectionPosition is a helper class used to
 * access information about individual positions in a selection.  These classes
 * refer to internal state within the gmx::SelectionCollection, and their
 * contents update automatically when the gmx::SelectionCollection is compiled
 * or evaluated.
 *
 * This module also provides gmx::SelectionOption and gmx::SelectionOptionInfo
 * classes for declaring options that evaluate to selections (see \ref
 * module_options for general explanation of the options mechanism).  These
 * classes provide the main interface to obtain gmx::Selection objects in
 * trajectory analysis using gmx::TrajectoryAnalysisModule.
 * To use these classes outside the trajectory analysis framework,
 * a gmx::SelectionOptionManager needs to be created to serve as a bridge
 * between the selection option classes and the gmx::SelectionCollection
 * object.
 * \if libapi
 * gmx::SelectionFileOption can be used to implement generic file input for
 * selection options (done internally in the trajectory analysis framework).
 *
 * The selection module contains some lower-level functionality that is
 * currently internal to it (centerofmass.h, indexutil.h, poscalc.h,
 * position.h), but could possibly be useful also outside the module.
 * It should be considered whether they should be moved somewhere else.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for selection handling.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_H
#define GMX_SELECTION_H

#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionmanager.h"

#endif
