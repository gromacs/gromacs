/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \defgroup module_selection Parsing and Evaluation of Analysis Selections
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for initializing and evaluating selections.
 *
 * The core of the selection engine is accessed through SelectionCollection,
 * which manages a set of selections.  Documentation for that class explains
 * the general selection mechanisms.
 *
 * For each selection that is parsed using a SelectionCollection, a \ref
 * Selection handle is returned and can be used to access information about
 * that selection.  SelectionPosition is a helper class used to access
 * information about individual positions in a selection.  These classes refer
 * to internal state within the SelectionCollection, and their contents update
 * automatically when the SelectionCollection is compiled or evaluated.
 *
 * This module also provides SelectionOption and SelectionOptionInfo classes
 * for declaring options that evaluate to selections (see \ref module_options
 * for general explanation of the options mechanism).  These classes provide
 * the main interface to obtain \ref Selection objects in trajectory analysis
 * using TrajectoryAnalysisModule.
 *
 * \if libapi
 * The selection module contains some lower-level functionality that is
 * currently internal to it (centerofmass.h, indexutil.h, poscalc.h,
 * position.h), but could possibly be useful also outside the module.
 * It should be considered whether they should be moved somewhere else.
 * \endif
 *
 * \if internal
 * Implementation details of different parts of the module are discussed on
 * separate pages:
 *   - \ref page_module_selection_custom
 *   - \ref page_module_selection_parser
 *   - \ref page_module_selection_compiler
 *   - \ref page_module_selection_insolidangle
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for selection handling.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_H
#define GMX_SELECTION_H

#include "selection/selection.h"
#include "selection/selectioncollection.h"
#include "selection/selectionoption.h"
#include "selection/selectionoptioninfo.h"

#endif
