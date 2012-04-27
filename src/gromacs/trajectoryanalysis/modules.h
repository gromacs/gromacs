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
/*! \file
 * \brief
 * Generic interface for creation of trajectory analysis modules.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_H
#define GMX_TRAJECTORYANALYSIS_MODULES_H

#include "analysismodule.h"

namespace gmx
{

/*! \brief
 * Creates a TrajectoryAnalysisModule object corresponding to a name.
 *
 * \param[in]  name  Name of the module to create (recognized names are
 *      defined in modules.h).
 * \returns  An allocated TrajectoryAnalysisModule object.
 * \throws   InvalidInputError if \p name is not recognized.
 *
 * This function should be used to instantiate analysis methods defined in the
 * library.
 *
 * In addition to recognizing exact matches on \p name, the function also
 * identifies cases where \p name is a prefix of exactly one recognized name
 * (exact matches always take precedence).
 *
 * \inpublicapi
 */
TrajectoryAnalysisModulePointer
createTrajectoryAnalysisModule(const char *name);

namespace analysismodules
{

static const char * const angle    = "angle";
static const char * const distance = "distance";
static const char * const select   = "select";

} // namespace modules

} // namespace gmx

#endif
