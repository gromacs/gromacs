/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal
 * \file
 * \brief
 * Declares a driver routine for lmfit.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_CORRELATION_FUNCTIONS_GMX_LMCURVE_H
#define GMX_CORRELATION_FUNCTIONS_GMX_LMCURVE_H
#include "gromacs/correlationfunctions/expfit.h"
/*! \brief function type for passing to fitting routine */
typedef double (*t_lmcurve)(double x, const double* a);
/*! \brief lmfit_exp supports fitting of different functions
 *
 * This routine calls the Levenberg-Marquardt non-linear fitting
 * routine for fitting a data set with errors to a target function.
 * Fitting routines included in gromacs in src/external/lmfit.
 */
bool lmfit_exp(int          nfit,
               const double x[],
               const double y[],
               const double dy[],
               double       parm[],
               bool         bVerbose,
               int          eFitFn,
               int          nfix);

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
extern t_lmcurve lmcurves[effnNR + 1];

#endif
