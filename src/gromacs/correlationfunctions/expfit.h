/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/*! \libinternal
 * \file
 * \brief
 * Declares routine for fitting a data set to a curve
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_EXPFIT_H
#define GMX_EXPFIT_H

#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Enum to select fitting functions
 */
enum {
    effnNONE, effnEXP1, effnEXP2, effnEXPEXP,
    effnEXP5, effnEXP7, effnEXP9,
    effnVAC,  effnERF,  effnERREST, effnPRES, effnNR
};

/*! \brief
 * Short name of each function type.
 * This is exported for now in order to use when
 * calling parse_common_args.
 */
extern const char *s_ffn[effnNR+2];

/*! \brief
 * Returns  description corresponding to the enum above, or NULL if out of range
 * \param[in] effn Index
 * \return Description or NULL
 */
const char *effnDescription(int effn);

/*! \brief
 * Returns  number of function parameters associated with a fitting function.
 * \param[in] effn Index
 * \return number or -1 if index out of range
 */
int effnNparams(int effn);

/*! \brief
 * Returns  corresponding to the selected enum option in sffn
 * \param[in] sffn Two dimensional string array coming from parse_common_args
 * \return the ffn enum
 */
int sffn2effn(const char **sffn);

/*! \brief
 * Returns the value of fit function eFitFn at x
 * \param[in] eFitFn the index to the fitting function (0 .. effnNR)
 * \param[in] parm Array of parameters, the length of which depends on eFitFn
 * \param[in] x The value of x
 * \return the value of the fit
 */
double fit_function(const int eFitFn, const double parm[], const double x);

/*! \brief
 * Use Levenberg-Marquardt method to fit to a nfitparm parameter exponential
 * or to a transverse current autocorrelation function.
 *
 * If x == NULL, the timestep dt will be used to create a time axis.
 * \param[in] ndata Number of data points
 * \param[in] c1 The data points
 * \param[in] sig The standard deviation in the points (can be NULL)
 * \param[in] dt The time step
 * \param[in] x The X-axis (may be NULL, see above)
 * \param[in] begintimefit Starting time for fitting
 * \param[in] endtimefit Ending time for fitting
 * \param[in] oenv Output formatting information
 * \param[in] bVerbose Should the routine write to console?
 * \param[in] eFitFn Fitting function (0 .. effnNR)
 * \param[inout] fitparms[] Fitting parameters, see printed manual for a
 * detailed description. Note that in all implemented cases the parameters
 * corresponding to time constants will be generated with increasing values.
 * Such input parameters should therefore be provided in increasing order.
 * If this is not the case or if subsequent time parameters differ by less than
 * a factor of 2, they will be modified to ensure tau_i+1 >= 2 tau_i.
 * \param[in] fix Constrains fit parameter i at it's starting value, when the i'th bit
 * of fix is set. This works only when the N last parameters are fixed
 * but not when a parameter somewhere in the middle needs to be fixed.
 * \param[in] fn_fitted If not NULL file to print the data and fitted curve to
 * \return integral.
 */
real do_lmfit(int ndata, real c1[], real sig[], real dt, real *x,
              real begintimefit, real endtimefit, const output_env_t oenv,
              gmx_bool bVerbose, int eFitFn, double fitparms[], int fix,
              const char *fn_fitted);

/*! \brief
 * Fit an autocorrelation function to a pre-defined functional form
 *
 * \todo check parameters
 * \param[in] ncorr
 * \param[in] fitfn Fitting function (0 .. effnNR)
 * \param[in] oenv Output formatting information
 * \param[in] bVerbose Should the routine write to console?
 * \param[in] tbeginfit Starting time for fitting
 * \param[in] tendfit Ending time for fitting
 * \param[in] dt The time step
 * \param[in] c1 The data points
 * \param[inout] fit The fitting parameters
 * \return the integral over the autocorrelation function?
 */
real fit_acf(int ncorr, int fitfn, const output_env_t oenv, gmx_bool bVerbose,
             real tbeginfit, real tendfit, real dt, real c1[], real *fit);

#ifdef __cplusplus
}
#endif

#endif
