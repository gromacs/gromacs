/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * \defgroup module_correlationfunctions Correlation functions
 * \ingroup group_analysismodules
 * \brief
 * Compute correlation functions and fit analytical functions to the result.
 */
/*! \libinternal
 * \file
 * \brief
 * Declares routine for computing autocorrelation functions
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_AUTOCORR_H
#define GMX_AUTOCORR_H

#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Normal correlation f(t)*f(t+dt) */
#define eacNormal (1<<0)
/*! \brief Cosine correlation cos(f(t)-f(t+dt)) */
#define eacCos    (1<<1)
/*! \brief Vector correlation f(t).f(t+dt) */
#define eacVector (1<<2)
/*! \brief Norm of cross product |f(t) (x) f(t+dt)| */
#define eacRcross (1<<3  | eacVector)
/*! \brief Vector with Legendre polynomial of order 0 (same as vector) */
#define eacP0     (1<<4  | eacVector)
/*! \brief Vector with Legendre polynomial of order P_1(f(t).f(t+dt)) */
#define eacP1     (1<<5  | eacVector)
/*! \brief Vector with Legendre polynomial of order P_2(f(t).f(t+dt)) */
#define eacP2     (1<<6  | eacVector)
/*! \brief Vector with Legendre polynomial of order P_3(f(t).f(t+dt)) */
#define eacP3     (1<<7  | eacVector)
/*! \brief Vector with Legendre polynomial of order P_4(f(t).f(t+dt)) */
#define eacP4     (1<<8  | eacVector)
/*! \brief Binary identy correlation (f(t) == f(t+dt)) */
#define eacIden   (1<<9) //Not supported for multiple cores

/*! \brief
 * Add commandline arguments related to autocorrelations to the existing array.
 * *npargs must be initialised to the number of elements in pa,
 * it will be incremented appropriately.
 *
 * \param npargs The number of arguments before and after (is changed in this function)
 * \param[in] pa The initial argument list
 * \return the new array
 */
t_pargs *add_acf_pargs(int *npargs, t_pargs *pa);

/*! \brief
 * Returns the number of points to output from a correlation function.
 * Works only AFTER do_auto_corr has been called!
 * \return the output length for the correlation function
 */
int get_acfnout(void);

/*! \brief
 * Returns the fitting function selected.
 * Works only AFTER do_auto_corr has been called!
 * \return the fit function type.
 */
int get_acffitfn(void);

/*! \brief
 * Calls low_do_autocorr (see below). add_acf_pargs has to be called before this
 * can be used.
 * \param[in] fn File name for xvg output (may this be NULL)?
 * \param[in] oenv The output environment information
 * \param[in] title is the title in the output file
 * \param[in] nframes is the number of frames in the time series
 * \param[in] nitem is the number of items
 * \param[inout] c1 is an array of dimension [ 0 .. nitem-1 ] [ 0 .. nframes-1 ]
 *          on output, this array is filled with the correlation function
 *          to reduce storage
 * \param[in] dt is the time between frames
 * \param[in] mode Different types of ACF can be done, see above
 * \param[in] bAver    If set, all ndih C(t) functions are averaged into a single
 *          C(t)
 */
void do_autocorr(const char *fn, const output_env_t oenv,
                 const char *title,
                 int nframes, int nitem, real **c1,
                 real dt, unsigned long mode, gmx_bool bAver);

/*! \brief
 * Low level computation of autocorrelation functions
 *
 * do_autocorr calculates autocorrelation functions for many things.
 * It takes a 2 d array containing nitem arrays of length nframes
 * for each item the ACF is calculated.
 *
 * A number of "modes" exist for computation of the ACF controlled
 * by variable mode, with the following meaning.
 *
 * Mode       | Function
 * -----------|------------
 * eacNormal  | C(t) = < X (tau) * X (tau+t) >
 * eacCos     | C(t) = < cos (X(tau) - X(tau+t)) >
 * eacIden    | C(t) = < (X(tau) == X(tau+t)) > (not fully supported yet)
 * eacVector  | C(t) = < X(tau) * X(tau+t)
 * eacP1      | C(t) = < cos (X(tau) * X(tau+t) >
 * eacP2      | C(t) = 1/2 * < 3 cos (X(tau) * X(tau+t) - 1 >
 * eacRcross  | C(t) = < ( X(tau) * X(tau+t) )^2 >
 *
 * For modes eacVector, eacP1, eacP2 and eacRcross the input should be
 * 3 x nframes long, where each triplet is taken as a 3D vector
 *
 * For mode eacCos inputdata must be in radians, not degrees!
 *
 * Other parameters are:
 *
 * \param[in] fn is output filename (.xvg) where the correlation function(s) are printed
 * \param[in] oenv controls output file properties
 * \param[in] title is the title in the output file
 * \param[in] nframes is the number of frames in the time series
 * \param[in] nitem is the number of items
 * \param[in] nout
 * \param[inout] c1 is an array of dimension [ 0 .. nitem-1 ] [ 0 .. nframes-1 ]
 *          on output, this array is filled with the correlation function
 *          to reduce storage
 * \param[in] dt is the time between frames
 * \param[in] mode Different types of ACF can be done, see above
 * \param[in] nrestart     is the number of steps between restarts for direct ACFs
 *              (i.e. without FFT) When set to 1 all points are used as
 *              time origin for averaging
 * \param[in] bAver    If set, all ndih C(t) functions are averaged into a single
 *          C(t)
 * \param[in] bNormalize   If set, all ACFs will be normalized to start at 0
 * \param[in] bVerbose If set output to console will be generated
 * \param[in] tbeginfit Time to start fitting to the ACF
 * \param[in] tendfit Time to end fitting to the ACF
 * \param[in] nfitparm Number of fitting parameters in a multi-exponential fit
 */
void low_do_autocorr(const char *fn, const output_env_t oenv,
                     const char *title, int  nframes, int nitem,
                     int nout, real **c1, real dt, unsigned long mode,
                     int nrestart, gmx_bool bAver, gmx_bool bNormalize,
                     gmx_bool bVerbose, real tbeginfit, real tendfit,
                     int nfitparm);

#ifdef __cplusplus
}
#endif

#endif
