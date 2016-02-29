/*
 * $Id: slater_low.h,v 1.2 2008/11/26 11:06:22 spoel Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 *
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _slater_low_h
#define _slater_low_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCLN
/* slater_integrals.cpp (c) 2008 Paul J. van Maaren and David van der Spoel */
#include <cln/cln.h>

using namespace cln;
#define PRECISION 80

static cl_R           ZERO      = "0.0_80";
static cl_R           ONE       = "1.0_80";
static cl_R           TWO       = "2.0_80";
static cl_R           THREE     = "3.0_80";
static cl_R           FOUR      = "4.0_80";
static cl_R           FIVE      = "5.0_80";
static cl_R           SIX       = "6.0_80";
static cl_R           SEVEN     = "7.0_80";
static cl_R           EIGHT     = "8.0_80";
static cl_R           NINE      = "9.0_80";
static float_format_t precision = float_format(80);

extern cl_R Power(cl_R a, int b);

extern cl_R Slater_1S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_1S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_2S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_3S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_4S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_5S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Slater_6S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_1S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_1S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_2S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_2S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_3S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_3S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_4S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_4S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_5S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_5S(cl_R r, cl_R xi, cl_R xj);

extern cl_R DSlater_6S_6S(cl_R r, cl_R xi, cl_R xj);

extern cl_R Nuclear_1S(cl_R r, cl_R xi);

extern cl_R Nuclear_2S(cl_R r, cl_R xi);

extern cl_R Nuclear_3S(cl_R r, cl_R xi);

extern cl_R Nuclear_4S(cl_R r, cl_R xi);

extern cl_R Nuclear_5S(cl_R r, cl_R xi);

extern cl_R Nuclear_6S(cl_R r, cl_R xi);

extern cl_R DNuclear_1S(cl_R r, cl_R xi);

extern cl_R DNuclear_2S(cl_R r, cl_R xi);

extern cl_R DNuclear_3S(cl_R r, cl_R xi);

extern cl_R DNuclear_4S(cl_R r, cl_R xi);

extern cl_R DNuclear_5S(cl_R r, cl_R xi);

extern cl_R DNuclear_6S(cl_R r, cl_R xi);

#endif

#endif
