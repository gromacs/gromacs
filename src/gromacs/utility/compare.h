/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Utilities for comparing data structures (for gmx check).
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_COMPARE_H
#define GMX_UTILITY_COMPARE_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

//! Compares two real values for equality.
gmx_bool equal_real(real i1, real i2, real ftol, real abstol);
//! Compares two float values for equality.
gmx_bool equal_float(float i1, float i2, float ftol, float abstol);
//! Compares two double values for equality.
gmx_bool equal_double(double i1, double i2, real ftol, real abstol);

//! Compares two integers and prints differences.
void cmp_int(FILE *fp, const char *s, int index, int i1, int i2);

//! Compares two 64-bit integers and prints differences.
void cmp_int64(FILE *fp, const char *s, gmx_int64_t i1, gmx_int64_t i2);

//! Compares two unsigned short values and prints differences.
void cmp_us(FILE *fp, const char *s, int index, unsigned short i1, unsigned short i2);

//! Compares two unsigned char values and prints differences.
void cmp_uc(FILE *fp, const char *s, int index, unsigned char i1, unsigned char i2);

//! Compares two boolean values and prints differences, and returns whether both are true.
gmx_bool cmp_bool(FILE *fp, const char *s, int index, gmx_bool b1, gmx_bool b2);

//! Compares two strings and prints differences.
void cmp_str(FILE *fp, const char *s, int index, const char *s1, const char *s2);

//! Compares two reals and prints differences.
void cmp_real(FILE *fp, const char *s, int index, real i1, real i2, real ftol, real abstol);

//! Compares two floats and prints differences.
void cmp_float(FILE *fp, const char *s, int index, float i1, float i2, float ftol, float abstol);

//! Compares two doubles and prints differences.
void cmp_double(FILE *fp, const char *s, int index, double i1, double i2, double ftol, double abstol);

#endif
