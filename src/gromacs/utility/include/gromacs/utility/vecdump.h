/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#ifndef GMX_UTILITY_VECDUMP_H
#define GMX_UTILITY_VECDUMP_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/vectypes.h"


namespace gmx
{

template<typename>
class ArrayRef;
class TextWriter;

//! \{
/*! \brief Dump the formatted \c values to \c writer along with a description and/or indices. */
void dumpIntArrayRef(TextWriter* writer, const char* description, ArrayRef<const int> values, bool showIndices);
void dumpFloatArrayRef(TextWriter* writer, const char* description, ArrayRef<const float> values, bool showIndices);
void dumpDoubleArrayRef(TextWriter* writer, const char* description, ArrayRef<const double> values, bool showIndices);
void dumpRealArrayRef(TextWriter* writer, const char* description, ArrayRef<const real> values, bool showIndices);
void dumpRvecArrayRef(TextWriter* writer, const char* description, ArrayRef<const RVec> values);
//! \}

} // namespace gmx

void pr_ivec(FILE* fp, int indent, const char* title, const int vec[], int n, gmx_bool bShowNumbers);
void pr_ivec_block(FILE* fp, int indent, const char* title, const int vec[], int n, gmx_bool bShowNumbers);
void pr_ivecs(FILE* fp, int indent, const char* title, const ivec vec[], int n, gmx_bool bShowNumbers);
void pr_rvec(FILE* fp, int indent, const char* title, const real vec[], int n, gmx_bool bShowNumbers);
void pr_fvec(FILE* fp, int indent, const char* title, const float vec[], int n, gmx_bool bShowNumbers);
void pr_dvec(FILE* fp, int indent, const char* title, const double vec[], int n, gmx_bool bShowNumbers);
void pr_rvecs(FILE* fp, int indent, const char* title, const rvec vec[], int n);

//! New interface receiving an ArrayRef as input
void prRVecs(FILE* fp, int indent, const char* title, gmx::ArrayRef<const gmx::RVec> vec);

#endif
