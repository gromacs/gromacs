/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief
 * Tests for functionality to do a table load for 2xnn kernels.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "base.h"

namespace SIMDTests
{

#ifdef GMX_NBNXN_SIMD_2XNN

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_2xnn_outer_header.h"

#include "load_table.h"

//! Typedef for the test fixture
typedef SimdFunctionDoingTableLoad_f SimdFunctionDoing2xnnTableLoad_f;

TEST_F(SimdFunctionDoing2xnnTableLoad_f, load_table_f_Works)
{
    Tester(ReferenceFunctions::load_table_f,
           TestFunctions::load_table_f);
}

//! Typedef for the test fixture
typedef SimdFunctionDoingTableLoad_f_v SimdFunctionDoing2xnnTableLoad_f_v;

TEST_F(SimdFunctionDoing2xnnTableLoad_f_v, load_table_f_v_Works)
{

    Tester(ReferenceFunctions::load_table_f_v,
           TestFunctions::load_table_f_v);
}

#endif

} // namespace
