/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

#ifndef _nbnxn_kernel_simd_common_code_h
#define _nbnxn_kernel_simd_common_code_h

#if   defined(GMX_NBNXN_SIMD_2XNN)

#include "simd_2xnn/nbnxn_kernel_simd_2xnn_includes.h"

#define NBK_FN(elec, ljcomb, enertype) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _ ## enertype

#elif defined(GMX_NBNXN_SIMD_4XN)

#include "simd_4xn/nbnxn_kernel_simd_4xn_includes.h"

#define NBK_FN(elec, ljcomb, enertype) nbnxn_kernel_simd_4xn_ ## elec ## _comb_ ## ljcomb ## _ ## enertype

#else

#error "Unknown NBNXN SIMD type"

#endif

/*! \brief Kinds of electrostatic treatments in SIMD Verlet kernels
 */
enum {
    coultRF, coultTAB, coultTAB_TWIN, coultEWALD, coultEWALD_TWIN, coultNR
};

static p_nbk_func_ener p_nbk_ener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom, ener), NBK_FN(rf, lb, ener), NBK_FN(rf, none, ener) },
  { NBK_FN(tab, geom, ener), NBK_FN(tab, lb, ener), NBK_FN(tab, none, ener) },
  { NBK_FN(tab_twin, geom, ener), NBK_FN(tab_twin, lb, ener), NBK_FN(tab_twin, none, ener) },
  { NBK_FN(ewald, geom, ener), NBK_FN(ewald, lb, ener), NBK_FN(ewald, none, ener) },
  { NBK_FN(ewald_twin, geom, ener), NBK_FN(ewald_twin, lb, ener), NBK_FN(ewald_twin, none, ener) } };

static p_nbk_func_ener p_nbk_energrp[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom, energrp), NBK_FN(rf, lb, energrp), NBK_FN(rf, none, energrp) },
  { NBK_FN(tab, geom, energrp), NBK_FN(tab, lb, energrp), NBK_FN(tab, none, energrp) },
  { NBK_FN(tab_twin, geom, energrp), NBK_FN(tab_twin, lb, energrp), NBK_FN(tab_twin, none, energrp) },
  { NBK_FN(ewald, geom, energrp), NBK_FN(ewald, lb, energrp), NBK_FN(ewald, none, energrp) },
  { NBK_FN(ewald_twin, geom, energrp), NBK_FN(ewald_twin, lb, energrp), NBK_FN(ewald_twin, none, energrp) } };

static p_nbk_func_noener p_nbk_noener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom, noener), NBK_FN(rf, lb, noener), NBK_FN(rf, none, noener) },
  { NBK_FN(tab, geom, noener), NBK_FN(tab, lb, noener), NBK_FN(tab, none, noener) },
  { NBK_FN(tab_twin, geom, noener), NBK_FN(tab_twin, lb, noener), NBK_FN(tab_twin, none, noener) },
  { NBK_FN(ewald, geom, noener), NBK_FN(ewald, lb, noener), NBK_FN(ewald, none, noener) },
  { NBK_FN(ewald_twin, geom, noener), NBK_FN(ewald_twin, lb, noener), NBK_FN(ewald_twin, none, noener) } };

#undef NBK_FN

#endif
