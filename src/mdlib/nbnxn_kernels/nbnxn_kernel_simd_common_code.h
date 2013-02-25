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

/*! \brief Kinds of electrostatic treatments in SIMD Verlet kernels
 */
enum {
    coultRF, coultTAB, coultTAB_TWIN, coultEWALD, coultEWALD_TWIN, coultNR
};

#if   defined(GMX_NBNXN_SIMD_2XNN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _ener
#elif defined(GMX_NBNXN_SIMD_4XN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_4xn_ ## elec ## _comb_ ## ljcomb ## _ener
#else
#error "Unknown NBNXN SIMD type"
#endif

static p_nbk_func_ener p_nbk_ener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN

#if   defined(GMX_NBNXN_SIMD_2XNN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _energrp
#elif defined(GMX_NBNXN_SIMD_4XN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_4xn_ ## elec ## _comb_ ## ljcomb ## _energrp
#else
#error "Unknown NBNXN SIMD type"
#endif

static p_nbk_func_ener p_nbk_energrp[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN

#if   defined(GMX_NBNXN_SIMD_2XNN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_2xnn_ ## elec ## _comb_ ## ljcomb ## _noener
#elif defined(GMX_NBNXN_SIMD_4XN)
#define NBK_FN(elec, ljcomb) nbnxn_kernel_simd_4xn_ ## elec ## _comb_ ## ljcomb ## _noener
#else
#error "Unknown NBNXN SIMD type"
#endif

static p_nbk_func_noener p_nbk_noener[coultNR][ljcrNR] =
{ { NBK_FN(rf, geom), NBK_FN(rf, lb), NBK_FN(rf, none) },
  { NBK_FN(tab, geom), NBK_FN(tab, lb), NBK_FN(tab, none) },
  { NBK_FN(tab_twin, geom), NBK_FN(tab_twin, lb), NBK_FN(tab_twin, none) },
  { NBK_FN(ewald, geom), NBK_FN(ewald, lb), NBK_FN(ewald, none) },
  { NBK_FN(ewald_twin, geom), NBK_FN(ewald_twin, lb), NBK_FN(ewald_twin, none) } };
#undef NBK_FN

#endif
