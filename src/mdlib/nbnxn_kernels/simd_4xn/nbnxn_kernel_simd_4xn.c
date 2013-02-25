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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

#ifdef GMX_NBNXN_SIMD_4XN

/* Include the half-width SIMD macros */
#ifdef GMX_NBNXN_HALF_WIDTH_SIMD
#define GMX_USE_HALF_WIDTH_SIMD_HERE
#endif
#include "gmx_simd_macros.h"
#include "gmx_simd_vec.h"

/* Declare the extern function this file defines */
#include "nbnxn_kernel_simd_4xn.h"

#if !(GMX_SIMD_WIDTH_HERE == 2 || GMX_SIMD_WIDTH_HERE == 4 || GMX_SIMD_WIDTH_HERE == 8)
#error "unsupported SIMD width"
#endif

/* Declare all the kernel functions */
#include "nbnxn_kernel_simd_4xn_rf_comb_geom_ener.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_geom_energrp.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_geom_noener.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_lb_ener.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_lb_energrp.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_lb_noener.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_none_ener.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_none_energrp.h"
#include "nbnxn_kernel_simd_4xn_rf_comb_none_noener.h"

#include "nbnxn_kernel_simd_4xn_tab_comb_geom_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_geom_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_geom_noener.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_lb_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_lb_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_lb_noener.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_none_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_none_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_comb_none_noener.h"

#include "nbnxn_kernel_simd_4xn_tab_twin_comb_geom_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_geom_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_geom_noener.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_lb_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_lb_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_lb_noener.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_none_ener.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_none_energrp.h"
#include "nbnxn_kernel_simd_4xn_tab_twin_comb_none_noener.h"

#include "nbnxn_kernel_simd_4xn_ewald_comb_geom_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_geom_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_geom_noener.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_lb_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_lb_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_lb_noener.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_none_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_none_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_comb_none_noener.h"

#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_geom_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_geom_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_geom_noener.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_lb_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_lb_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_lb_noener.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_none_ener.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_none_energrp.h"
#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_none_noener.h"

/* Define the macros to support constructing the kernel function
 * pointer lookup tables. */
#define NBK_FN(elec, ljcomb, enertype) nbnxn_kernel_simd_4xn ## _ ## elec ## _comb_ ## ljcomb ## _ ## enertype

/* Is there a better way to define/name this? Does it belong in
 * gmx_simd_macros.h? */
#define GMX_SIMD_J_UNROLL_SIZE 1

#else

#define GMX_NBNXN_SIMD_STUB_DECLARATION_ONLY

#endif /* GMX_NBNXN_SIMD_4XN */

/* Define the macro for the kernel function name, whether it will be a
 * stub or not. */
#define NBK_FN_ROOT nbnxn_kernel_simd_4xn

#include "../simd_common/nbnxn_kernel_simd_common.c"

#undef NBK_FN_ROOT
#undef NBK_FN
#undef GMX_SIMD_J_UNROLL_SIZE
#undef GMX_NBNXN_SIMD_STUB_DECLARATION_ONLY
