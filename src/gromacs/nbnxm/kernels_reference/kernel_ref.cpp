/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#include "gmxpre.h"

#include "kernel_ref.h"

#include <cassert>
#include <cmath>

#include <algorithm>

#include "../nbnxm_geometry.h"

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/simd_energy_accumulator.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* Analytical reaction-field kernels */
#define CALC_COUL_RF
#define LJ_CUT
#include "kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
#undef CALC_COUL_RF


/* Tabulated exclusion interaction electrostatics kernels */
#define CALC_COUL_TAB
#define LJ_CUT
#include "kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
/* Twin-range cut-off kernels */
#define VDW_CUTOFF_CHECK
#define LJ_CUT
#include "kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
#undef VDW_CUTOFF_CHECK
#undef CALC_COUL_TAB
