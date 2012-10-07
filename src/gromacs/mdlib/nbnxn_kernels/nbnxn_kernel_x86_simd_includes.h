/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */

/* This files includes all x86 SIMD kernel flavors.
 * Only the Electrostatics type and optionally the VdW cut-off check
 * need to be set before including this file.
 */

/* Include the force+energy kernels */
#define CALC_ENERGIES
#define LJ_COMB_GEOM
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
#undef CALC_ENERGIES

/* Include the force+energygroups kernels */
#define CALC_ENERGIES
#define ENERGY_GROUPS
#define LJ_COMB_GEOM
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
#undef ENERGY_GROUPS
#undef CALC_ENERGIES

/* Include the force only kernels */
#define LJ_COMB_GEOM
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_GEOM
#define LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
#undef LJ_COMB_LB
#include "nbnxn_kernel_x86_simd_outer.h"
