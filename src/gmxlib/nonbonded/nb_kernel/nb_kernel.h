/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
 *
 * $Id$
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _NB_KERNEL_GENERIC_H_
#define _NB_KERNEL_GENERIC_H_

/** Portable standard standard nonbonded kernels.
 * 
 *  \file  nb_kernel.h
 *
 *  \internal
 */

#include <stdio.h>

#include "../nb_kerneltype.h"

void
nb_kernel_setup(FILE *fplog,nb_kernel_t **list);

#endif

