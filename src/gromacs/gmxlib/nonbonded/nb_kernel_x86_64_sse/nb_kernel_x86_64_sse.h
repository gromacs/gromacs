/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
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
#ifndef _NB_KERNEL_X86_64_SSE_H_
#define _NB_KERNEL_X86_64_SSE_H_

/*! \file  nb_kernel_x86_64_sse.h
 *  \brief x86_64 SSE-optimized level2 nonbonded kernels.
 *
 *  \internal
 */

#include <stdio.h>

#include <types/simple.h>

#include "../nb_kerneltype.h"
#include "nb_kernel_allvsall_sse2_single.h"
#include "nb_kernel_allvsallgb_sse2_single.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void
nb_kernel_setup_x86_64_sse(FILE *log,nb_kernel_t **list);

#ifdef __cplusplus
}
#endif

#endif /* _NB_KERNEL_X86_64_SSE_H_ */
