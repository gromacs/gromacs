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
#ifndef _NB_KERNEL_IA64_SINGLE_H_
#define _NB_KERNEL_IA64_SINGLE_H_

/*! \file  nb_kernel_ia64_double.h
 *  \brief ia64 double precision assembly nonbonded kernels.
 *
 *  \internal
 */

#include <stdio.h>

#include <gmx_types.h>
#include <gmx_neighborlist.h>
#include <gmx_nonbonded.h>

/* Include kernel headers in local directory.
 * We can only have one routine in each file due to a bug
 * in the intel assembler program...
 */
#include "nb_kernel010_ia64_double.h"
#include "nb_kernel010nf_ia64_double.h"
#include "nb_kernel030_ia64_double.h"
#include "nb_kernel030nf_ia64_double.h"
#include "nb_kernel100_ia64_double.h"
#include "nb_kernel100nf_ia64_double.h"
#include "nb_kernel110_ia64_double.h"
#include "nb_kernel110nf_ia64_double.h"
#include "nb_kernel200_ia64_double.h"
#include "nb_kernel200nf_ia64_double.h"
#include "nb_kernel210_ia64_double.h"
#include "nb_kernel210nf_ia64_double.h"
#include "nb_kernel300_ia64_double.h"
#include "nb_kernel300nf_ia64_double.h"
#include "nb_kernel310_ia64_double.h"
#include "nb_kernel310nf_ia64_double.h"
#include "nb_kernel330_ia64_double.h"
#include "nb_kernel330nf_ia64_double.h"
#include "nb_kernel400_ia64_double.h"
#include "nb_kernel400nf_ia64_double.h"
#include "nb_kernel410_ia64_double.h"
#include "nb_kernel410nf_ia64_double.h"
#include "nb_kernel430_ia64_double.h"
#include "nb_kernel430nf_ia64_double.h"




#endif

