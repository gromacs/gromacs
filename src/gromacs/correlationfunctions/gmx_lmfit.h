/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
/*! \libinternal
 * \file
 * \brief
 * Declares the interface to an lmfit library
 *
 * We might link to either an external or internal lmfit library, and
 * the latter case we choose to rename the symbols exported by that
 * code, so that code linking to GROMACS will not have symbol clashes
 * with its own possible lmfit dependency. This header file hides
 * those details from the GROMACS source files that depend on lmfit
 * functionality. The build system needs to define HAVE_EXTERNAL_LMFIT
 * to 0/1 appropriately, and set the system include path as needed.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_CORRELATION_FUNCTIONS_GMX_LMFIT_H
#define GMX_CORRELATION_FUNCTIONS_GMX_LMFIT_H

#if HAVE_EXTERNAL_LMFIT

#include <lmmin.h>
#include <lmstruct.h>

#define gmx_lm_control_double lm_control_double
#define gmx_lm_control_float lm_control_float
#define gmx_lm_infmsg lm_infmsg
#define gmx_lmmin lmmin

#else

#include "external/lmfit/gmx_lmmin.h"
#include "external/lmfit/gmx_lmstruct.h"

#endif

#endif
