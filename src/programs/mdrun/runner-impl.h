/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Declares the private implementation for running the inetgrators.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_RUNNER_IMPL_H
#define GMX_MDLIB_RUNNER_IMPL_H

#include <cstdio>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_filenm;

namespace gmx
{
// Declare overload of the mdrunner() defined in the header file. Documented
// at the definition below.
/*! \cond internal
 * \brief Declare implementation of the runner.
 */
int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
             int nstglobalcomm,
             ivec ddxyz, int dd_rank_order, int npme, real rdd, real rconstr,
             const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt, int nstlist_cmdline,
             gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
             int gmx_unused nmultisim,
             const ReplicaExchangeParameters &replExParams,
             real pforce, real cpt_period, real max_hours,
             int imdport, unsigned long Flags);
/// \endcond

}      // namespace gmx

#endif // GMX_MDLIB_RUNNER_IMPL_H
