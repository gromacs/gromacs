/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief Declares the routine running the inetgrators.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_RUNNER_H
#define GMX_MDLIB_RUNNER_H

#include <cstdio>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;
struct t_commrec;
struct t_filenm;

namespace gmx
{

/*! \brief Driver routine, that calls the different methods
 *
 * \param[in] hw_opt   Hardware detection structure
 * \param[in] fplog    File pointer for log file
 * \param[in] cr       Communication data
 * \param[in] nfile    Number of files
 * \param[in] fnm      Array of filenames and file properties
 * \param[in] oenv     Output variables for storing xvg files etc.
 * \param[in] bVerbose Verbose output or not
 * \param[in] nstglobalcomm Number of steps between global communication
 * \param[in] ddxyz    Division of sub-boxes over processors for
 *                     use in domain decomposition parallellization
 * \param[in] dd_rank_order Ordering of the PP and PME ranks
 * \param[in] npme     The number of separate PME ranks requested, -1 = auto
 * \param[in] rdd      The maximum distance for bonded interactions with DD (nm)
 * \param[in] rconstr  Maximum distance for P-LINCS (nm)
 * \param[in] dddlb_opt File name for debugging
 * \param[in] dlb_scale File name for debugging
 * \param[in] ddcsx     File name for debugging
 * \param[in] ddcsy     File name for debugging
 * \param[in] ddcsz     File name for debugging
 * \param[in] nbpu_opt  Type of nonbonded processing unit
 * \param[in] nstlist_cmdline  Override neighbor search frequency
 * \param[in] nsteps_cmdline   Override number of simulation steps
 * \param[in] nstepout     How often to write to the console
 * \param[in] resetstep    Reset the step counter
 * \param[in] nmultisim    Number of parallel simulations to run
 * \param[in] repl_ex_nst  Number steps between replica exchange attempts
 * \param[in] repl_ex_nex  Number of replicas in REMD
 * \param[in] repl_ex_seed The seed for Monte Carlo swaps
 * \param[in] pforce       Minimum force for printing (for debugging)
 * \param[in] cpt_period    How often to checkpoint the simulation
 * \param[in] max_hours     Maximume length of the simulation (wall time)
 * \param[in] imdport       Interactive MD port (socket)
 * \param[in] Flags         More command line options
 */
int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, struct t_commrec *cr, int nfile,
             const t_filenm fnm[], const gmx_output_env_t *oenv, gmx_bool bVerbose,
             int nstglobalcomm, ivec ddxyz, int dd_rank_order, int npme,
             real rdd, real rconstr, const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt, int nstlist_cmdline,
             gmx_int64_t nsteps_cmdline, int nstepout, int resetstep,
             int nmultisim, int repl_ex_nst, int repl_ex_nex,
             int repl_ex_seed, real pforce, real cpt_period, real max_hours,
             int imdport, unsigned long Flags);


}      // namespace gmx

#endif // GMX_MDLIB_RUNNER_H
