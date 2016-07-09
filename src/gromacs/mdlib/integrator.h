/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
/*! \brief Declares the integrator type for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_INTEGRATOR_H
#define GMX_MDLIB_INTEGRATOR_H

#include <cstdio>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_output_env_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;

namespace gmx
{
/*! \brief Integrator algorithm implementation.
 *
 * \param[in] fplog               Log file for output
 * \param[in] cr                  Communication record
 * \param[in] nfile               Number of files
 * \param[in] fnm                 Filename structure array
 * \param[in] oenv                Output information
 * \param[in] bVerbose            Verbose output or not
 * \param[in] nstglobalcomm       How often global communication is done
 * \param[in] vsite               Virtual site information
 * \param[in] constr              Constraint information
 * \param[in] stepout             How often we writen to the console
 * \param[in] inputrec            Input record with mdp options
 * \param[in] top_global          Molecular topology for the whole system
 * \param[in] fcd                 Force and constraint data
 * \param[in] state_global        The state (x, v, f, box etc.) of the whole system
 * \param[in] mdatoms             Structure containing atom information
 * \param[in] nrnb                Accounting for floating point operations
 * \param[in] wcycle              Wall cycle timing information
 * \param[in] ed                  Essential dynamics sampling information
 * \param[in] fr                  Force record with cut-off information and more
 * \param[in] repl_ex_nst         How often we do replica exchange (in steps)
 * \param[in] repl_ex_nex         How many replicas we have
 * \param[in] repl_ex_seed        The seed for Monte Carlo swaps
 * \param[in] membed              Membrane embedding data structure
 * \param[in] cpt_period          How often to checkpoint the simulation
 * \param[in] max_hours           Maximume length of the simulation (wall time)
 * \param[in] imdport             Interactive MD port (socket)
 * \param[in] Flags               Flags to control mdrun
 * \param[in] walltime_accounting More timing information
 */
typedef double integrator_t (FILE *fplog, t_commrec *cr,
                             int nfile, const t_filenm fnm[],
                             const gmx_output_env_t *oenv, gmx_bool bVerbose,
                             int nstglobalcomm,
                             gmx_vsite_t *vsite, gmx_constr_t constr,
                             int stepout,
                             t_inputrec *inputrec,
                             gmx_mtop_t *top_global, t_fcdata *fcd,
                             t_state *state_global,
                             t_mdatoms *mdatoms,
                             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                             gmx_edsam_t ed,
                             t_forcerec *fr,
                             int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                             gmx_membed_t gmx_unused * membed,
                             real cpt_period, real max_hours,
                             int imdport,
                             unsigned long Flags,
                             gmx_walltime_accounting_t walltime_accounting);

}      // namespace gmx

#endif // GMX_MDLIB_INTEGRATOR_H
