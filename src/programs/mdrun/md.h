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
 * \brief Declares the integrators for molecular dynamics simulations
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_MD_H
#define GMX_MDLIB_MD_H

#include "gromacs/mdlib/integrator.h"

namespace gmx
{

//! MD simulations
integrator_t do_md;

// I think we can end up with a consolidated factory method at some point, but
// in the short term there is other reformulation of the simulation loop taking
// place, and ultimately I would like to see this functionality in libgromacs
// itself, so for now I'm just trying to encapsulate do_md a bit. I would also
// like to move to encapsulated arguments, such as struct mdrunner_arglist...
/*! \brief Implementation of integrator factory.
 *
 * Initial implementation can only produce CliIntegrators.
 */
class IntegratorFactoryImpl : public IIntegratorFactory
{
    /*! Implement IIntegratorFactory interface */
    std::shared_ptr<ICliIntegrator> get_cli_integrator(FILE *fplog, t_commrec *cr, const gmx::MDLogger &mdlog,
        int nfile, const t_filenm fnm[],
        const gmx_output_env_t *oenv, gmx_bool bVerbose,
        int nstglobalcomm,
        gmx_vsite_t *vsite, gmx_constr_t constr,
        int stepout, gmx::IMDOutputProvider *outputProvider,
        t_inputrec *ir,
        gmx_mtop_t *top_global,
        t_fcdata *fcd,
        t_state *state_global,
        ObservablesHistory *observablesHistory,
        t_mdatoms *mdatoms,
        t_nrnb *nrnb, gmx_wallcycle_t wcycle,
        gmx_edsam_t ed, t_forcerec *fr,
        int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
        gmx_membed_t *membed,
        real cpt_period, real max_hours,
        int imdport,
        unsigned long Flags,
        gmx_walltime_accounting_t walltime_accounting) override;
};

}      // namespace gmx

#endif // GMX_MDLIB_MD_H
