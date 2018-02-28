/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018, by the GROMACS development team, led by
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
/*! \internal
 * \brief Declares the integrator interface for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_INTEGRATOR_H
#define GMX_MDRUN_INTEGRATOR_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class energyhistory_t;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct MdrunOptions;
struct ObservablesHistory;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_filenm;
struct t_inputrec;
struct t_nrnb;
class t_state;

namespace gmx
{

class Constraints;
class IMDOutputProvider;
class MDLogger;
class MDAtoms;

//! Function type for integrator code.
using IntegratorFunctionType = void();

/*! \internal
 * \brief Struct to handle setting up and running the different "integrators".
 *
 * This struct is a mere aggregate of parameters to pass to evaluate an
 * energy, so that future changes to names and types of them consume
 * less time when refactoring other code.
 *
 * Aggregate initialization is used, for which the chief risk is that
 * if a member is added at the end and not all initializer lists are
 * updated, then the member will be value initialized, which will
 * typically mean initialization to zero.
 *
 * Having multiple integrators as member functions isn't a good
 * design, and we definitely only intend one to be called, but the
 * goal is to make it easy to change the names and types of members
 * without having to make identical changes in several places in the
 * code. Once many of them have become modules, we should change this
 * approach.
 *
 * Note that the presence of const reference members means that the
 * default constructor would be implicitly deleted. But we only want
 * to make one of these when we know how to initialize these members,
 * so that is perfect. To ensure this remains true even if we would
 * remove those members, we explicitly delete this constructor.
 * Other constructors, copies and moves are OK. */
struct Integrator
{
    //! Handles logging.
    FILE                            *fplog;
    //! Handles communication.
    t_commrec                       *cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t            *ms;
    //! Handles logging.
    const MDLogger                  &mdlog;
    //! Count of input file options.
    int                              nfile;
    //! Content of input file options.
    const t_filenm                  *fnm;
    //! Handles writing text output.
    const gmx_output_env_t          *oenv;
    //! Contains command-line options to mdrun.
    const MdrunOptions              &mdrunOptions;
    //! Handles virtual sites.
    gmx_vsite_t                     *vsite;
    //! Handles constraints.
    Constraints                     *constr;
    //! Handles writing output files.
    IMDOutputProvider               *outputProvider;
    //! Contains user input mdp options.
    t_inputrec                      *inputrec;
    //! Full system topology.
    gmx_mtop_t                      *top_global;
    //! Helper struct for force calculations.
    t_fcdata                        *fcd;
    //! Full simulation state (only non-nullptr on master rank).
    t_state                         *state_global;
    //! History of simulation observables.
    ObservablesHistory              *observablesHistory;
    //! Atom parameters for this domain.
    MDAtoms                         *mdAtoms;
    //! Manages flop accounting.
    t_nrnb                          *nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle                   *wcycle;
    //! Parameters for force calculations.
    t_forcerec                      *fr;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters &replExParams;
    //! Parameters for membrane embedding.
    gmx_membed_t                    *membed;
    //! Manages wall time accounting.
    gmx_walltime_accounting         *walltime_accounting;
    //! Implements the normal MD integrators.
    IntegratorFunctionType           do_md;
    //! Implements steepest descent EM.
    IntegratorFunctionType           do_steep;
    //! Implements conjugate gradient energy minimization
    IntegratorFunctionType           do_cg;
    //! Implements onjugate gradient energy minimization using the L-BFGS algorithm
    IntegratorFunctionType           do_lbfgs;
    //! Implements normal mode analysis
    IntegratorFunctionType           do_nm;
    //! Implements test particle insertion
    IntegratorFunctionType           do_tpi;
    /*! \brief Function to run the correct IntegratorFunctionType,
     * based on the .mdp integrator field. */
    void run(unsigned int ei);
    //! We only intend to construct such objects with an initializer list.
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9)
    // Aspects of the C++11 spec changed after GCC 4.8.5, and
    // compilation of the initializer list construction in runner.cpp
    // fails in GCC 4.8.5.
    Integrator() = delete;
#endif
};

}      // namespace gmx

#endif // GMX_MDRUN_INTEGRATOR_H
