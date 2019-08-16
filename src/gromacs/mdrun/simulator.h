/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Declares the simulator interface for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_SIMULATOR_H
#define GMX_MDRUN_SIMULATOR_H

#include <cstdio>

#include <memory>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class energyhistory_t;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct ObservablesHistory;
struct pull_t;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_filenm;
struct t_inputrec;
struct t_nrnb;
struct t_swap;
class t_state;

namespace gmx
{

enum class StartingBehavior;
class BoxDeformation;
class Constraints;
class PpForceWorkload;
class IMDOutputProvider;
class ImdSession;
class MDLogger;
class MDAtoms;
class StopHandlerBuilder;
struct MdrunOptions;

//! Function type for simulator code.
using SimulatorFunctionType = void();

/*! \internal
 * \brief Struct to handle setting up and running the different simulation types.
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
 * Having multiple simulation types as member functions isn't a good
 * design, and we definitely only intend one to be called, but the
 * goal is to make it easy to change the names and types of members
 * without having to make identical changes in several places in the
 * code. Once many of them have become modules, we should change this
 * approach.
 *
 * Use a braced initializer list to construct one of these. */
struct Simulator
{
    //! Handles logging.
    FILE                               *fplog;
    //! Handles communication.
    t_commrec                          *cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t               *ms;
    //! Handles logging.
    const MDLogger                     &mdlog;
    //! Count of input file options.
    int                                 nfile;
    //! Content of input file options.
    const t_filenm                     *fnm;
    //! Handles writing text output.
    const gmx_output_env_t             *oenv;
    //! Contains command-line options to mdrun.
    const MdrunOptions                 &mdrunOptions;
    //! Whether the simulation will start afresh, or restart with/without appending.
    StartingBehavior                    startingBehavior;
    //! Handles virtual sites.
    gmx_vsite_t                        *vsite;
    //! Handles constraints.
    Constraints                        *constr;
    //! Handles enforced rotation.
    gmx_enfrot                         *enforcedRotation;
    //! Handles box deformation.
    BoxDeformation                     *deform;
    //! Handles writing output files.
    IMDOutputProvider                  *outputProvider;
    //! Contains user input mdp options.
    t_inputrec                         *inputrec;
    //! The Interactive Molecular Dynamics session.
    ImdSession                         *imdSession;
    //! The pull work object.
    pull_t                             *pull_work;
    //! The coordinate-swapping session.
    t_swap                             *swap;
    //! Full system topology.
    gmx_mtop_t                         *top_global;
    //! Helper struct for force calculations.
    t_fcdata                           *fcd;
    //! Full simulation state (only non-nullptr on master rank).
    t_state                            *state_global;
    //! History of simulation observables.
    ObservablesHistory                 *observablesHistory;
    //! Atom parameters for this domain.
    MDAtoms                            *mdAtoms;
    //! Manages flop accounting.
    t_nrnb                             *nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle                      *wcycle;
    //! Parameters for force calculations.
    t_forcerec                         *fr;
    //! Data for energy output.
    gmx_enerdata_t                     *enerd;
    //! Schedule of force-calculation work each step for this task.
    PpForceWorkload                    *ppForceWorkload;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters    &replExParams;
    //! Parameters for membrane embedding.
    gmx_membed_t                       *membed;
    //! Manages wall time accounting.
    gmx_walltime_accounting            *walltime_accounting;
    //! Registers stop conditions
    std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder;
    //! Implements the normal MD simulations.
    SimulatorFunctionType               do_md;
    //! Implements the rerun functionality.
    SimulatorFunctionType               do_rerun;
    //! Implements steepest descent EM.
    SimulatorFunctionType               do_steep;
    //! Implements conjugate gradient energy minimization
    SimulatorFunctionType               do_cg;
    //! Implements onjugate gradient energy minimization using the L-BFGS algorithm
    SimulatorFunctionType               do_lbfgs;
    //! Implements normal mode analysis
    SimulatorFunctionType               do_nm;
    //! Implements test particle insertion
    SimulatorFunctionType               do_tpi;
    //! Implements MiMiC QM/MM workflow
    SimulatorFunctionType               do_mimic;
    /*! \brief Function to run the correct SimulatorFunctionType,
     * based on the .mdp integrator field. */
    void run(unsigned int ei, bool doRerun);
};

}      // namespace gmx

#endif // GMX_MDRUN_SIMULATOR_H
