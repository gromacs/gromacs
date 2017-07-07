/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <memory>

#include "customMD.h"

#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/inputrec.h"

#include "statemanager.h"
#include "integratorFactory.h"
#include "sequenceFactory.h"
#include "elements/elementBase.h"
#include "elements/integratorElement.h"

/*! \libinternal
    \copydoc integrator_t (FILE *fplog, t_commrec *cr, const gmx::MDLogger &mdlog,
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
                           real cpt_period, real max_hours,
                           int imdport,
                           unsigned long Flags,
                           gmx_walltime_accounting_t walltime_accounting)
 */
double gmx::do_customMD(FILE *fplog, t_commrec *cr, const gmx::MDLogger &mdlog,
                  int nfile, const t_filenm fnm[],
                  const gmx_output_env_t *oenv, gmx_bool bVerbose,
                  int nstglobalcomm,
                  gmx_vsite_t *vsite, gmx_constr_t constr,
                  int stepout, t_inputrec *ir,
                  gmx_mtop_t *top_global,
                  t_fcdata *fcd,
                  t_state *state_global,
                  energyhistory_t *energyHistory,
                  t_mdatoms *mdatoms,
                  t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                  gmx_edsam_t ed, t_forcerec *fr,
                  int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                  gmx_membed_t *membed,
                  real cpt_period, real max_hours,
                  int imdport,
                  unsigned long Flags,
                  gmx_walltime_accounting_t walltime_accounting)
{
    /* 
        A boolean that stores information whether or not a constraint
        will be applied. The constraint status cannot be obtained from
        only inputrec data structure, or the .mdp file, as constraints
        can also be defined in topology file (ex/ water with SETTLE).
        bConstrain is therefore needed in addition to the inputrec to
        pass information about the status of constraints to the sequence.
    */
    bool bConstrain = FALSE;
    
    if (constr)
    {
        bConstrain = TRUE;
    }
    
    /* 
        Create the sequences, which are lists of pointers to the element objects.
        The sequence factory creates a sequence using the information from mdp
        (ir) file and constraint status (bConstrain).
        There are three types of sequences:
        
        * preRunSequence : Is called once before the MD loop 
          (for elements that needs to be run only once before the md loop starts)
        * runSequence : Is called inside the integrator loop
          (this sequence is the sequence of elements that make up the integrator)
        * postRunSequence : Is called once after the MD loop
          (for elements that needs to be run only once, after conclusions of the md loop)
    */
    std::list<Element*> preRunSequence =
            SequenceFactory::createPreRunSequence(ir, bConstrain);
    std::list<Element*> runSequence =
            SequenceFactory::createRunSequence(ir, bConstrain);
    std::list<Element*> postRunSequence =
            SequenceFactory::createPostRunSequence(ir, bConstrain);
    
    /*
        Using the constructed preRunSequence and runSequence, the Integrator
        factory will create an integrator object.
    */
    std::unique_ptr<Integrator> integrator =
            IntegratorFactory::create(runSequence,
                                      preRunSequence,
                                      postRunSequence);
    
    /*
        The StateManager object holds all the data that will be needed for the
        integrator.

        *** !!! Temporary Solution: !!!***
        Currently this object holds all the data that the standard do_md()
        function in the GROMACS code holds. All elements in the integrator
        sequences (preRunSequence and runSequence) have full read and write
        access to all the data.

        A general discussion is needed to decide on the data structures,
        categorization and data access levels for the elements. Please see the
        accompanying document, section IV (Data reorganization) for more details.
    */
    StateManager state_manager(fplog, cr, mdlog, nfile, fnm,
                               oenv, bVerbose,
                               nstglobalcomm,
                               vsite, constr,
                               stepout, ir, top_global,
                               fcd, state_global, energyHistory,
                               mdatoms, nrnb, wcycle, ed, fr,
                               repl_ex_nst, repl_ex_nex, repl_ex_seed,
                               membed,
                               cpt_period, max_hours,
                               imdport,
                               Flags,
                               walltime_accounting);
    
    /*
        Initialization of the data in state_manager.

        *Final goal:* Should initialize _only_ data which is not dependent on
        the chosen integrator algorithm (e.g. initialize domain decomposition,
        output files, etc.)

        *Current status:* Is still performing some integrator-dependent
        initializations, which should be moved to the integrator->prerun() call.
    */
    state_manager.loopSetup();

    /*
        The initialize method gives data access (by giving a reference to the
        state_manager object) to all the elements in the integrator sequences.
    */
    integrator->initialize(state_manager);
    
    /*
        The prerun() method loops through the elements of the preRunSequence,
        calling their respective run() methods in order.
        This initializes any integrator-dependent variables that might need
        initialization (e.g. algorithms that need a force calculation before the
        start of the first integrator step)
     */
    integrator->prerun();

    /*
        The actual integrator loop, executed `nsteps` times.
     */
    while (!state_manager.bLastStep)
    {
        /*
            stepSetup() performs any changes in integrator-independent
            variables of the state manager occurring before every integrator
            step.
        */
        state_manager.stepSetup();

        /*
            The run() method loops through the elements of the runSequence,
            calling their respective run() methods in order.
            This call represents the actual integrator.
         */
        integrator->run();
        
        /*
            stepTeardown() performs any changes in integrator-independent
            variables of the state manager occurring after every integrator
            step.
        */
        state_manager.stepTeardown();
    }

    /*
        The postrun() method loops through the elements of the postRunSequence,
        calling their respective run() methods in order.
        This finalizes any integrator-dependent variables that might need
        finalization
     */
    integrator->postrun();

    /*
        Finalization of the data in state_manager.

        *Final goal:* Should finalize _only_ data which is not dependent on
        the chosen integrator algorithm.

        *Current status:* Is still performing some integrator-dependent
        finalizations, which should be moved to the integrator->postrun() call.
    */
    state_manager.loopTeardown();
    
    return 0;
}
