//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/domdec/domdec.h"
#include "computeglobalselement.h"

void gmx::ComputeGlobalsElement::run()
{
    /* ############## IF NOT VV, Calculate globals HERE  ############ */
    /* With Leap-Frog we can skip compute_globals at
     * non-communication steps, but we need to calculate
     * the kinetic energy one step before communication.
     */
    // Organize to do inter-simulation signalling on steps if
    // and when algorithms require it.
    bool doInterSimSignal = false;

    if (dataManager.bGStat || (do_per_step(dataManager.step+1, dataManager.nstglobalcomm)) || doInterSimSignal)
    {
        // Since we're already communicating at this step, we
        // can propagate intra-simulation signals. Note that
        // check_nstglobalcomm has the responsibility for
        // choosing the value of nstglobalcomm that is one way
        // bGStat becomes true, so we can't get into a
        // situation where e.g. checkpointing can't be
        // signalled.
        bool                doIntraSimSignal = true;
        SimulationSignaller signaller(&dataManager.signals, integrator.cr, integrator.ms, doInterSimSignal, doIntraSimSignal);

        compute_globals(integrator.fplog, dataManager.gstat, integrator.cr, integrator.inputrec, integrator.fr, dataManager.ekind, dataManager.state, dataManager.mdatoms, integrator.nrnb, dataManager.vcm,
                        integrator.wcycle, dataManager.enerd, dataManager.force_vir, dataManager.shake_vir, dataManager.total_vir, dataManager.pres, dataManager.mu_tot,
                        integrator.constr, &signaller,
                        dataManager.lastbox,
                        &dataManager.totalNumberOfBondedInteractions, &dataManager.bSumEkinhOld,
                        (dataManager.bGStat ? CGLO_GSTAT : 0)
                        | CGLO_ENERGY
                        | (dataManager.bStopCM ? CGLO_STOPCM : 0)
                        | CGLO_TEMPERATURE
                        | CGLO_PRESSURE
                        | CGLO_CONSTRAINT
                        | (dataManager.shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0)
                        );
        checkNumberOfBondedInteractions(integrator.fplog, integrator.cr, dataManager.totalNumberOfBondedInteractions,
                                        integrator.top_global, dataManager.top, dataManager.state,
                                        &dataManager.shouldCheckNumberOfBondedInteractions);
    }

    if (!dataManager.bGStat)
    {
        /* We will not sum ekinh_old,
         * so signal that we still have to do it.
         */
        dataManager.bSumEkinhOld = TRUE;
    }
}
