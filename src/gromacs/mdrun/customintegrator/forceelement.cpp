//
// Created by pascal on 4/20/18.
//

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "forceelement.h"

gmx::ForceElement::ForceElement(Integrator &integ, DataManager &dm) : Element(integ, dm)
{
    doPCoupling   = integrator.inputrec->epc != epcNO;
    nstPCouple    = integrator.inputrec->nstpcouple;
    nstCalcEnergy = integrator.inputrec->nstcalcenergy;

    if (DOMAINDECOMP(integrator.cr))
    {
        ddOpenBalanceRegion  = DdOpenBalanceRegionBeforeForceComputation::yes;
        ddCloseBalanceRegion = DdCloseBalanceRegionAfterForceComputation::yes;
    }
    else
    {
        ddOpenBalanceRegion  = DdOpenBalanceRegionBeforeForceComputation::no;
        ddCloseBalanceRegion = DdCloseBalanceRegionAfterForceComputation::no;
    }
}

void gmx::ForceElement::run()
{
    bool bCalcEner = dataManager.bLastStep || do_per_step(dataManager.step, nstCalcEnergy);
    bool bCalcVir  = bCalcEner ||
        (doPCoupling && do_per_step(dataManager.step, nstPCouple));

    // TODO: implement FEP
    bool bDoFEP = false;

    int  force_flags = ((dataManager.bNS ? GMX_FORCE_NS : 0) |
                        GMX_FORCE_STATECHANGED |
                        ((inputrecDynamicBox(integrator.inputrec)) ? GMX_FORCE_DYNAMICBOX : 0) |
                        GMX_FORCE_ALLFORCES |
                        (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                        (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                        (bDoFEP ? GMX_FORCE_DHDL : 0)
                        );

    do_force(integrator.fplog, integrator.cr, integrator.ms,
             integrator.inputrec, dataManager.step, integrator.nrnb,
             integrator.wcycle, dataManager.top, dataManager.groups,
             dataManager.state->box, dataManager.state->x,
             &dataManager.state->hist,
             dataManager.f, dataManager.force_vir,
             dataManager.mdatoms, dataManager.enerd,
             integrator.fcd,
             dataManager.state->lambda, dataManager.graph,
             integrator.fr, integrator.vsite,
             dataManager.mu_tot, dataManager.t, dataManager.ed,
             force_flags,
             ddOpenBalanceRegion, ddCloseBalanceRegion);
}
