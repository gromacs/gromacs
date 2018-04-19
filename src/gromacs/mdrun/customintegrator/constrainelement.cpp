//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/mdlib/update.h"
#include "gromacs/topology/topology.h"
#include "constrainelement.h"


void gmx::ConstrainElement::run()
{
    constrain_coordinates(dataManager.step, &dataManager.dvdl_constr,
                          integrator.inputrec, dataManager.mdatoms,
                          dataManager.state,
                          integrator.fr->bMolPBC,
                          &dataManager.top->idef, dataManager.shake_vir,
                          integrator.cr, integrator.ms,
                          integrator.nrnb, integrator.wcycle,
                          dataManager.upd, integrator.constr,
                          dataManager.bCalcVir, dataManager.do_log,
                          dataManager.do_ene);


    // TODO: Find a better way - needs to be handled by data manager
    finish_update(integrator.inputrec, dataManager.mdatoms,
                  dataManager.state, dataManager.graph,
                  integrator.nrnb, integrator.wcycle,
                  dataManager.upd, integrator.constr);
}
