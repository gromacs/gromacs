//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/state.h"
#include "updatefusedelement.h"

void gmx::UpdateFusedElement::run()
{
    /* Box is changed in update() when we do pressure coupling,
     * but we should still use the old box for energy corrections and when
     * writing it to the energy file, so it matches the trajectory files for
     * the same timestep above. Make a copy in a separate array.
     */
    // TODO: Find a better way - needs to be handled by data manager
    copy_mat(dataManager.state->box, dataManager.lastbox);

    wallcycle_start(integrator.wcycle, ewcUPDATE);
    update_coords(dataManager.step, integrator.inputrec, dataManager.mdatoms, dataManager.state, dataManager.f, integrator.fcd,
                  dataManager.ekind, dataManager.M, dataManager.upd, etrtPOSITION, integrator.cr, integrator.constr);
    wallcycle_stop(integrator.wcycle, ewcUPDATE);

    // TODO: Find a better way - needs to be handled by data manager
    if (!integrator.constr)
    {
        finish_update(integrator.inputrec, dataManager.mdatoms,
                      dataManager.state, dataManager.graph,
                      integrator.nrnb, integrator.wcycle, dataManager.upd, integrator.constr);
    }
}
