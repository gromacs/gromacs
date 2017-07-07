#include "updatePositionsElement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/update.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void UpdatePositions::initialize(StateManager& dataRef)
{
    _data = &dataRef;
    dt = step_fraction*_data->ir->delta_t;
}

void UpdatePositions::run()
{
    /*
     * x(t+dt) = x(t) + v(t)*dt
     * 
    */

    /* 
     * Note for copy_mat call
     * lastbox is used by computeglobals element, 
     * the computeglobal operates on the boxvectors to calculate the pressure.
     * and is set before the position update to the state->box variable
     * This is done incase the update position changes the box vector.
     * Ex/ When doing pressure coupling the update position will change the box vector
     * 
    */
    
    copy_mat(_data->state->box, _data->lastbox);
    
    save_x(_data->mdatoms, _data->upd, _data->state);
    set_vprime(_data->mdatoms, _data->upd, _data->state);
    set_xprime(_data->mdatoms, _data->upd, _data->state);
    progress_xprime(_data->mdatoms, _data->ir, dt, _data->upd);
    set_x(_data->mdatoms, _data->upd, _data->state);
}
