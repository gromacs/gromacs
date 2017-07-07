#include "updateVelocitiesElement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/update.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void UpdateVelocities::initialize(StateManager& dataRef)
{
    _data = &dataRef;
    dt = step_fraction*_data->ir->delta_t;
}

void UpdateVelocities::run()
{
    /*
     * v(t+dt) = v(t) + (dt/m)*f(t)
     * 
    */

    set_vprime(_data->mdatoms, _data->upd, _data->state);
    set_fprime(_data->mdatoms, _data->upd, &_data->f);
    progress_vprime(_data->mdatoms, _data->ir, dt, _data->upd);
    set_v(_data->mdatoms, _data->upd, _data->state);
}
