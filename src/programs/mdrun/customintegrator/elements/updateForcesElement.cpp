#include "updateForcesElement.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"

#include "programs/mdrun/customintegrator/statemanager.h"

void UpdateForces::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void UpdateForces::run()
{
    clear_mat(_data->force_vir);
    do_force(_data->fplog, _data->cr, _data->ir, _data->step, _data->nrnb, _data->wcycle, _data->top, _data->groups,
             _data->state->box, &_data->state->x, &_data->state->hist,
             &_data->f, _data->force_vir, _data->mdatoms, _data->enerd, _data->fcd,
             &_data->state->lambda, _data->graph,
             _data->fr, _data->vsite, _data->mu_tot, _data->t, _data->ed, _data->bBornRadii,
             (_data->bNS ? GMX_FORCE_NS : 0) | _data->force_flags);
}