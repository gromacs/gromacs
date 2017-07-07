#include "constrainElement.h"

#include "gromacs/mdlib/update.h"
#include "gromacs/utility/basedefinitions.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void InitialConstrain::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void InitialConstrain::run()
{
   /* constrain the current position */
    // set xp to x, to constrain the current position x
    set_vprime(_data->mdatoms, _data->upd, _data->state);
    set_xprime(_data->mdatoms, _data->upd, _data->state);
    constrain_positions(_data->step,
                       &_data->dvdl_constr,
                       _data->ir,
                       _data->mdatoms,
                       _data->state,
                       _data->upd,
                       _data->fr->bMolPBC,
                       &_data->top->idef,
                       _data->shake_vir,
                       _data->cr,
                       _data->nrnb,
                       _data->constr,
                       FALSE,
                       0,
                       nullptr);
    set_x(_data->mdatoms, _data->upd, _data->state);

    if (halfStep)
    {
        real scl = -1.0;
        // Constrain halfStep velocity
        set_xprime(_data->mdatoms, _data->upd, _data->state);
        // scale velocities by -1
        scale_vprime(_data->mdatoms, _data->ir, _data->upd, scl);
        // progress positions to get x(-dt) 
        progress_xprime(_data->mdatoms, _data->ir, _data->ir->delta_t, _data->upd);
        // Apply constrain x(-t)= xp using x(0)=x as reference
        set_v(_data->mdatoms, _data->upd, _data->state);
        constrain_positions(_data->step,
                           &_data->dvdl_constr,
                           _data->ir,
                           _data->mdatoms,
                           _data->state,
                           _data->upd,
                           _data->fr->bMolPBC,
                           &_data->top->idef,
                           _data->shake_vir,
                           _data->cr,
                           _data->nrnb,
                           _data->constr,
                           FALSE, -1,
                           as_rvec_array(_data->state->v.data()));
        
        // xp = x(-t) xavg = x(0)
        set_vprime(_data->mdatoms, _data->upd, _data->state);
        //scale velocities by -1
        scale_vprime(_data->mdatoms, _data->ir, _data->upd, scl);
    }
    else
    {
        /* constrain the fullstep inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        constrain_velocities(_data->step,
                             &_data->dvdl_constr,
                             _data->ir,
                             _data->mdatoms,
                             _data->state,
                             _data->upd,
                             _data->fr->bMolPBC,
                             &_data->top->idef,
                             _data->shake_vir,
                             _data->cr,
                             _data->nrnb,
                             _data->constr,
                             _data->bCalcVir, 0);

    }

    set_v(_data->mdatoms, _data->upd, _data->state);
}

void ConstrainPositionAndVelocity::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void ConstrainPositionAndVelocity::run()
{
    load_x(_data->mdatoms, _data->upd, _data->state);
    // Uses x= x(t) and xp = x(t+dt)
    constrain_positions(_data->step,
                       &_data->dvdl_constr,
                       _data->ir,
                       _data->mdatoms,
                       _data->state,
                       _data->upd,
                       _data->fr->bMolPBC,
                       &_data->top->idef,
                       _data->shake_vir,
                       _data->cr,
                       _data->nrnb,
                       _data->constr,
                       _data->bCalcVir, 1,
                       as_rvec_array(_data->state->v.data()));
    set_x(_data->mdatoms, _data->upd, _data->state);
}

void ConstrainVelocity::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void ConstrainVelocity::run()
{
    set_vprime(_data->mdatoms, _data->upd, _data->state);
    // Assumes that x and v are on the same time step
    constrain_velocities(_data->step,
                         &_data->dvdl_constr,
                         _data->ir,
                         _data->mdatoms,
                         _data->state,
                         _data->upd,
                         _data->fr->bMolPBC,
                         &_data->top->idef,
                         _data->shake_vir,
                         _data->cr,
                         _data->nrnb,
                         _data->constr,
                         _data->bCalcVir, 1);
    set_v(_data->mdatoms, _data->upd, _data->state);
}