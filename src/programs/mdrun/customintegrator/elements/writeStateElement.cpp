#include "writeStateElement.h"

// mdrun.h is dependent on vec.h, thus needs to be included before mdrun.h
// This is not a good way, may report to gromacs developers
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void WriteState::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void WriteState::run()
{
     do_md_trajectory_writing(_data->fplog, _data->cr, _data->nfile, _data->fnm, _data->step, _data->step_rel, _data->t,
                             _data->ir, _data->state, _data->state_global, _data->energyHistory,
                             _data->top_global, _data->fr,
                             _data->outf, _data->mdebin, _data->ekind, &_data->f,
                             &_data->nchkpt,
                             _data->bCPT, _data->bLastStep, (_data->Flags & MD_CONFOUT),
                             _data->bSumEkinhOld);
}