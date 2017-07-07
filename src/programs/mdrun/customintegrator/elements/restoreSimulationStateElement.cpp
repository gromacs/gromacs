#include "restoreSimulationStateElement.h"

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/state.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void RestoreSimulationState::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void RestoreSimulationState::run()
{
	*(_data->cr) = _data->simstate_map.at(label).cr;
	*(_data->state) = _data->simstate_map.at(label).state;
}