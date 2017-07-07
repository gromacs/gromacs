#include "saveEnergyElement.h"

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdtypes/commrec.h"
#include "programs/mdrun/customintegrator/statemanager.h"
#include <cstdio>

void SaveEnergy::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void SaveEnergy::run()
{
	// store the gmx_enerdata_t structure
	_data->simstate_map[label].enerd = *(_data->enerd);
}