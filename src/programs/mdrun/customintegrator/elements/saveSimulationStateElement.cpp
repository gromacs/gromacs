#include "saveSimulationStateElement.h"

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdtypes/commrec.h"
#include "programs/mdrun/customintegrator/statemanager.h"
#include <cstdio>

void SaveSimulationState::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void SaveSimulationState::run()
{
	/*
	Copy value of state and communication record structures to simulationState structure. 
	Store them in a map with string labels as key and the simulationState structure as value
	The string label is given as a variable upon construction of the element.
	
	**Future**: Instead of saving the whole structures only save the chaged variables in the struture

	To store info required for reproducible parallel runs store the cr struct
	More specifically the following might be stored from cr, instead of the whole cr.
	gatindex : array of indices from local to global atoms
	gatindex_nalloc: the length of the array
	nat_home : number of home atoms, Is this the same as mdatoms->homenr
	If it is the same, mdatoms->homenr might be redundant nat_home can be used instead.
	*/
	// store the state_t structure
	_data->simstate_map[label].state = *(_data->state);
	// store the t_commrec structure
	_data->simstate_map[label].cr = *(_data->cr);
}