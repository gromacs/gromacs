#include "metropolisAcceptanceTestElement.h"
#include "restoreSimulationStateElement.h"

#include "gromacs/math/units.h"
#include "gromacs/utility/real.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "programs/mdrun/customintegrator/statemanager.h"

void MetropolisAcceptanceTest::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void MetropolisAcceptanceTest::run()
{
	/*
	deltaH = TotalEnerg(B) - TotalEnergy(A)
	A = exp(-beta*deltaH)
	if Rand() > A
		RejectState
	else:
		AcceptState
	*/

	// Create a default random number generator using Gromacs random library
	gmx::DefaultRandomEngine rng(unsigned(_data->ir->ld_seed));
    gmx::UniformRealDistribution<real> uniform;
	real A;
    
    /* Note: the temperature information is obtained from mdp
    * It will take the first temperature group's temperature to 
    * draw velocities for the whole system.
    * It might be worth re-thinking what the code should do if 
    * there are multiple temperature groups.
    * Or introducing another temperature field in the mdp for HMC integrators
    * may be more useful*/

	real kT = real(BOLTZ*_data->ir->opts.ref_t[0]);
    real TotalEnergyA = real(_data->simstate_map.at(label1).enerd.term[F_ETOT]);
    real TotalEnergyB = real(_data->simstate_map.at(label2).enerd.term[F_ETOT]);
    real deltaH = TotalEnergyB - TotalEnergyA;
    A = exp(-deltaH/kT);
    
    // Throw a Dice with probability P = min(1,A)
    if (uniform(rng) > A)
    {
    	// Reject the current state and restore old state
    	printf("State rejected, restoring state A!\n");
    	restoreState(label1);
    }
    else
    {
        // Accept state
    	restoreState(label2);
    }
}

void MetropolisAcceptanceTest::restoreState(std::string label)
{
	// To reduce code duplication, this function executes 
	// the run method of the RestoreSimulationState Element

	RestoreSimulationState restoreSimState(label);
	restoreSimState.initialize(*_data);
	restoreSimState.run();
}
