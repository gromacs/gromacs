#ifndef _updateForcesElement_h
#define _updateForcesElement_h

#include "elementBase.h"
#include "gromacs/utility/basedefinitions.h"

class Element;
class StateManager;

// An element that calculates the forces, potential energy, and virial given current positions
class UpdateForces : public Element
{
public:
    UpdateForces(int nstrun, bool runAt):Element(nstrun,runAt){};
    UpdateForces(){};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();

private:
    // pointer to a StateManager class
    StateManager* _data;
};

#endif