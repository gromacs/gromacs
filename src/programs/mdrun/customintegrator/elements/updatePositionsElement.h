#ifndef _updatePositionsElement_h
#define _updatePositionsElement_h

#include "elementBase.h"

#include "gromacs/utility/real.h"

class Element;
class StateManager;

/*
    An element that progresses the positions by a given fraction of timestep
    The fraction of the time step is given during the construction
*/
class UpdatePositions : public Element
{
public:
    UpdatePositions(real step_fraction = 1.0):step_fraction(step_fraction){};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager class
    StateManager* _data;
    real step_fraction, dt;
};

#endif