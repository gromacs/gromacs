#ifndef _drawVelocitiesFromGaussianElement_h
#define _drawVelocitiesFromGaussianElement_h

#include "elementBase.h"
#include <string>

class Element;
class StateManager;

/* An Element that draws velocities from gaussian distribution
with mean=0 and sigma=sqrt(kT/m)

Given the same seed, the random number generations is reproducible even with domain decomposition

*/
class DrawVelocitiesFromGaussian : public Element
{
public:
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager class
    StateManager* _data;
};

#endif