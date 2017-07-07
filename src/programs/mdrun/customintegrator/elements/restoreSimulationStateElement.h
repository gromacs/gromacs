#ifndef _restoreSimulationStateElement_h
#define _restoreSimulationStateElement_h

#include "elementBase.h"
#include <string>

class Element;
class StateManager;
/* An Element that restores a simulationState given by a label during
constructiom of the element
*/
class RestoreSimulationState : public Element
{
public:
	// Constructor that enables the element to be run at every nstrun step
	RestoreSimulationState(int nstrun, std::string label):Element(nstrun), label(label){};
	// Constructor than makes the element run every step
	RestoreSimulationState(std::string label):label(label){};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager class
    StateManager* _data;
    std::string label;
};

#endif