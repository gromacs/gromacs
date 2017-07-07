#ifndef _saveSimulationStateElement_h
#define _saveSimulationStateElement_h

#include "elementBase.h"
#include <string>

class Element;
class StateManager;

/* 
An Element that saves the simulationState and attaches a label to it, allowing
it to be restored by using the label as a key. If the label already exists, 
it will overwrite the state stored by that label.
*/
class SaveSimulationState : public Element
{
public:
    // Constructor that enables the element to be run at every nstrun step
    SaveSimulationState(int nstrun, std::string label):Element(nstrun), label(label){};
    // Constructor than makes the element run every step
    SaveSimulationState(std::string label):label(label){};
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