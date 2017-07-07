#ifndef _saveEnergyElement_h
#define _saveEnergyElement_h

#include "elementBase.h"
#include <string>

class Element;
class StateManager;

/* 
An Element that saves the enerd structure to a simulationState structure with a given label,
which enables to restore it by the label
If the label already exists it will overwrite the enerd struct value of that label
*/
class SaveEnergy : public Element
{
public:
	// Constructor that enables the element to be run at every nstrun step
	SaveEnergy(int nstrun, std::string label):Element(nstrun), label(label){};
	// Constructor than makes the element run every step
	SaveEnergy(std::string label):label(label){};
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