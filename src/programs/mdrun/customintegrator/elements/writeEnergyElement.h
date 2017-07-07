#ifndef _writeEnergyElement_h
#define _writeEnergyElement_h

#include "elementBase.h"

class Element;
class StateManager;

// An element that writes out the energies to log and edr files.
// It uses the mdebin structure to write out the energy related information
class WriteEnergy : public Element
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