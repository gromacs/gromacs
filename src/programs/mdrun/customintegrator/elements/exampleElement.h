#ifndef _exampleElement_h
#define _exampleElement_h

#include "elementBase.h"

class StateManager;

// An example Element class, that can be modified to write your own element
class Example : public Element
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
