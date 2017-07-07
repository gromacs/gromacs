#ifndef _helloWorldElement_h
#define _helloWorldElement_h

#include "elementBase.h"

class StateManager;

// An element printing "Hello World! (step = , dt = )" to log file
// Useless on a long run, very useful in debugging
class HelloWorld : public Element
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