#ifndef _metropolisAcceptanceTestElement_h
#define _metropolisAcceptanceTestElement_h

#include "elementBase.h"
#include <string>

class Element;
class StateManager;

/*
   An element that performs a Metropolis acceptance test

   Compares a saved state (determined by `label1`) to the current state (denoted
   by `label2`) using the Metropolis criterion. If the acceptance test fails,
   the element is restoring the saved state.
*/
class MetropolisAcceptanceTest : public Element
{
public:
    // Constructor that makes the element run every step
    MetropolisAcceptanceTest(std::string label1, std::string label2):
                            label1(label1),
                            label2(label2) {};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
    // Restore state given by the label
    // This should be moved to a base AcceptanceTest class
    void restoreState(std::string label);
private:
    // pointer to a StateManager class
    StateManager* _data;
    std::string label1, label2;
};

#endif