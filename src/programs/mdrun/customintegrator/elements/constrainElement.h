#ifndef _constrainElement_h
#define _constrainElement_h

#include "elementBase.h"

class Element;
class StateManager;

/* 
An element that constrains the initial positions and velocities
Depending on the time step of the velocities (half or fullstep)
the constraining of velocities needs to be done in a different way
*/
class InitialConstrain : public Element
{
public:
    InitialConstrain(bool halfStep):halfStep(halfStep){};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager class
    StateManager* _data;
    bool halfStep;
};

// An element that constrains fullstep positions and halfstep velocities
class ConstrainPositionAndVelocity : public Element
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

// An element that constrains the fullstep velocities only
class ConstrainVelocity : public Element
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