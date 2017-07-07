#ifndef _writeStateElement_h
#define _writeStateElement_h

#include "elementBase.h"

class Element;
class StateManager;

/* 
    An element that writes out the state to xtc, trr, gro files
    
    ***Currently***: 
    The element calls the "do_md_trajectory_writing" that is present in the 
    standard Gromacs code. Thus depending on the user inputs, this element
    can write multiple information to various files such as xtc, trr, gro
    and checkpoint files

    ***Future***:
    We might not have an element as part of the integrator that does IO
    instead we will have the "save" elements that will save relevant information
    which the IO function can use to write to files. This IO could be done at the
    end of the integrator step as an independent element to the integrator algorithm.
*/
class WriteState : public Element
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