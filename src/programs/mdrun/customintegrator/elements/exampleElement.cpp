#include "exampleElement.h"

#include "programs/mdrun/customintegrator/statemanager.h"

void Example::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void Example::run()
{
    // Write the code in this run funtion 
    // that you want the element to execute during integrator block
}
