#include "helloWorldElement.h"

#include "gromacs/utility/logger.h"

#include "programs/mdrun/customintegrator/statemanager.h"

void HelloWorld::initialize(StateManager& dataRef)
{
    _data = &dataRef;
}

void HelloWorld::run()
{
    GMX_LOG(_data->mdlog->info).appendTextFormatted("Hello World! (step = %d, dt = %f)",
                                                    _data->step, _data->ir->delta_t);
}