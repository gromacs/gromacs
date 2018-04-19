//
// Created by pascal on 4/20/18.
//

#ifndef GROMACS_ELEMENT_H
#define GROMACS_ELEMENT_H

#include "gromacs/mdrun/integrator.h"
#include "gromacs/mdrun/customintegrator/datamanager.h"

namespace gmx
{
class Element
{
    public:
        Element(Integrator &integ, DataManager &dm) :
            integrator(integ), dataManager(dm) {};

        virtual void loopSetup()    = 0;
        virtual void run()          = 0;
        virtual void loopTeardown() = 0;

        virtual ~Element() = default;

    protected:
        Integrator  &integrator;
        DataManager &dataManager;
};
}
#endif //GROMACS_ELEMENT_H
