//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_WRITEENERGYELEMENT_H
#define GROMACS_WRITEENERGYELEMENT_H

#include "element.h"

namespace gmx
{

class WriteEnergyElement : public Element
{
    public:
        WriteEnergyElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}


#endif //GROMACS_WRITEENERGYELEMENT_H
