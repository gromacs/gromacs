//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_DDELEMENT_H
#define GROMACS_DDELEMENT_H

#include "element.h"

namespace gmx
{

class DDElement : public Element
{
    public:
        DDElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}

#endif //GROMACS_DDELEMENT_H
