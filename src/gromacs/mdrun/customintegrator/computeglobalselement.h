//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_COMPUTEGLOBALSELEMENT_H
#define GROMACS_COMPUTEGLOBALSELEMENT_H

#include "element.h"

namespace gmx
{

class ComputeGlobalsElement : public Element
{
    public:
        ComputeGlobalsElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}

#endif //GROMACS_COMPUTEGLOBALSELEMENT_H
