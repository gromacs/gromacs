//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_CONSTRAINELEMENT_H
#define GROMACS_CONSTRAINELEMENT_H

#include "element.h"

namespace gmx
{

class ConstrainElement : public Element
{
    public:
        ConstrainElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}


#endif //GROMACS_CONSTRAINELEMENT_H
