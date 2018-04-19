//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_UPDATEFUSEDELEMENT_H
#define GROMACS_UPDATEFUSEDELEMENT_H

#include "element.h"

namespace gmx
{

class UpdateFusedElement : public Element
{
    public:
        UpdateFusedElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}


#endif //GROMACS_UPDATEFUSEDELEMENT_H
