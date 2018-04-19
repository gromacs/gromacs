//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_PMEELEMENT_H
#define GROMACS_PMEELEMENT_H

#include "element.h"

namespace gmx
{

class PMELoadBalanceElement : public Element
{
    public:
        PMELoadBalanceElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}

        void loopSetup() override;
        void run() override;
        void loopTeardown() override;
};

}


#endif //GROMACS_PMEELEMENT_H
