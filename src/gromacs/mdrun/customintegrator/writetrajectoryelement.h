//
// Created by Pascal Merz on 4/21/18.
//

#ifndef GROMACS_WRITETRAJECTORYELEMENT_H
#define GROMACS_WRITETRAJECTORYELEMENT_H

#include "element.h"

namespace gmx
{

class WriteTrajectoryElement : public Element
{
    public:
        WriteTrajectoryElement(Integrator &integ, DataManager &dm) : Element(integ, dm) {}
        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}
};

}


#endif //GROMACS_WRITETRAJECTORYELEMENT_H
