//
// Created by pascal on 4/20/18.
//

#include "element.h"

#ifndef GROMACS_FORCEELEMENT_H
#define GROMACS_FORCEELEMENT_H

namespace gmx
{
class ForceElement : public Element
{
    public:
        ForceElement(Integrator &integ, DataManager &dm);

        void loopSetup() override {}
        void run() override;
        void loopTeardown() override {}

    private:
        bool doPCoupling;
        int  nstPCouple;
        int  nstCalcEnergy;
        DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion;
        DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion;

};
}

#endif //GROMACS_FORCEELEMENT_H
