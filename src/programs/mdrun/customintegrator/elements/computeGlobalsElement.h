#ifndef _computeGlobalsElement_h
#define _computeGlobalsElement_h

#include "elementBase.h"
#include "gromacs/utility/basedefinitions.h"

class Element;
class StateManager;

// An example Element class, that can be modified to write your own element
class ComputeGlobals : public Element
{
public:
    ComputeGlobals(gmx_bool bEnergy,
                   gmx_bool bReadEkin,
                   gmx_bool bScaleEkin,
                   gmx_bool bEkinAveVel,
                   gmx_bool bStopCM,
                   gmx_bool bTemp,
                   gmx_bool bPres,
                   gmx_bool bConstraint,
                   bool bInterSimSignal,
                   bool bIntraSimSignal);

    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager class
    StateManager* _data;
    int cglo_flags;
    bool bIntraSimSignal, bInterSimSignal;
};

#endif