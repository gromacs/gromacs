#include "phase.h"

static const char *phases[epNR] = { "gas", "liquid", "solid", "plasma" };
    
std::string phase2string(ePhase ep)
{
    return phases[ep];
}

ePhase string2phase(std::string phase)
{
    for(int i = 0; (i<epNR); i++)
    {
        if (0 == phase.compare(phases[i]))
        {
            return (ePhase) i;
        }
    }
    return epNR;
}


