#ifndef _ouoperatorElement_h
#define _ouoperatorElement_h

#include "elementBase.h"

#include <vector>
#include <gromacs/utility/real.h>
#include "gromacs/utility/basedefinitions.h"

class Element;
class StateManager;
// An element implementing the Ornstein-Uehlenbeck operator
class OUOperator : public Element
{
public:
    OUOperator(real step_fraction = 1.0) :
            _step_fraction(step_fraction)
    {};
    // initialize the element with the data in StateManager object
    void initialize(StateManager& dataRef);
    // executes the element's function
    void run();
private:
    // pointer to a StateManager object
    StateManager* _data;
    real _step_fraction;
    std::vector<real> _a;
    std::vector<real> _sqrt_one_a2_beta;
};

#endif