#ifndef _integratorFactory_h
#define _integratorFactory_h

#include <memory>

#include "elements/elementBase.h"
#include "elements/integratorElement.h"

class Integrator;

class IntegratorFactory
{
public:
    /*
    Given the runSequence and preRunSequence, this method creates an
    Integrator object. Both sequences are lists of pointers of elements, which
    are added to the respective sequence of the returned Integrator object.

    The preRunSequence and runSequence fully define the integrator.
    */
    static std::unique_ptr<Integrator> create(std::list<Element*> runSequence,
                                              std::list<Element*> preRunSequence,
                                              std::list<Element*> postRunSequence);
};
#endif
