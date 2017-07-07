#include "integratorFactory.h"

/*
    Given the runSequence and preRunSequence, this method creates an
    Integrator object. Both sequences are lists of pointers of elements, which
    are added to the respective sequence of the returned Integrator object.

    The preRunSequence and runSequence fully define the integrator.
*/
std::unique_ptr<Integrator>
IntegratorFactory::create(std::list<Element*> runSequence,
                          std::list<Element*> preRunSequence,
                          std::list<Element*> postRunSequence)
{
    std::unique_ptr<Integrator> integrator (new Integrator);
    
    for (auto&& el : runSequence)
    {
        integrator->addRunElement(el);
    }
    for (auto&& el : preRunSequence)
    {
        integrator->addPreRunElement(el);
    }
    for (auto&& el : postRunSequence)
    {
        integrator->addPostRunElement(el);
    }

    return integrator;
}
