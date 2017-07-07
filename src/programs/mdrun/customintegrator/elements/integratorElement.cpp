#include "integratorElement.h"

#include "programs/mdrun/customintegrator/statemanager.h"

Integrator::~Integrator()
{
   for (auto&& e : _element)
    {
        delete e;
        e = nullptr;
    } 
    _element.clear();

    for (auto&& e : _preRunElement)
    {
        delete e;
        e = nullptr;
    } 
    _preRunElement.clear();
}

void Integrator::run()
{
    bool runAt;
    int nstrun;
    for (auto&& e : _element)
    {
        runAt = e->get_runAt();
        nstrun = e->get_nstrun();
        if (!runAt)
        {
            if (counter % nstrun == 0)
            {
                e->run();
            }
        }
        else
        {
            if (counter == nstrun)
            {
                e->run();
            }
        }

    }
    counter += 1;
}

void Integrator::prerun()
{
    for (auto&& e : _preRunElement)
    {
        e->run();
    }
}

void Integrator::postrun()
{
    for (auto&& e : _postRunElement)
    {
        e->run();
    }
}

void Integrator::addRunElement(Element *element)
{
    _element.push_back(element);
}

void Integrator::addPreRunElement(Element *element)
{
    _preRunElement.push_back(element);
}

void Integrator::addPostRunElement(Element *element)
{
    _postRunElement.push_back(element);
}

void Integrator::initialize(StateManager& dataRef)
{
    for (auto&& e : _element)
    {
        e->initialize(dataRef);
    }

    for (auto&& e : _preRunElement)
    {
        e->initialize(dataRef);
    }

    for (auto&& e : _postRunElement)
    {
        e->initialize(dataRef);
    }
}
