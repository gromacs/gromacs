#ifndef _integratorElement_h
#define _integratorElement_h

#include <list>
#include "gromacs/utility/basedefinitions.h"
#include "elementBase.h"

class Element;
class StateManager;

/* 
Integrator Element
The integrator is a *composite element* that has a list of pointers to other
elements.
The composite element can store pointers to other composite elements to allow
for nested (multi-timestep) integrators (or, in principle, it could store a
pointer to itself - maybe we should catch that?)

The Integrator Composite Element has three lists:

    * _element       : Elements executed every integrator step
    * _preRunElement : Elements executed once before the start of the integrator
    * _postRunElement: Elements executed once after the end of the integrator
*/
class Integrator : public Element
{
public:
    Integrator():counter(0){};

    // iterates through the list of prerun elements and executes their
    // respective run functions
    void prerun();
    // iterates through the list of elements and executes their run function
    void run();
    // iterates through the list of postrun elements and executes their
    // respective run functions
    void postrun();
    // add a new element to the list
    void addRunElement(Element* element);
    // add a new element to the list
    void addPreRunElement(Element* element);
    // add a new element to the list
    void addPostRunElement(Element* element);
    // initialize elements in the list
    void initialize(StateManager& dataRef);
    // destructor
    ~Integrator();

private:
    // list of pointers to element objects
    std::list<Element*> _element;
    std::list<Element*> _preRunElement;
    std::list<Element*> _postRunElement;
    int counter;
};

#endif