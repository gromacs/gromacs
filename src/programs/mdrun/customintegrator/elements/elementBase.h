#ifndef _elementBase_h
#define _elementBase_h


class StateManager;

/*
An Abstract class for all element classes.
*/
class Element
{
public:
    Element():nstrun(1), runAt(false){};
    Element(int nstrun):nstrun(nstrun), runAt(false){};
    Element(int nstrun, bool runAt):nstrun(nstrun), runAt(runAt){};

    // initialize the element with the data in StateManager object
    virtual void initialize(StateManager& dataRef)=0;
    // executes the element's function
    virtual void run()=0;
   
    int get_nstrun(){return nstrun;}
    bool get_runAt(){return runAt;}
    virtual ~Element(){};
private:
    // nstrun : the stride of the element with respect to the integrator step
    int nstrun;
    // if runAt is True, nstrun can be interpreted as the
    // step at whic the element should be executed
    bool runAt;
};
#endif