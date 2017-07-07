#ifndef _sequenceFactory_h
#define _sequenceFactory_h

#include <list>
#include "elements/elementBase.h"

struct t_inputrec;

/*
    The SequenceFactory creates sequences of Elements from given input

    Three types of sequences are generated:
    * pre-run, ran once before the integrator loop,
    * run , repeatedly ran within the integrator, and
    * post-run, ran once after the end of the integrator loop.

    The sequence of elements is determined by the input options (t_inputrec plus
    additional bConstrain boolean).

    ***In the future***, we envision to offer a general integrator parsing, which
    would allow users to determine an integration algorithm by listing the
    required elements without the need of code-changes or recompiling.

    ***Currently***: the create{PreRun,Run,PostRun}Sequence() functions return
    sequences of elements essentially hard-coded based on the value of the
    `custom-type` flag and some additional mdp options.
*/

class SequenceFactory
{
public:
    static std::list<Element*> createRunSequence(t_inputrec *ir,
                                                 bool bConstrain);
    static std::list<Element*> createPreRunSequence(t_inputrec *ir,
                                                    bool bConstrain);
    static std::list<Element*> createPostRunSequence(t_inputrec *ir,
                                                     bool bConstrain);
};
#endif