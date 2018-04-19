//
// Created by pascal on 4/20/18.
//

#ifndef GROMACS_LOOPELEMENT_H
#define GROMACS_LOOPELEMENT_H

#include <vector>
#include "element.h"
#include "datamanager.h"

namespace gmx
{
class LoopElement : public Element
{
    public:
        LoopElement(Integrator &integ, DataManager &dm,
                    int reps = 1, bool propagateStep = false) :
            Element(integ, dm), repetitions(reps), inner_loop(propagateStep) {}

        void loopSetup() override {}

        void run() override
        {
            for (auto &element : elements)
            {
                element->loopSetup();
            }
            for (int r = 0; r < repetitions; r++)
            {
                if (!inner_loop || r > 0)
                {
                    // if it's an inner loop, the data manager has been updated
                    // before entering
                    dataManager.preStep();
                }
                for (auto &element : elements)
                {
                    element->run();
                }
                if (!inner_loop || r < repetitions - 1)
                {
                    // if it's an inner loop, the data manager will be updated
                    // after exiting
                    dataManager.postStep();
                }
            }
            for (auto &element : elements)
            {
                element->loopTeardown();
            }
        }

        void loopTeardown() override {}

        ~LoopElement() override
        {
            // TODO: is this the place to delete elements?
            // Need to think about ownership and lifetime of elements
            // Currently we're definitely leaking memory
            // for (auto &element : elements)
            // {
            //     delete element;
            // }
        }

        void add(Element* el)
        {
            GMX_ASSERT(el != this, "Cannot recursively add loop elements.");
            elements.push_back(el);
        }

    protected:
        int                   repetitions;
        bool                  inner_loop;
        std::vector<Element*> elements;

};
}
#endif //GROMACS_LOOPELEMENT_H
