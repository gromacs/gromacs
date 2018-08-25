//
// Created by Eric Irrgang on 8/17/18.
//

#ifndef GROMACS_MD_IMPL_H
#define GROMACS_MD_IMPL_H

#include "md.h"

namespace gmx
{

class MDIntegrator::Impl: public IntegratorParamsContainer
{
    public:
        template<class... Args> Impl(Args&&...args) :
        IntegratorParamsContainer{std::forward<Args>(args)...} {};

        void run();
};

} // end namespace gmx

#endif //GROMACS_MD_IMPL_H
