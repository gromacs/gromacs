#ifndef GMXAPI_SYSTEM_IMPL_H
#define GMXAPI_SYSTEM_IMPL_H

#include "gmxapi/system.h"

namespace gmxapi
{

class Atoms;
class MDInput;
class IMDRunner;

class System::Impl
{
    public:
        Impl();
        ~Impl();
//        std::unique_ptr<Atoms> atoms();
//        void setAtoms(const Atoms& atoms);
        std::shared_ptr<MDEngine> md();
        void md(std::shared_ptr<MDEngine> md);

        std::shared_ptr<IMDRunner> &runner();

        void runner(std::shared_ptr<IMDRunner> runner);

    private:
        std::shared_ptr<MDEngine>  md_;
        std::shared_ptr<IMDRunner> runner_;
//        std::unique_ptr<Atoms>  atoms_;

};

}      // end namespace gmxapi

#endif // header guard
