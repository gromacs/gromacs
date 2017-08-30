#ifndef GMXAPI_SYSTEM_BUILDER_H
#define GMXAPI_SYSTEM_BUILDER_H

#include "gmxapi/system.h"
#include "gmxapi/md.h"

namespace gmxapi
{

class MDInput;
class Topology;

class System::Builder
{
    public:
        Builder();
        ~Builder();

        // Allow an appropriate default Context to be determined and configured.
//        Builder &defaultContext(const MDInput &inputrec);

        // Use the information in the input record to configure an appropriate runner.
        Builder &runner(std::shared_ptr<IMDRunner> runner);

//        Builder &structure(const MDInput &inputrec);

        Builder &mdEngine(std::shared_ptr<MDEngine> md);

//        Builder &topology(std::shared_ptr<Topology> topology);

        // Pass ownership of the assembled System.
        std::unique_ptr<System> build();
    private:
        std::unique_ptr<System> system_;
};

}      // end namespace gmxapi

#endif // header guard
