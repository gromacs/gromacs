#include "gmxapi/runner.h"
#include "gmxapi/md/runnerstate.h"
#include "gmxapi/md.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/TpxState.h"

#include "gmxapi/system.h"
#include "gmxapi/context.h"
#include <gtest/gtest.h>

namespace
{

class DummyMD : public gmxapi::MDEngine
{
};

TEST(ApiRunner, Build)
{
    auto md            = std::make_shared<gmxapi::MDEngine>();
    auto runnerBuilder = gmxapi::UninitializedMDRunnerState::Builder();
    runnerBuilder.mdEngine(md);
    runnerBuilder.tpxState(std::make_shared<gmx::TpxState>());
    auto runner = runnerBuilder.build();
//    auto session = runner->initialize(gmxapi::defaultContext());
//    ASSERT_TRUE(session != nullptr);
//    auto status = session->run();
    // Just make sure we made it this far...
//    ASSERT_TRUE(!status.success());
}

TEST(ApiRunner, BasicMD)
{
    // Need to set up a test fixture...
    const std::string filename = "${CMAKE_CURRENT_BINARY_DIR}/topol.tpr";

    auto              system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
        ASSERT_TRUE(context != nullptr);
        ASSERT_TRUE(system != nullptr);
        ASSERT_TRUE(system->runner() != nullptr);
        auto           runner  = system->runner();
        auto           session = runner->initialize(context);
        ASSERT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
//        ASSERT_NO_THROW(session->run(1000));
        ASSERT_TRUE(status.success());
    }

}

}
