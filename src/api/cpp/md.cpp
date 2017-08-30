#include <memory>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "gmxapi/runner.h"

#include "gromacs/compat/make_unique.h"
#include "md-impl.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetree.h"

namespace gmxapi
{

MDEngine::~MDEngine() = default;

const std::string MDEngine::info() const
{
    return "Generic MDEngine object";
}

std::unique_ptr<MDBuilder> MDEngine::builder()
{
//    std::unique_ptr<MDBuilder> builder = std::make_unique<MDEngineBuilder>();
//    return builder;
    class DummyBuilder : public MDBuilder
    {
        public:
            virtual ~DummyBuilder() override = default;

            virtual std::unique_ptr<MDEngine> build() override
            {
                std::unique_ptr<MDEngine> proxy = gmx::compat::make_unique<MDProxy>();
                return proxy;
            }
    };
    std::unique_ptr<MDBuilder> builder = gmx::compat::make_unique<DummyBuilder>();
    return builder;
}

void MDEngine::bind(IMDRunner* runner)
{
    auto builder = this->builder();
    runner->registerMDBuilder(std::move(builder));
}

MDInput::MDInput() :
    inputRecord_ {gmx::compat::make_unique<t_inputrec>()},
state_ {
    gmx::compat::make_unique<t_state>()
},
topology_ {
    gmx::compat::make_unique<gmx_mtop_t>()
}
{
}

MDInput::MDInput(std::unique_ptr<t_inputrec> &&inputRecord,
                 std::unique_ptr<t_state>    &&state,
                 std::unique_ptr<gmx_mtop_t> &&topology) :
    inputRecord_ {std::move(inputRecord)},
state_ {
    std::move(state)
},
topology_ {
    std::move(topology)
}
{
    set_state_entries(state_.get(), inputRecord_.get());
}

int MDInput::nAtoms() const
{
    if (state_ != nullptr)
    {
        return state_->natoms;
    }
    else
    {
        return 0;
    }
}

gmx::KeyValueTreeObject MDInput::params() const
{
    return *inputRecord_->params;
}

const t_state* MDInput::state() const
{
    return state_.get();
}

std::unique_ptr<MDInput> MDInput::from_tpr_file(std::string filename)
{
    auto inputRecord = gmx::compat::make_unique<t_inputrec>();
    auto state       = gmx::compat::make_unique<t_state>();
    auto topology    = gmx::compat::make_unique<gmx_mtop_t>();

    // Don't know much about read_tpx_state right now, so I'd rather see the
    // exceptions raised...
    //try
    //{
    // Fill the output parameters.
    read_tpx_state(filename.c_str(), inputRecord.get(), state.get(), topology.get());
    //}
    //catch (std::exception& e)
    //{
    //    return nullptr;
    //}

    return gmx::compat::make_unique<MDInput>(std::move(inputRecord), std::move(state), std::move(topology));
}


// This function currently creates data structures from the TPR file, implying a unique
// simulation state. Typically a proxy object would not uniquely own such an instance,
// but in this case it is the data from which to instantiate a runner and MD Engine.
// However, this makes the proxy object non-trivially copyable.
// \todo Data structures used to hold API object state need to be copyable and so should
// be small and/or references to data managed elsewhere.
//std::unique_ptr<MDProxy> mdFromTpr(const std::string filename)
//{
//    using gmx::compat::make_unique;
//    auto tprInput = MDInput::fromTprFile(filename);
//
//    // Transfer ownership of input to new state object
//    auto newState = make_unique<MDStateFromMDInput>(std::move(tprInput), filename);
//    // Get new proxy object and transfer ownership of state.
//    auto md = make_unique<MDProxy>();
//    md->setState(std::move(newState));
//    return md;
//}

std::unique_ptr<MDProxy> mdFromTpr(const std::string filename)
{
    auto newState = gmx::compat::make_unique<MDStatePlaceholder>(filename);
    auto md       = gmx::compat::make_unique<MDProxy>();
    md->setState(std::move(newState));
    return md;
}

} //end namespace gmxapi
