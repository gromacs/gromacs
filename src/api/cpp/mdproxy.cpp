//
// Created by Eric Irrgang on 7/31/17.
//

#include <memory>
#include <functional>
#include <cassert>

#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"

#include "gromacs/compat/make_unique.h"
#include "md-impl.h"
#include "api/cpp/gmxapi/md/runnerstate.h"

namespace gmxapi
{

/*!
 * \brief Represent empty or new MDEngine instance.
 *
 * Not necessary while MDEngine is instantiable.
 */
class MDNull : public MDEngine
{
};

MDStateFromMDInput::MDStateFromMDInput() : input_ {nullptr}
{}

MDStateFromMDInput::MDStateFromMDInput(std::unique_ptr<MDInput> input) :
    MDStateFromMDInput {std::move(input), std::string()}
{}

MDStateFromMDInput::MDStateFromMDInput(std::unique_ptr<MDInput> input, std::string metadata) :
    input_ {std::move(input)},
metadata_ {
    metadata
}
{}

MDStateFromMDInput::~MDStateFromMDInput() = default;


MDProxy::MDProxy() :
    instanceState_ {std::make_shared<MDEngine>()}
{};
MDProxy::MDProxy(MDProxy &&) noexcept            = default;
MDProxy &MDProxy::operator=(MDProxy &&) noexcept = default;

void MDProxy::setState(std::shared_ptr<MDEngine> state)
{
    instanceState_ = std::move(state);
}

const std::string MDProxy::info() const
{
    assert(instanceState_ != nullptr);
    return instanceState_->info();
}

std::unique_ptr<MDBuilder> MDStateFromMDInput::builder()
{
    return nullptr;
}

const std::string MDStateFromMDInput::info() const
{
    if (input_ == nullptr)
    {
        return "uninitialized MDStateFromMDInput";
    }
    else
    {
        std::string output("MDStateFromMDInput initialized");
        if (!metadata_.empty())
        {
            output.append(" with metadata: ");
            output.append(metadata_);
        }
        return output;
    }
}

std::unique_ptr<MDBuilder> MDProxy::builder()
{
    assert(instanceState_ != nullptr);
    return instanceState_->builder();
}

} // end namespace gmxapi
