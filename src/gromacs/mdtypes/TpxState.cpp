//
// Created by Eric Irrgang on 8/27/17.
//

#include "TpxState.h"

#include <memory>
#include <string>

#include "gromacs/compat/make_unique.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

TpxState::TpxState() :
    filename_ {},
inputrecInstance_ {
    std::make_shared<t_inputrec>()
},
stateInstance_ {
    std::make_shared<t_state>()
},
mtop_ {
    nullptr
},
initialized_ {
    false
},
dirty_ {
    false
}
{
    snew(mtop_, 1);
}

TpxState::~TpxState()
{
    if (mtop_ != nullptr)
    {
        sfree(mtop_);
    }
};

std::unique_ptr<TpxState> TpxState::initializeFromFile(const char* filename)
{
    std::string arg;
    arg = filename;
    return initializeFromFile(arg);
}

std::unique_ptr<TpxState> TpxState::initializeFromFile(const std::string &filename)
{
    auto newState = gmx::compat::make_unique<TpxState>();
    read_tpx_state(filename.c_str(), newState->inputrecInstance_.get(), newState->stateInstance_.get(), newState->mtop_);

    newState->filename_    = filename;
    newState->initialized_ = true;
    return newState;
}

t_inputrec *TpxState::getRawInputrec()
{
    dirty_ = true;
    return inputrecInstance_.get();
}

gmx_mtop_t *TpxState::getRawMtop()
{
    dirty_ = true;
    return mtop_;
}

t_state *TpxState::getRawState()
{
    dirty_ = true;
    return stateInstance_.get();
}

bool TpxState::isInitialized() const
{
    return initialized_;
}

std::unique_ptr<TpxState>
TpxState::initializeFromWrappers(std::unique_ptr<t_inputrec> inputRecord, std::unique_ptr<t_state> state,
                                 std::unique_ptr<gmx_mtop_t> mtop)
{
    auto newState = gmx::compat::make_unique<TpxState>();
    newState->inputrecInstance_ = std::move(inputRecord);
    newState->stateInstance_    = std::move(state);
    newState->mtop_             = mtop.release();
    return newState;
}

TpxState::TpxState(TpxState &&source) noexcept
{
    if (this != &source)
    {
        std::lock_guard<std::mutex> lock(source.exclusive_);
        filename_         = std::move(source.filename_);
        inputrecInstance_ = std::move(source.inputrecInstance_);
        stateInstance_    = std::move(source.stateInstance_);
        mtop_             = source.mtop_;
        initialized_.store(source.initialized_.load());
        dirty_.store(source.dirty_.load());

        // Make an effort to invalidate the old object in case there are outstanding
        // handles (which constitute a bug that we can make noisier in future revisions.)
        source.mtop_        = nullptr;
        source.initialized_ = false;
        source.dirty_       = true;
    }
}

TpxState &TpxState::operator=(TpxState &&source) noexcept
{
    std::lock_guard<std::mutex> lockDestination(this->exclusive_);
    if (this != &source)
    {
        std::lock_guard<std::mutex> lockSource(source.exclusive_);
        filename_         = std::move(source.filename_);
        inputrecInstance_ = std::move(source.inputrecInstance_);
        stateInstance_    = std::move(source.stateInstance_);
        mtop_             = source.mtop_;
        initialized_.store(source.initialized_.load());
        dirty_.store(source.dirty_.load());

        // Make an effort to invalidate the old object in case there are outstanding
        // handles (which constitute a bug that we can make noisier in future revisions.)
        source.mtop_        = nullptr;
        source.initialized_ = false;
        source.dirty_       = true;
        // Todo: Notifications to subscribed objects
    }

    return *this;
}

bool TpxState::isDirty() const
{
    return dirty_;
}

void TpxState::markClean()
{
    dirty_ = false;
}

const char* TpxState::filename() const
{
    return filename_.c_str();
}

} // end namespace gmx
