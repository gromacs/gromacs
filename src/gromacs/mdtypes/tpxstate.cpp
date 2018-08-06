/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
//
// Created by Eric Irrgang on 8/27/17.
//

#include "tpxstate.h"

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
    std::make_shared<gmx_mtop_t>()
},
initialized_ {
    false
},
dirty_ {
    false
}
{
}

TpxState::~TpxState()
{

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
    return mtop_.get();
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
