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

#ifndef GMX_MDTYPES_TPXSTATE_H
#define GMX_MDTYPES_TPXSTATE_H

#include <atomic>
#include <memory>
#include <mutex>
#include <string>

struct t_inputrec;
class t_state;
struct gmx_mtop_t;

namespace gmx
{

/*!
 * \brief Container for microstate data.
 *
 * This is a dead end and needs to be removed.
 *
 * This class was intended as the type to be returned by read_tpx_state() and to
 * hold sufficient information to provide the starting point (or continuation
 * point) for a simulation trajectory. TpxState objects may (but are not
 * required to) own a t_state global state, a gmx_mtop_t molecular topology, and
 * a t_inputrec input record instance. Ownership of a TpxState may be shared by
 * a gmx::Mdrunner and the code owning the gmx::Mdrunner.
 */
class TpxState final
{
    private:
        // TpxState is currently always file-backed.
        std::string filename_;
        std::shared_ptr<t_inputrec>                     inputrecInstance_;
        std::shared_ptr<t_state>                        stateInstance_;
        std::shared_ptr<gmx_mtop_t>                     mtop_;
        std::atomic<bool>                               initialized_;
        std::atomic<bool>                               dirty_;
        mutable std::mutex                              exclusive_;
    public:
        TpxState();
        ~TpxState();

        // Copy semantics TBD. In addition to unclear copy semantics of members, probably need to use setters to allow
        // for notifications of data changs.
        TpxState(const TpxState&) = delete;
        TpxState               &operator=(const TpxState &) = delete;

        // Move should be okay
        TpxState(TpxState && old) noexcept;
        TpxState &operator=(TpxState &&) noexcept;

        static std::unique_ptr<TpxState> initializeFromFile(const char* filename);
        static std::unique_ptr<TpxState> initializeFromFile(const std::string &filename);

        // Takes ownership of arguments to be members of new object.
        static std::unique_ptr<TpxState>
        initializeFromWrappers(std::unique_ptr<t_inputrec> inputRecord,
                               std::unique_ptr<t_state>    state,
                               std::unique_ptr<gmx_mtop_t> mtop);

        t_inputrec* getRawInputrec();
        gmx_mtop_t* getRawMtop();
        t_state*    getRawState();

        /// \returns Whether data has been loaded into the object
        bool isInitialized() const;

        /// \returns true if we do not have a guarantee that the object is in a self-consistent state
        bool isDirty() const;

        /// \brief Allow caller to assert the validity of an instance.
        void markClean();

        const char* filename() const;
};


}      // end namespace gmx

#endif //GROMACS_TPXSTATE_H
