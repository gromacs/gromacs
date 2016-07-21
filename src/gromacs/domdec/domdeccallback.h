/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::DomDecCallBack
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDECCALLBACK_H
#define GMX_DOMDECCALLBACK_H

#include <memory>
#include <vector>

#include "gromacs/utility/classhelpers.h"

struct gmx_ga2la_t;
struct gmx_localtop_t;
struct t_mdatoms;

namespace gmx
{

/*! \libinternal \brief Abstract base class for all domain decomposition callbacks.
 *
 * \ingroup module_domdec */
class DomDecCallBack
{
    friend class DomDecCallBackContainer;
    private:
        /*! \brief Triggered from DomDecCallBackContainer after domain decomposition routine is finished.
         *
         * \param[in] ga2la global atom to local atom lookup
         * \param[in] top_local the node-local topology
         * \param[in] mdatoms the md atoms
         */
        virtual void domDecDone(const gmx_ga2la_t *ga2la, const gmx_localtop_t *top_local, const t_mdatoms *mdatoms) = 0;
};

/*! \libinternal \brief
 * Hold handles to all functions that shall be triggered during domain decomposition.
 *
 * \ingroup module_domdec
 */
class DomDecCallBackContainer
{
    public:

        typedef std::unique_ptr<DomDecCallBack> DomDecCallBackHandle;
        DomDecCallBackContainer();

        /*! \brief Triggered after domain decomposition routine is finished.
         *
         * \param[in] ga2la look-up to identify local atoms
         * \param[in] top_local local topology data
         * \param[in] mdatoms local md atom data
         * */
        void triggerCallBackDomDecDone(const gmx_ga2la_t *ga2la, const gmx_localtop_t *top_local, const t_mdatoms *mdatoms);

        /*! \brief The DomDecCallBackContainer handles all callback to domain decomposition during a simulation.
         * \param[in] callbackhandle Handle to the callback.
         */
        void add(DomDecCallBackHandle &&callbackhandle);

    private:
        class Impl;
        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
