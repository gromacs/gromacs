/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Declares the topology class for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */


#ifndef GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H
#define GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H

#include <vector>

#include "modularsimulatorinterfaces.h"

struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_vsite_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;

namespace gmx
{
class Constraints;
class MDAtoms;

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Object holding the topology
 *
 * Clients can register to get an updated local topology whenever there
 * is a change (infrequent, only due to domdec currently).
 */
class TopologyHolder final
{
public:
    //! Constructor
    TopologyHolder(const gmx_mtop_t& globalTopology,
                   const t_commrec*  cr,
                   const t_inputrec* inputrec,
                   t_forcerec*       fr,
                   MDAtoms*          mdAtoms,
                   Constraints*      constr,
                   gmx_vsite_t*      vsite);

    //! Get global topology
    const gmx_mtop_t& globalTopology() const;

    //! Register topology client
    void registerClient(ITopologyHolderClient* client);

    //! Allow domdec to update local topology
    friend class DomDecHelper;

private:
    //! Constant reference to the global topolgy
    const gmx_mtop_t& globalTopology_;
    //! Pointer to the currently valid local topology
    std::unique_ptr<gmx_localtop_t> localTopology_;

    //! List of clients to be updated if local topology changes
    std::vector<ITopologyHolderClient*> clients_;

    //! Update local topology
    void updateLocalTopology();
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_TOPOLOGYHOLDER_H
