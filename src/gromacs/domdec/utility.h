/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief Declares utility functions used in the domain decomposition module.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_UTILITY_H
#define GMX_DOMDEC_DOMDEC_UTILITY_H

#include <cstddef>

#include <array>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "domdec_internal.h"

class t_state;
struct t_forcerec;

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \brief Returns true if the DLB state indicates that the balancer is on. */
static inline bool isDlbOn(const DlbState& dlbState)
{
    return (dlbState == DlbState::onCanTurnOff || dlbState == DlbState::onUser);
};

/*! \brief Returns true if the DLB state indicates that the balancer is off/disabled.
 */
static inline bool isDlbDisabled(const DlbState& dlbState)
{
    return (dlbState == DlbState::offUser || dlbState == DlbState::offForever);
};

/*! \brief Returns the character, x/y/z, corresponding to dimension dim */
char dim2char(int dim);

/*! \brief Sets matrix to convert from Cartesian to lattice coordinates */
void make_tric_corr_matrix(int npbcdim, const matrix box, matrix tcm);

/*! \brief Ensure box obeys the screw restrictions, fatal error if not */
void check_screw_box(const matrix box);

/*! \brief Return the atom information flags for atom a */
static inline int ddGetAtomInfo(gmx::ArrayRef<const gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock,
                                int                                                   a)
{
    size_t index = 0;
    while (a >= atomInfoForEachMoleculeBlock[index].indexOfLastAtomInMoleculeBlock)
    {
        index++;
    }
    const gmx::AtomInfoWithinMoleculeBlock& atomInfoOfMoleculeBlock = atomInfoForEachMoleculeBlock[index];

    return atomInfoOfMoleculeBlock.atomInfo[(a - atomInfoOfMoleculeBlock.indexOfFirstAtomInMoleculeBlock)
                                            % atomInfoOfMoleculeBlock.atomInfo.size()];
};

/*! \brief Returns the number of MD steps for which load has been recorded */
static inline int dd_load_count(const gmx_domdec_comm_t* comm)
{
    return (comm->ddSettings.eFlop ? comm->flop_n : comm->cycl_n[ddCyclF]);
}

/*! \brief Ensure fr and state can hold numAtoms atoms
 *
 * \param[in]  fr        Force record
 * \param[in]  state     Current state
 * \param[out] numAtoms  Number of atoms
 */
void dd_resize_atominfo_and_state(t_forcerec* fr, t_state* state, int numAtoms);

/*! \brief Returns a domain-to-domain cutoff distance given an atom-to-atom cutoff */
static inline real atomToAtomIntoDomainToDomainCutoff(const DDSystemInfo& systemInfo, real cutoff)
{
    if (systemInfo.useUpdateGroups)
    {
        cutoff += 2 * systemInfo.maxUpdateGroupRadius;
    }

    return cutoff;
}

/*! \brief Returns an atom-to-domain cutoff distance given a domain-to-domain cutoff */
static inline real domainToDomainIntoAtomToDomainCutoff(const DDSystemInfo& systemInfo, real cutoff)
{
    if (systemInfo.useUpdateGroups)
    {
        cutoff -= systemInfo.maxUpdateGroupRadius;
    }

    return cutoff;
}

#endif
