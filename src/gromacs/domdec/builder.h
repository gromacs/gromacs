/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief This file declares a builder class for the manager
 * of domain decomposition
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_BUILDER_H
#define GMX_DOMDEC_BUILDER_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_localtop_t;
struct t_commrec;
struct t_inputrec;
class t_state;

namespace gmx
{
class MDLogger;
class LocalAtomSetManager;
class RangePartitioning;
struct DomdecOptions;
struct MdrunOptions;
struct MDModulesNotifiers;
class ObservablesReducerBuilder;

template<typename T>
class ArrayRef;

/*! \libinternal
 * \brief Builds a domain decomposition management object
 *
 * This multi-phase construction needs first a decision about the
 * duty(s) of each rank, and then perhaps to be advised of GPU streams
 * for transfer operations. */
class DomainDecompositionBuilder
{
public:
    //! Constructor
    DomainDecompositionBuilder(const MDLogger&                   mdlog,
                               t_commrec*                        cr,
                               const DomdecOptions&              options,
                               const MdrunOptions&               mdrunOptions,
                               const gmx_mtop_t&                 mtop,
                               const t_inputrec&                 ir,
                               const MDModulesNotifiers&         notifiers,
                               const matrix                      box,
                               ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
                               bool                              useUpdateGroups,
                               real                              maxUpdateGroupRadius,
                               ArrayRef<const RVec>              xGlobal,
                               bool                              useGpuForNonbonded,
                               bool                              useGpuForPme,
                               bool                              useGpuForUpdate,
                               bool                              useGpuDirectHalo,
                               bool                              canUseGpuPmeDecomposition);
    //! Destructor
    ~DomainDecompositionBuilder();
    //! Build the resulting DD manager
    std::unique_ptr<gmx_domdec_t> build(LocalAtomSetManager*       atomSets,
                                        const gmx_localtop_t&      localTopology,
                                        const t_state&             localState,
                                        ObservablesReducerBuilder* observablesReducerBuilder);

private:
    class Impl;
    //! Pimpl to hide implementation details
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
