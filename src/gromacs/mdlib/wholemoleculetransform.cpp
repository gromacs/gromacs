/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include "gmxpre.h"

#include "wholemoleculetransform.h"

#include <algorithm>

#include "gromacs/domdec/ga2la.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

WholeMoleculeTransform::WholeMoleculeTransform(const gmx_mtop_t& mtop,
                                               const PbcType     pbcType,
                                               const bool        useAtomReordering) :
    pbcType_(pbcType)
{
    gmx_localtop_t localTop(mtop.ffparams);

    gmx_mtop_generate_local_top(mtop, &localTop, false);

    graph_ = mk_graph(localTop.idef, mtop.natoms);

    wholeMoleculeCoordinates_.resize(mtop.natoms);

    if (useAtomReordering)
    {
        // Store a copy of the global edges
        globalEdgeAtomBegin_       = graph_.edgeAtomBegin;
        graphGlobalAtomOrderEdges_ = graph_.edges;
        // Reset the graph range to cover the whole system
        graph_.edgeAtomBegin = 0;
        graph_.edgeAtomEnd   = graph_.shiftAtomEnd;

        // Resize the edge color list for potential addition of non-connected atoms
        graph_.edgeColor.resize(graph_.edgeAtomEnd);
    }
}

void WholeMoleculeTransform::updateAtomOrder(ArrayRef<const int> globalAtomIndices, const gmx_ga2la_t& ga2la)
{
    GMX_RELEASE_ASSERT(!graphGlobalAtomOrderEdges_.empty(),
                       "We need the edges in the global atom order");

    GMX_RELEASE_ASSERT(globalAtomIndices.ssize() == graph_.shiftAtomEnd,
                       "The number of indices should match the number of atoms we shift");

    const int globalEdgeAtomEnd = globalEdgeAtomBegin_ + graphGlobalAtomOrderEdges_.size();
    graph_.edges.clear();
    for (const int globalAtomIndex : globalAtomIndices)
    {
        if (globalAtomIndex >= globalEdgeAtomBegin_ && globalAtomIndex < globalEdgeAtomEnd)
        {
            // Get the list of global atoms linked to this atom and add their local indices
            const ArrayRef<const int> globalList =
                    graphGlobalAtomOrderEdges_[globalAtomIndex - globalEdgeAtomBegin_];
            graph_.edges.pushBackListOfSize(globalList.size());
            ArrayRef<int> lastList = graph_.edges.back();
            std::transform(globalList.begin(), globalList.end(), lastList.begin(), [&ga2la](int i) -> int {
                return ga2la.find(i)->la;
            });
        }
        else
        {
            // The atom has no edges, push back an empty list
            graph_.edges.pushBackListOfSize(0);
        }
    }

    GMX_RELEASE_ASSERT(int(graph_.edges.size()) == graph_.shiftAtomEnd,
                       "We should have as many lists of edges as the system (shift) size");
}

void WholeMoleculeTransform::updateForAtomPbcJumps(ArrayRef<const RVec> x, const matrix box)
{
    mk_mshift(nullptr, &graph_, pbcType_, box, as_rvec_array(x.data()));
}

ArrayRef<const RVec> WholeMoleculeTransform::wholeMoleculeCoordinates(ArrayRef<const RVec> x, const matrix box)
{
    shift_x(&graph_, box, as_rvec_array(x.data()), as_rvec_array(wholeMoleculeCoordinates_.data()));

    return wholeMoleculeCoordinates_;
}

} // namespace gmx
