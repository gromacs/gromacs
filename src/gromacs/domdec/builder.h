/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;

namespace gmx
{
class MDLogger;
class LocalAtomSetManager;
struct DomdecOptions;
struct MdrunOptions;

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
    DomainDecompositionBuilder(const MDLogger&      mdlog,
                               t_commrec*           cr,
                               const DomdecOptions& options,
                               const MdrunOptions&  mdrunOptions,
                               bool                 prefer1DAnd1Pulse,
                               const gmx_mtop_t&    mtop,
                               const t_inputrec&    ir,
                               const matrix         box,
                               ArrayRef<const RVec> xGlobal);
    //! Destructor
    ~DomainDecompositionBuilder();
    //! Build the resulting DD manager
    gmx_domdec_t* build(LocalAtomSetManager* atomSets);

private:
    class Impl;
    //! Pimpl to hide implementation details
    PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
