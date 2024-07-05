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
 * \brief This file declares functions used by the domdec module
 * for (bounding) box and pbc information generation.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_BOX_H
#define GMX_DOMDEC_BOX_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{
template<typename>
class ArrayRef;
}
struct gmx_ddbox_t;
struct gmx_domdec_t;
struct t_commrec;
struct t_inputrec;
enum class DDRole;

/*! \brief Set the box and PBC data in \p ddbox */
void set_ddbox(const gmx_domdec_t&            dd,
               bool                           mainRankHasTheSystemState,
               const matrix                   box,
               bool                           calculateUnboundedSize,
               gmx::ArrayRef<const gmx::RVec> x,
               gmx_ddbox_t*                   ddbox);

/*! \brief Set the box and PBC data in \p ddbox */
void set_ddbox_cr(DDRole                         ddRole,
                  MPI_Comm                       communicator,
                  const gmx::IVec*               numDomains,
                  const t_inputrec&              ir,
                  const matrix                   box,
                  gmx::ArrayRef<const gmx::RVec> x,
                  gmx_ddbox_t*                   ddbox);

/*! \brief Computes and returns a domain decomposition box */
gmx_ddbox_t get_ddbox(const gmx::IVec&               numDomains,
                      const t_inputrec&              ir,
                      const matrix                   box,
                      gmx::ArrayRef<const gmx::RVec> x);

#endif
