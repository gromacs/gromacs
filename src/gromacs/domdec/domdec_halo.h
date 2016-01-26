/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 *  \brief
 * Halo communication for the force calculation.
 *
 * This file contains functions to set up the halo communication
 * for the eighth shell domain decomposition and to execute
 * the halo communication.
 * Only nearest neighbor communication is supported, which means that
 * each rank/domain sends and receives from 1, 3 and 7 ranks with 1D, 2D
 * and 3D decomposition respectively. All communication is non-blocking.
 *
 * \inlibraryapi
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_HALO_H
#define GMX_DOMDEC_DOMDEC_HALO_H

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"

struct gmx_ddbox_t;

/*! \brief Set up the eighth-shell halo coordinate/force communcation
 *
 * Set up the eighth-shell halo communcation for non-bonded + bonded
 * interactions
 * \p bCellsChanged indicates if the domain decomposition cells changed,
 * either due to dynamic load balancing or pressure scaling.
 *
 * \param[in,out] dd     Pointer to the domain decomposition setup.
 * \param[in]     box    Unit-cell matrix.
 * \param[in]     ddbox  Domain decomposition unit-cell and PBC data.
 * \param[in]     fr     Pointr to the force record struct.
 * \param[in]     bCellsChanged Tells if the domain decomposition cells changed.
 */
void setup_halo_communication(gmx_domdec_t *dd,
                              const matrix box, const gmx_ddbox_t *ddbox,
                              t_forcerec *fr,
                              gmx_bool bCellsChanged);

#endif
