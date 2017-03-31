/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/// \cond
// Adding this file to the doxygen mesh isn't essential right now.
/*! \libinternal \file
 * \ingroup module_trajectory Classes for handling trajectory data
 * \brief defines t_trxframe frame data and operations
 */

/* The gmx_bools indicate whether a field was read from the trajectory.
 * Do not try to use a pointer when its gmx_bool is FALSE, as memory might
 * not be allocated.
 */
#ifndef GMX_TRAJECTORY_TRX_H
#define GMX_TRAJECTORY_TRX_H

#include <cstdio>

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
// required to be able to delete t_atoms
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/smalloc.h"

/*! \libinternal \brief
 * C structure of basic data read from a trajectory
 *
 * \ingroup module_trajectory
 * \ingroup libraryapi
 * \{
 */
struct t_trxframe
{
    int             not_ok;    ///< integrity flags
    gmx_bool        bDouble;   ///< Double precision?
    int             natoms;    ///< number of atoms (atoms, x, v, f, index) */
    gmx_bool        bTitle;    ///< has title
    const char     *title;     ///< title of the frame
    gmx_bool        bStep;     ///< has step
    gmx_int64_t     step;      ///< MD step number
    gmx_bool        bTime;     ///< has time
    real            time;      ///< time of the frame
    gmx_bool        bLambda;   ///< has lambda
    gmx_bool        bFepState; ///< does it contain fep_state?
    real            lambda;    ///< free energy perturbation lambda
    int             fep_state; ///< which fep state are we in?
    gmx_bool        bAtoms;    ///< has atoms
    t_atoms        *atoms;     ///< atoms struct (natoms)
    gmx_bool        bPrec;     ///< has prec
    real            prec;      ///< precision of x, fraction of 1 nm
    gmx_bool        bX;        ///< has x
    rvec           *x;         ///< coordinates (natoms)
    gmx_bool        bV;        ///< has v
    rvec           *v;         ///< velocities (natoms)
    gmx_bool        bF;        ///< has f
    rvec           *f;         ///< forces (natoms)
    gmx_bool        bBox;      ///< has box
    matrix          box;       ///< the 3 box vectors
    gmx_bool        bPBC;      ///< has pbc
    int             ePBC;      ///< the type of pbc
    gmx_bool        bIndex;    ///< has index
    int            *index;     ///< atom indices of contained coordinates
};
typedef struct t_trxframe t_trxframe;
/// \}

// Request explicit instantiation.
template class std::shared_ptr<t_trxframe>;
template class std::unique_ptr<t_trxframe>;

/*! \libinternal \brief Compare trajectory frames
 */
void comp_frame(FILE *fp, t_trxframe *fr1, t_trxframe *fr2,
                gmx_bool bRMSD, real ftol, real abstol);

/*! \libinternal \brief Frees memory for a trajectory frame
 *
 */
void done_frame(t_trxframe *frame);

namespace gmx
{
namespace trajectory
{

/*! \brief Deleter for t_trxframe
 *
 * Used for objects created with frame_copy. Not sure if this tackles the same
 * exact task as done_frame, introduced recently.
 * \internal
 */
void trxframe_deleter(t_trxframe* f);

/*! \brief make a deep copy of a trajectory frame
 *
 * Allocates heap memory for structure and member arrays.
 * Returns nullptr if unable to allocate memory.
 * A lambda deleter would be nice, but linkage gets tricky.
 * \param frame trx frame to copy
 * \returns trxframe_ptr
 * \internal
 */
std::unique_ptr<t_trxframe, void(*)(t_trxframe*)> trxframe_copy(const t_trxframe &frame);

} // end namespace trajectory
} //end namespace gmx

#endif
/// \endcond
