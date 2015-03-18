/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015, by the GROMACS development team, led by
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
 * \brief This file declares functions for inter-rank signalling by mdrun
 *
 * This handles details of responding to termination conditions,
 * coordinating checkpoints, and coordinating multi-simulations.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_MDRUN_SIGNALLING_H
#define GMX_MDLIB_MDRUN_SIGNALLING_H

#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/utility/real.h"

struct t_commrec;

namespace gmx
{

template <typename T>
class ArrayRef;

}

/*! \brief Simulation conditions to transmit.
 *
 * Keep in mind that they are transmitted to other ranks through an
 * MPI_Reduce after casting them to a real (so the signals can be sent
 * together with other data). This means that the only meaningful
 * values are positive, negative or zero. */
enum {
    eglsCHKPT, eglsSTOPCOND, eglsRESETCOUNTERS, eglsNR
};

/*! \internal
 * \brief Object used by mdrun ranks to signal to each other
 *
 * Note that xlc on BG/Q requires sig to be of size char (see unit tests
 * of ArrayRef for details). */
struct gmx_signalling_t {
    int  nstms;             /**< The frequency for inter-simulation communication */
    char sig[eglsNR];       /**< The signal set by this rank in do_md */
    char set[eglsNR];       /**< The communicated signal, equal for all ranks once communication has occurred */
    real mpiBuffer[eglsNR]; /**< Buffer for communication */
};

/*! \brief Construct a struct gmx_signalling_t */
void init_global_signals(struct gmx_signalling_t *gs,
                         const t_commrec         *cr,
                         const t_inputrec        *ir,
                         int                      repl_ex_nst);

/*! \brief Fill the array of reals in which inter- and
 * intra-simulation signals will be communicated
 * with the signal values to be sent. */
gmx::ArrayRef<real>
prepareSignalBuffer(struct gmx_signalling_t *gs);

/*! \brief Handle intra- and inter-simulation signals recieved
 *
 * If a multi-simulation signal should be handled, communicate between
 * simulation-master ranks, then propagate from the masters to the
 * rest of the ranks for each simulation.
 *
 * Then, set the flags that mdrun will use to respond to the signals
 * received. */
void
handleSignals(struct gmx_signalling_t  *gs,
              const t_commrec          *cr,
              bool                      bInterSimGS);

#endif
