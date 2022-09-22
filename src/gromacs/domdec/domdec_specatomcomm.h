/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2005- The GROMACS Authors
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
 * \brief This file declares functions for domdec to use
 * while managing communication of atoms required for special purposes
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_SPECATOMCOMM_H
#define GMX_DOMDEC_DOMDEC_SPECATOMCOMM_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_domdec_t;
struct t_commrec;

namespace gmx
{
template<typename T>
class HashedMap;
} // namespace gmx

/*! \internal \brief The communication setup along a single dimension */
struct gmx_specatsend_t
{
    std::vector<int> a;     /**< The indices of atoms to send */
    int              nrecv; /**< The number of atoms to receive */
};

/*! \internal \brief Struct with setup and buffers for special atom communication */
struct gmx_domdec_specat_comm_t
{
    /* The number of indices to receive during the setup */
    int nreq[DIM][2][2] = { { { 0 } } }; /**< The nr. of atoms requested, per DIM, direction and direct/indirect */
    /* The atoms to send */
    gmx_specatsend_t  spas[DIM][2]; /**< The communication setup per DIM, direction */
    std::vector<bool> sendAtom;     /**< Work buffer that tells if spec.atoms should be sent */

    /* Send buffers */
    std::vector<int>       ibuf;  /**< Integer send buffer */
    std::vector<gmx::RVec> vbuf;  /**< rvec send buffer */
    std::vector<gmx::RVec> vbuf2; /**< rvec send buffer */
    /* The range in the local buffer(s) for received atoms */
    int at_start; /**< Start index of received atoms */
    int at_end;   /**< End index of received atoms */
};

/*! \brief Communicates the force for special atoms, the shift forces are reduced with \p fshift != NULL */
void dd_move_f_specat(const gmx_domdec_t* dd, gmx_domdec_specat_comm_t* spac, gmx::RVec* f, gmx::RVec* fshift);

/*! \brief Communicates the coordinates for special atoms
 *
 * \param[in]     dd         Domain decomposition struct
 * \param[in]     spac       Special atom communication struct
 * \param[in]     box        Box, used for pbc
 * \param[in,out] x0         Vector to communicate
 * \param[in,out] x1         Vector to communicate, when != NULL
 * \param[in]     bX1IsCoord Tells is \p x1 is a coordinate vector (needs pbc)
 */
void dd_move_x_specat(const gmx_domdec_t*       dd,
                      gmx_domdec_specat_comm_t* spac,
                      const matrix              box,
                      gmx::RVec*                x0,
                      gmx::RVec*                x1,
                      bool                      bX1IsCoord);

/*! \brief Sets up the communication for special atoms
 *
 * \param[in]     dd           Domain decomposition struct
 * \param[in,out] ireq         List of requested atom indices, updated due to aggregation
 * \param[in,out] spac         Special atom communication struct
 * \param[out]    ga2la_specat Global to local special atom index
 * \param[in]     at_start     Index in local state where to start storing communicated atoms
 * \param[in]     vbuf_fac     Buffer factor, 1 or 2 for communicating 1 or 2 vectors
 * \param[in]     specat_type  Name of the special atom, used for error message
 * \param[in]     add_err      Text to add at the end of error message when atoms can't be found
 */
int setup_specat_communication(gmx_domdec_t*             dd,
                               std::vector<int>*         ireq,
                               gmx_domdec_specat_comm_t* spac,
                               gmx::HashedMap<int>*      ga2la_specat,
                               int                       at_start,
                               int                       vbuf_fac,
                               const char*               specat_type,
                               const char*               add_err);

#endif
