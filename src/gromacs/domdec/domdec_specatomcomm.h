/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2012,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

typedef struct {
    int  nsend;
    int *a;
    int  a_nalloc;
    int  nrecv;
} gmx_specatsend_t;

typedef struct {
    int *ind;
    int  nalloc;
    int  n;
} ind_req_t;

/*! \internal \brief Struct with setup and buffers for special atom communication */
typedef struct gmx_domdec_specat_comm {
    /* The number of indices to receive during the setup */
    int              nreq[DIM][2][2];  /**< The nr. of atoms requested, per DIM, direction and direct/indirect */
    /* The atoms to send */
    gmx_specatsend_t spas[DIM][2];     /**< The communication setup per DIM, direction */
    gmx_bool        *bSendAtom;        /**< Work buffer that tells if spec.atoms should be sent */
    int              bSendAtom_nalloc; /**< Allocation size of \p bSendAtom */
    /* Send buffers */
    int             *ibuf;             /**< Integer send buffer */
    int              ibuf_nalloc;      /**< Allocation size of \p ibuf */
    rvec            *vbuf;             /**< rvec send buffer */
    int              vbuf_nalloc;      /**< Allocation size of \p vbuf */
    rvec            *vbuf2;            /**< rvec send buffer */
    int              vbuf2_nalloc;     /**< Allocation size of \p vbuf2 */
    /* The range in the local buffer(s) for received atoms */
    int              at_start;         /**< Start index of received atoms */
    int              at_end;           /**< End index of received atoms */

    /* The atom indices we need from the surrounding cells.
     * We can gather the indices over nthread threads.
     */
    int        nthread;                /**< Number of threads used for spec.atom communication */
    ind_req_t *ireq;                   /**< Index request buffer per thread, allocation size \p nthread */
} gmx_domdec_specat_comm_t;

/*! \brief Communicates the force for special atoms, the shift forces are reduced with \p fshift != NULL */
void dd_move_f_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                      rvec *f, rvec *fshift);

/*! \brief Communicates the coordinates for special atoms
 *
 * \param[in]     dd         Domain decomposition struct
 * \param[in]     spac       Special atom communication struct
 * \param[in]     box        Box, used for pbc
 * \param[in,out] x0         Vector to communicate
 * \param[in,out] x1         Vector to communicate, when != NULL
 * \param[in]     bX1IsCoord Tells is \p x1 is a coordinate vector (needs pbc)
 */
void dd_move_x_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                      matrix box,
                      rvec *x0,
                      rvec *x1, gmx_bool bX1IsCoord);

/*! \brief Sets up the communication for special atoms
 *
 * \param[in]     dd         Domain decomposition struct
 * \param[in]     ireq List of requested atom indices
 * \param[in,out] spac   Special atom communication struct
 * \param[out]    ga2la_specat Global to local special atom index
 * \param[in]     at_start     Index in local state where to start storing communicated atoms
 * \param[in]     vbuf_fac     Buffer factor, 1 or 2 for communicating 1 or 2 vectors
 * \param[in]     specat_type  Name of the special atom, used for error message
 * \param[in]     add_err      Text to add at the end of error message when atoms can't be found
 */
int setup_specat_communication(gmx_domdec_t             *dd,
                               ind_req_t                *ireq,
                               gmx_domdec_specat_comm_t *spac,
                               gmx_hash_t                ga2la_specat,
                               int                       at_start,
                               int                       vbuf_fac,
                               const char               *specat_type,
                               const char               *add_err);

/*! \brief Initialize a special communication struct */
gmx_domdec_specat_comm_t *specat_comm_init(int nthread);

#endif
