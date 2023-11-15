/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2008- The GROMACS Authors
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
 * \brief This file declares functions for (mostly) the domdec module
 * to use MPI functionality
 *
 * \todo Wrap the raw dd_bcast in md.cpp into a higher-level function
 * in the domdec module, then this file can be module-internal.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_NETWORK_H
#define GMX_DOMDEC_DOMDEC_NETWORK_H

#include "gromacs/math/vectypes.h"

struct gmx_domdec_t;

namespace gmx
{
template<typename>
class ArrayRef;
}
/* \brief */
enum
{
    dddirForward,
    dddirBackward
};

/*! \brief Move a view of T values in the communication region one
 * cell along the domain decomposition
 *
 * Moves in the dimension indexed by ddDimensionIndex, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
template<typename T>
void ddSendrecv(const gmx_domdec_t* dd,
                int                 ddDimensionIndex,
                int                 direction,
                gmx::ArrayRef<T>    sendBuffer,
                gmx::ArrayRef<T>    receiveBuffer);

//! Extern declaration for int specialization
extern template void ddSendrecv<int>(const gmx_domdec_t* dd,
                                     int                 ddDimensionIndex,
                                     int                 direction,
                                     gmx::ArrayRef<int>  sendBuffer,
                                     gmx::ArrayRef<int>  receiveBuffer);

//! Extern declaration for real specialization
extern template void ddSendrecv<real>(const gmx_domdec_t* dd,
                                      int                 ddDimensionIndex,
                                      int                 direction,
                                      gmx::ArrayRef<real> sendBuffer,
                                      gmx::ArrayRef<real> receiveBuffer);

//! Extern declaration for gmx::RVec specialization
extern template void ddSendrecv<gmx::RVec>(const gmx_domdec_t*      dd,
                                           int                      ddDimensionIndex,
                                           int                      direction,
                                           gmx::ArrayRef<gmx::RVec> sendBuffer,
                                           gmx::ArrayRef<gmx::RVec> receiveBuffer);

/*! \brief Move revc's in the comm. region one cell along the domain decomposition
 *
 * Moves in dimension indexed by ddimind, simultaneously in the forward
 * and backward directions.
 */
void dd_sendrecv2_rvec(const struct gmx_domdec_t* dd,
                       int                        ddimind,
                       rvec*                      buf_s_fw,
                       int                        n_s_fw,
                       rvec*                      buf_r_fw,
                       int                        n_r_fw,
                       rvec*                      buf_s_bw,
                       int                        n_s_bw,
                       rvec*                      buf_r_bw,
                       int                        n_r_bw);


/* The functions below perform the same operations as the MPI functions
 * with the same name appendices, but over the domain decomposition
 * nodes only.
 * The DD main node is the coordinator for these operations.
 */

/*! \brief Broadcasts \p nbytes from \p data on \p DDMAINRANK to all PP ranks */
void dd_bcast(const gmx_domdec_t* dd, int nbytes, void* data);

/*! \brief Scatters \p nbytes from \p src on \p DDMAINRANK to all PP ranks, received in \p dest */
void dd_scatter(const gmx_domdec_t* dd, int nbytes, const void* src, void* dest);

/*! \brief Gathers \p nbytes from \p src on all PP ranks, received in \p dest on \p DDMAINRANK */
void dd_gather(const gmx_domdec_t* dd, int nbytes, const void* src, void* dest);

/*! \brief Scatters \p scounts elements of type \p T from \p src on \p DDMAINRANK to all PP ranks, receiving \p rcount * \p sc_scattervSize elements in \p dest.
 *
 * \tparam T  Data type, can only be int or real
 *
 * See man MPI_Scatterv for details of how to construct scounts and disps.
 * If rcount==0, rbuf is allowed to be nullptr.
 */
template<typename T>
void dd_scatterv(const gmx_domdec_t*      dd,
                 gmx::ArrayRef<const int> scounts,
                 gmx::ArrayRef<const int> disps,
                 const T*                 sbuf,
                 int                      rcount,
                 T*                       rbuf);

//! Instantiation of \p dd_scatterv for type int
extern template void dd_scatterv(const gmx_domdec_t*      dd,
                                 gmx::ArrayRef<const int> scounts,
                                 gmx::ArrayRef<const int> disps,
                                 const int*               sbuf,
                                 int                      rcount,
                                 int*                     rbuf);

//! Instantiation of \p dd_scatterv for type real
extern template void dd_scatterv(const gmx_domdec_t*      dd,
                                 gmx::ArrayRef<const int> scounts,
                                 gmx::ArrayRef<const int> disps,
                                 const gmx::RVec*         sbuf,
                                 int                      rcount,
                                 gmx::RVec*               rbuf);

/*! \brief Gathers \p rcount elements of type \p T from \p src on all PP ranks, received in \p scounts elements in \p dest on \p DDMAINRANK.
 *
 * \tparam T  Data type, can only be int or real
 *
 * See man MPI_Gatherv for details of how to construct rcounts and disps.
 *
 * sendBuffer is allowed to be empty.
 */
template<typename T>
void dd_gatherv(const gmx_domdec_t&      dd,
                gmx::ArrayRef<const T>   sendBuffer,
                gmx::ArrayRef<const int> rcounts,
                gmx::ArrayRef<const int> disps,
                gmx::ArrayRef<T>         receiveBuffer);

//! Instantiation of \p dd_gatherv for type int
extern template void dd_gatherv(const gmx_domdec_t&      dd,
                                gmx::ArrayRef<const int> sendBuffer,
                                gmx::ArrayRef<const int> rcounts,
                                gmx::ArrayRef<const int> disps,
                                gmx::ArrayRef<int>       receiveBuffer);

//! Instantiation of \p dd_gatherv for type real
extern template void dd_gatherv(const gmx_domdec_t&            dd,
                                gmx::ArrayRef<const gmx::RVec> sendBuffer,
                                gmx::ArrayRef<const int>       rcounts,
                                gmx::ArrayRef<const int>       disps,
                                gmx::ArrayRef<gmx::RVec>       receiveBuffer);

#endif
