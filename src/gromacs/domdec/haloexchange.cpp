/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Implements the HaloExchange class for halo communication of coordinates and forces
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "haloexchange.h"

#include "config.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/template_mp.h"

#include "domainpaircomm.h"
#include "domdec_internal.h"

namespace gmx
{

HaloExchange::HaloExchange(const PbcType pbcType) : pbcType_(pbcType) {}

HaloExchange::~HaloExchange() = default;

namespace
{

/*! \brief Initiates a non-blocking send to another domain
 *
 * \param[in] send         The domain pair communication setup
 * \param[in] sendBuffer   The data to send
 * \param[in] tag          The MPI tag
 * \param[in] mpiRequests  List of requests where the send will be appended to
 */
template<typename DomainPairComm, typename T>
void ddIsendDomain(const DomainPairComm&     send,
                   ArrayRef<T>               sendBuffer,
                   const HaloMpiTag          tag,
                   std::vector<MPI_Request>* mpiRequests)
{
#if GMX_MPI
    if (send.numAtoms() > 0)
    {
        mpiRequests->emplace_back();

        MPI_Isend(sendBuffer.data(),
                  send.numAtoms() * sizeof(T),
                  MPI_BYTE,
                  send.rank(),
                  static_cast<int>(tag),
                  send.mpiCommAll(),
                  &mpiRequests->back());
    }
#else
    GMX_UNUSED_VALUE(send);
    GMX_UNUSED_VALUE(sendBuffer);
    GMX_UNUSED_VALUE(tag);
    GMX_UNUSED_VALUE(mpiRequests);
#endif
}

/*! \brief Initiates a non-blocking receive from another domain
 *
 * \param[in] receive        The domain pair communication setup
 * \param[in] receiveBuffer  Buffer to receive the data in
 * \param[in] tag            The MPI tag
 * \param[in] mpiRequests    List of requests where the send will be appended to
 */
template<typename DomainPairComm, typename T>
void ddIreceiveDomain(const DomainPairComm&     receive,
                      ArrayRef<T>               receiveBuffer,
                      const HaloMpiTag          tag,
                      std::vector<MPI_Request>* mpiRequests)
{
    GMX_ASSERT(receiveBuffer.ssize() >= receive.numAtoms(),
               "Receive buffer should be sufficiently large");

#if GMX_MPI
    if (receive.numAtoms() > 0)
    {
        mpiRequests->emplace_back();

        MPI_Irecv(receiveBuffer.data(),
                  receive.numAtoms() * sizeof(T),
                  MPI_BYTE,
                  receive.rank(),
                  static_cast<int>(tag),
                  receive.mpiCommAll(),
                  &mpiRequests->back());
    }
#else
    GMX_UNUSED_VALUE(receive);
    GMX_UNUSED_VALUE(receiveBuffer);
    GMX_UNUSED_VALUE(tag);
    GMX_UNUSED_VALUE(mpiRequests);
#endif
}

//! Templated version of \p packCoordinateSendBuffer()
template<bool commOverPbc, bool usesScrewPbc>
void packCoordinatesTemplated(const DomainCommBackward& domainComm,
                              const matrix              box,
                              const RVec&               shiftVec,
                              ArrayRef<const RVec>      x,
                              ArrayRef<RVec>            sendBuffer)
{
    const int numAtomsPerCell = domainComm.numAtomsPerCell();

    int j = 0;
    for (const DomainCommBackward::ColumnInfo& columnInfo : domainComm.columnsToSend())
    {
        for (int cell : columnInfo.cellRange)
        {
            for (int i = 0; i < numAtomsPerCell; i++)
            {
                sendBuffer[j] = x[cell * numAtomsPerCell + i];

                if constexpr (usesScrewPbc)
                {
                    sendBuffer[j][XX] += shiftVec[XX];
                    /* Rotate y and z.
                     * This operation requires a special shift force
                     * treatment, which is performed in calc_vir.
                     */
                    sendBuffer[j][YY] = box[YY][YY] - sendBuffer[j][YY];
                    sendBuffer[j][ZZ] = box[YY][ZZ] - sendBuffer[j][ZZ];
                }
                else if constexpr (commOverPbc)
                {
                    sendBuffer[j] += shiftVec;
                }

                j++;
            }
        }
    }

    GMX_ASSERT(j == domainComm.numAtoms(), "We should have handled all atoms to send");
}

//! Wrapper for MPI_Waitall that takes ArrayRefs
void mpiWaitall(ArrayRef<MPI_Request> mpiRequests, ArrayRef<MPI_Status> mpiStatuses)
{
#if GMX_MPI
    GMX_ASSERT(mpiStatuses.size() >= mpiRequests.size(), "We need sufficients statuses");

    int ret = MPI_Waitall(gmx::ssize(mpiRequests), mpiRequests.data(), mpiStatuses.data());

    GMX_RELEASE_ASSERT(ret == MPI_SUCCESS, "MPI_Waitall failed");

    GMX_UNUSED_VALUE(ret);
#else  // GMX_MPI
    GMX_UNUSED_VALUE(mpiRequests);
    GMX_UNUSED_VALUE(mpiStatuses);
#endif // GMX_MPI
}

} // namespace

void DomainCommBackward::packCoordinateSendBuffer(const matrix         box,
                                                  ArrayRef<const RVec> x,
                                                  ArrayRef<RVec>       sendBuffer) const
{
    RVec shiftVec = { 0.0_real, 0.0_real, 0.0_real };

    if (commOverPbc_)
    {
        for (int d = 0; d < DIM; d++)
        {
            for (int e = 0; e < DIM; e++)
            {
                shiftVec[e] += pbcCoordinateShift_[d] * box[d][e];
            }
        }
    }

    dispatchTemplatedFunction(
            [&](auto commOverPbc, auto usesScrewPbc) {
                packCoordinatesTemplated<commOverPbc, usesScrewPbc>(*this, box, shiftVec, x, sendBuffer);
            },
            commOverPbc_,
            usesScrewPbc_);
}

void HaloExchange::initiateReceiveX(ArrayRef<RVec> x)
{
    auto& mpiRequests = mpiCoordinateRequests_.receive;

    /* Post all the non-blocking receives */
    for (DomainPairComm& dpc : domainPairComm_)
    {
        DomainCommForward& receive = dpc.forward();

        if (receive.numAtoms() > 0)
        {
            ddIreceiveDomain(receive,
                             x.subArray(*receive.atomRange().begin(), receive.atomRange().size()),
                             HaloMpiTag::X,
                             &mpiRequests);
        }
    }
}

void HaloExchange::initiateSendX(const matrix box, ArrayRef<RVec> x)
{
    auto& mpiRequests = mpiCoordinateRequests_.send;

    GMX_ASSERT(mpiRequests.empty(),
               "All MPI Requests should have been handled before initiating sendX");

    for (DomainPairComm& dpc : domainPairComm_)
    {
        DomainCommBackward& send = dpc.backward();

        if (send.numAtoms() > 0)
        {
            send.packCoordinateSendBuffer(box, x, send.rvecBuffer());

            // Post the non-blocking send
            ddIsendDomain(send, send.rvecBuffer(), HaloMpiTag::X, &mpiRequests);
        }
    }
}

void HaloExchange::completeReceiveX()
{
    HaloMpiRequests& mpiRequests = mpiCoordinateRequests_;

    /* Wait for all non-blocking receives to complete */
    mpiWaitall(mpiRequests.receive, mpiStatus_);

    mpiRequests.receive.clear();
}

void HaloExchange::completeSendX()
{
    HaloMpiRequests& mpiRequests = mpiCoordinateRequests_;

    /* Wait for all non-blocking sends to complete */
    mpiWaitall(mpiRequests.send, mpiStatus_);

    mpiRequests.send.clear();
}

void HaloExchange::moveX(const matrix box, ArrayRef<RVec> x)
{
    initiateReceiveX(x);
    initiateSendX(box, x);
    completeReceiveX();
    completeSendX();
}

namespace
{

//! Templated version of \p packCoordinateSendBuffer()
template<bool usesScrewPbc, bool haveShiftForces>
void accumulateReceivedForcesTemplated(const DomainCommBackward& domainComm,
                                       ArrayRef<RVec>            forces,
                                       RVec gmx_unused*          shiftForce)
{
    const int numAtomsPerCell = domainComm.numAtomsPerCell();

    ArrayRef<const RVec> receivedForces = domainComm.rvecBuffer();

    int j = 0;
    for (const auto& columnInfo : domainComm.columnsToSend())
    {
        for (int cell : columnInfo.cellRange)
        {
            for (int i = 0; i < numAtomsPerCell; i++)
            {
                if constexpr (!usesScrewPbc)
                {
                    forces[cell * numAtomsPerCell + i] += receivedForces[j];
                }
                else
                {
                    // Accumulate the forces after rotating them
                    forces[cell * numAtomsPerCell + i][XX] += receivedForces[j][XX];
                    forces[cell * numAtomsPerCell + i][YY] -= receivedForces[j][YY];
                    forces[cell * numAtomsPerCell + i][ZZ] -= receivedForces[j][ZZ];
                }

                if constexpr (haveShiftForces)
                {
                    // Add this force to the shift force
                    *shiftForce += receivedForces[j];
                }

                j++;
            }
        }
    }

    GMX_ASSERT(j == domainComm.numAtoms(), "We should have handled all atoms to send");
}

} // namespace

void DomainCommBackward::accumulateReceivedForces(ArrayRef<RVec> forces, ArrayRef<RVec> shiftForces) const
{
    const bool haveShiftForces = (commOverPbc_ && !shiftForces.empty());

    RVec* shiftForce = (haveShiftForces ? &shiftForces[pbcForceShiftIndex_] : nullptr);

    dispatchTemplatedFunction(
            [&](auto usesScrewPbc, auto haveShiftForces) {
                accumulateReceivedForcesTemplated<usesScrewPbc, haveShiftForces>(*this, forces, shiftForce);
            },
            usesScrewPbc_,
            haveShiftForces);
}

void HaloExchange::initiateReceiveF()
{
    HaloMpiRequests& mpiRequests = mpiForceRequests_;

    /* Post all the non-blocking receives */
    for (DomainPairComm& dpc : domainPairComm_)
    {
        DomainCommBackward& receive = dpc.backward();

        if (receive.numAtoms() > 0)
        {
            /* We can reuse the send x buffer as the receive buffer for f,
             * since the received x need to be processed before f can be
             * calculated and communicated.
             */
            ddIreceiveDomain(receive, receive.rvecBuffer(), HaloMpiTag::F, &mpiRequests.receive);
        }
    }
}

void HaloExchange::initiateSendF(ArrayRef<const RVec> f)
{
    HaloMpiRequests& mpiRequests = mpiForceRequests_;

    GMX_ASSERT(mpiRequests.send.empty(),
               "All MPI Requests should have been handled before initiating sendF");

    /* Non-blocking send using direct force buffer pointers */
    for (DomainPairComm& dpc : domainPairComm_)
    {
        DomainCommForward& send = dpc.forward();

        if (send.numAtoms() > 0)
        {
            ddIsendDomain(send,
                          f.subArray(*send.atomRange().begin(), send.atomRange().size()),
                          HaloMpiTag::F,
                          &mpiRequests.send);
        }
    }
}

void HaloExchange::completeReceiveF(ArrayRef<RVec> forces, ArrayRef<RVec> shiftForces)
{
    HaloMpiRequests& mpiRequests = mpiForceRequests_;

    if (debug)
    {
        fprintf(debug, "Waiting for force receive: %d messages\n", int(mpiRequests.receive.size()));
    }

    /* Wait for all non-blocking communication to complete */
    mpiWaitall(mpiRequests.receive, mpiStatus_);

    mpiRequests.receive.clear();

    /* Reduce the received non-local forces with our local forces */
    for (DomainPairComm& dpc : domainPairComm_)
    {
        dpc.backward().accumulateReceivedForces(forces, shiftForces);
    }
}

// Complete the non-blocking sends of halo forces
void HaloExchange::completeSendF()
{
    HaloMpiRequests& mpiRequests = mpiForceRequests_;

    if (debug)
    {
        fprintf(debug, "Waiting for force send: %d messages\n", int(mpiRequests.send.size()));
    }

    /* Wait for all non-blocking sends to complete */
    mpiWaitall(mpiRequests.send, mpiStatus_);

    mpiRequests.send.clear();
}

void HaloExchange::moveF(ArrayRef<RVec> f, ArrayRef<RVec> shiftForces)
{
    initiateReceiveF();
    initiateSendF(f);
    completeReceiveF(f, shiftForces);
    completeSendF();
}

} // namespace gmx
