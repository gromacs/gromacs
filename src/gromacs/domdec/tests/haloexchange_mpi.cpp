/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief Tests for the halo exchange
 *
 *  The test sets up a 2D rank topology and performs a coordinate halo
 *  exchange (using the pre-existing CPU codepath), with 2 pulses in
 *  the first dimension and 1 pulse in the second. Each pulse involves
 *  a few non-contiguous indices. The sending rank, atom number and
 *  spatial 3D index are encoded in the x values, to allow correctness
 *  checking following the halo exchange.
 *
 * \todo Add more test variations
 * \todo Port to GPU codepath
 *
 * \author Alan Gray <alang@nvidia.com>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include <array>

#include <gtest/gtest.h>

#include "gromacs/domdec/atomdistribution.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/mdtypes/inputrec.h"

#include "testutils/mpitest.h"

namespace gmx
{
namespace
{

/*! \brief Get encoded numerical value for sending rank, atom number and spatial 3D index
 *
 * \param [in] sendRank         MPI rank of sender
 * \param [in] atomNumber       Atom number
 * \param [in] spatial3dIndex   Spatial 3D Index
 *
 * \returns                     Encoded value
 */
float encodedValue(const int sendRank, const int atomNumber, const int spatial3dIndex)
{
    return sendRank * 1000 + atomNumber * 100 + spatial3dIndex;
}

/*! \brief Initialize halo array
 *
 * \param [in] x              Atom coordinate data array
 * \param [in] numHomeAtoms   Number of home atoms
 * \param [in] numAtomsTotal  Total number of atoms, including halo
 */
void initHaloData(RVec* x, const int numHomeAtoms, const int numAtomsTotal)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < numAtomsTotal; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            x[i][j] = i < numHomeAtoms ? encodedValue(rank, i, j) : -1;
        }
    }
}

/*! \brief Define 2D rank topology with 4 MPI tasks
 *
 *    -----
 *   | 2 3 |
 *   | 0 1 |
 *    -----
 *
 * \param [in] dd  Domain decomposition object
 */
void define2dRankTopology(gmx_domdec_t* dd)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch (rank)
    {
        case 0:
            dd->neighbor[0][0] = 1;
            dd->neighbor[0][1] = 1;
            dd->neighbor[1][0] = 2;
            dd->neighbor[1][1] = 2;
            break;
        case 1:
            dd->neighbor[0][0] = 0;
            dd->neighbor[0][1] = 0;
            dd->neighbor[1][0] = 3;
            dd->neighbor[1][1] = 3;
            break;
        case 2:
            dd->neighbor[0][0] = 3;
            dd->neighbor[0][1] = 3;
            dd->neighbor[1][0] = 0;
            dd->neighbor[1][1] = 0;
            break;
        case 3:
            dd->neighbor[0][0] = 2;
            dd->neighbor[0][1] = 2;
            dd->neighbor[1][0] = 1;
            dd->neighbor[1][1] = 1;
            break;
    }
}

/*! \brief Define a 2D halo with 2 pulses in the first dimension
 *
 * \param [in] dd      Domain decomposition object
 * \param [in] indvec  Vector of index vectors
 */
void define2dHaloWith2PulsesInDim1(gmx_domdec_t* dd, std::vector<gmx_domdec_ind_t> indvec)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> indexvec;
    gmx_domdec_ind_t ind;

    dd->ndim  = 2;
    int nzone = 1;
    for (int dimIndex = 0; dimIndex < dd->ndim; dimIndex++)
    {

        // Set up indices involved in halo
        indexvec.clear();
        indvec.clear();

        dd->comm->cd[dimIndex].receiveInPlace = true;
        dd->dim[dimIndex]                     = 0;
        dd->ci[dimIndex]                      = rank;

        // First pulse involves (arbitrary) indices 1 and 3
        indexvec.push_back(1);
        indexvec.push_back(3);

        ind.index            = indexvec;
        ind.nsend[nzone + 1] = 2;
        ind.nrecv[nzone + 1] = 2;
        indvec.push_back(ind);

        if (dimIndex == 0) // Add another pulse with (arbitrary) indices 4,5,7
        {
            indexvec.clear();

            dd->comm->cd[dimIndex].ind = indvec;

            indexvec.push_back(4);
            indexvec.push_back(5);
            indexvec.push_back(7);

            ind.index            = indexvec;
            ind.nsend[nzone + 1] = 3;
            ind.nrecv[nzone + 1] = 3;
            indvec.push_back(ind);
        }

        dd->comm->cd[dimIndex].ind = indvec;

        nzone += nzone;
    }
}

/*! \brief Check results for above-defined 2D halo with 2 pulses in the first dimension
 *
 * \param [in] x             Atom coordinate data array
 * \param [in] dd            Domain decomposition object
 * \param [in] numHomeAtoms  Number of home atoms
 */
void checkResults2dHaloWith2PulsesInDim1(const RVec* x, const gmx_domdec_t* dd, const int numHomeAtoms)
{
    // Check results are expected from values encoded in x data
    for (int j = 0; j < DIM; j++)
    {
        // First Pulse in first dim: atoms 1 and 3 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms][j], encodedValue(dd->neighbor[0][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 1][j], encodedValue(dd->neighbor[0][0], 3, j));
        // Second Pulse in first dim: atoms 4,5,7 from forward horizontal neighbour
        EXPECT_EQ(x[numHomeAtoms + 2][j], encodedValue(dd->neighbor[0][0], 4, j));
        EXPECT_EQ(x[numHomeAtoms + 3][j], encodedValue(dd->neighbor[0][0], 5, j));
        EXPECT_EQ(x[numHomeAtoms + 4][j], encodedValue(dd->neighbor[0][0], 7, j));
        // First Pulse in second dim: atoms 1 and 3 from forward vertical neighbour
        EXPECT_EQ(x[numHomeAtoms + 5][j], encodedValue(dd->neighbor[1][0], 1, j));
        EXPECT_EQ(x[numHomeAtoms + 6][j], encodedValue(dd->neighbor[1][0], 3, j));
    }
}


TEST(HaloExchangeTest, Coordinates2dHaloWith2PulsesInDim1)
{
    GMX_MPI_TEST(4);

    // Set up atom data
    const int numHomeAtoms  = 10;
    const int numHaloAtoms  = 7;
    const int numAtomsTotal = numHomeAtoms + numHaloAtoms;
    RVec      x[numAtomsTotal];
    initHaloData(x, numHomeAtoms, numAtomsTotal);

    // Set up dd
    t_inputrec   ir;
    gmx_domdec_t dd(ir);
    dd.mpi_comm_all = MPI_COMM_WORLD;
    gmx_domdec_comm_t comm;
    dd.comm                      = &comm;
    dd.unitCellInfo.haveScrewPBC = false;

    DDAtomRanges atomRanges;
    atomRanges.setEnd(DDAtomRanges::Type::Home, numHomeAtoms);
    dd.comm->atomRanges = atomRanges;

    define2dRankTopology(&dd);

    std::vector<gmx_domdec_ind_t> indvec;
    define2dHaloWith2PulsesInDim1(&dd, indvec);

    // Perform halo exchange
    matrix box = { { 0., 0., 0. } };
    dd_move_x(&dd, box, static_cast<ArrayRef<RVec>>(x), nullptr);

    // Check results
    checkResults2dHaloWith2PulsesInDim1(x, &dd, numHomeAtoms);
}

} // namespace
} // namespace gmx
