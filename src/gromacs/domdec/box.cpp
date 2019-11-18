/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief This file defines functions used by the domdec module
 * for (bounding) box and pbc information generation.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "box.h"

#include "config.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nsgrid.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"

#include "domdec_internal.h"

/*! \brief Calculates the average and standard deviation in 3D of atoms */
static void calc_pos_av_stddev(gmx::ArrayRef<const gmx::RVec> x, rvec av, rvec stddev, const MPI_Comm* mpiCommunicator)
{
    dvec s1, s2;

    clear_dvec(s1);
    clear_dvec(s2);

    for (const gmx::RVec& coord : x)
    {
        for (int d = 0; d < DIM; d++)
        {
            s1[d] += coord[d];
            s2[d] += coord[d] * coord[d];
        }
    }

    /* With mpiCommunicator != nullptr, x.size() is the home atom count */
    int numAtoms = x.size();
#if GMX_MPI
    if (mpiCommunicator)
    {
        constexpr int c_bufSize = 7;
        double        sendBuffer[c_bufSize];
        double        receiveBuffer[c_bufSize];

        for (int d = 0; d < DIM; d++)
        {
            sendBuffer[d]       = s1[d];
            sendBuffer[DIM + d] = s2[d];
        }
        sendBuffer[6] = numAtoms;

        MPI_Allreduce(sendBuffer, receiveBuffer, c_bufSize, MPI_DOUBLE, MPI_SUM, *mpiCommunicator);

        for (int d = 0; d < DIM; d++)
        {
            s1[d] = receiveBuffer[d];
            s2[d] = receiveBuffer[DIM + d];
        }
        numAtoms = gmx::roundToInt(receiveBuffer[6]);
    }
#else  // GMX_MPI
    GMX_UNUSED_VALUE(mpiCommunicator);
#endif // GMX_MPI

    dsvmul(1.0 / numAtoms, s1, s1);
    dsvmul(1.0 / numAtoms, s2, s2);

    for (int d = 0; d < DIM; d++)
    {
        av[d]     = s1[d];
        stddev[d] = std::sqrt(s2[d] - s1[d] * s1[d]);
    }
}

/*! \brief Determines if dimensions require triclinic treatment and stores this info in ddbox */
static void set_tric_dir(const ivec* dd_nc, gmx_ddbox_t* ddbox, const matrix box)
{
    int   npbcdim, d, i, j;
    rvec *v, *normal;
    real  dep, inv_skew_fac2;

    npbcdim = ddbox->npbcdim;
    normal  = ddbox->normal;
    for (d = 0; d < DIM; d++)
    {
        ddbox->tric_dir[d] = 0;
        for (j = d + 1; j < npbcdim; j++)
        {
            if (box[j][d] != 0)
            {
                ddbox->tric_dir[d] = 1;
                if (dd_nc != nullptr && (*dd_nc)[j] > 1 && (*dd_nc)[d] == 1)
                {
                    gmx_fatal(FARGS,
                              "Domain decomposition has not been implemented for box vectors that "
                              "have non-zero components in directions that do not use domain "
                              "decomposition: ncells = %d %d %d, box vector[%d] = %f %f %f",
                              (*dd_nc)[XX], (*dd_nc)[YY], (*dd_nc)[ZZ], j + 1, box[j][XX],
                              box[j][YY], box[j][ZZ]);
                }
            }
        }

        /* Construct vectors v for dimension d that are perpendicular
         * to the triclinic plane of dimension d. Each vector v[i] has
         * v[i][i]=1 and v[i][d]!=0 for triclinic dimensions, while the third
         * component is zero. These are used for computing the distance
         * to a triclinic plane given the distance along dimension d.
         * Set the trilinic skewing factor that translates
         * the thickness of a slab perpendicular to this dimension
         * into the real thickness of the slab.
         */
        if (ddbox->tric_dir[d])
        {
            inv_skew_fac2 = 1;
            v             = ddbox->v[d];
            if (d == XX || d == YY)
            {
                /* Normalize such that the "diagonal" is 1 */
                svmul(1 / box[d + 1][d + 1], box[d + 1], v[d + 1]);
                for (i = 0; i < d; i++)
                {
                    v[d + 1][i] = 0;
                }
                inv_skew_fac2 += gmx::square(v[d + 1][d]);
                if (d == XX)
                {
                    /* Normalize such that the "diagonal" is 1 */
                    svmul(1 / box[d + 2][d + 2], box[d + 2], v[d + 2]);
                    /* Set v[d+2][d+1] to zero by shifting along v[d+1] */
                    dep = v[d + 2][d + 1] / v[d + 1][d + 1];
                    for (i = 0; i < DIM; i++)
                    {
                        v[d + 2][i] -= dep * v[d + 1][i];
                    }
                    inv_skew_fac2 += gmx::square(v[d + 2][d]);

                    cprod(v[d + 1], v[d + 2], normal[d]);
                }
                else
                {
                    /* cross product with (1,0,0) */
                    normal[d][XX] = 0;
                    normal[d][YY] = v[d + 1][ZZ];
                    normal[d][ZZ] = -v[d + 1][YY];
                }
                if (debug)
                {
                    fprintf(debug, "box[%d]  %.3f %.3f %.3f\n", d, box[d][XX], box[d][YY], box[d][ZZ]);
                    for (i = d + 1; i < DIM; i++)
                    {
                        fprintf(debug, "  v[%d]  %.3f %.3f %.3f\n", i, v[i][XX], v[i][YY], v[i][ZZ]);
                    }
                }
            }
            ddbox->skew_fac[d] = 1.0 / std::sqrt(inv_skew_fac2);
            /* Set the normal vector length to skew_fac */
            dep = ddbox->skew_fac[d] / norm(normal[d]);
            svmul(dep, normal[d], normal[d]);

            if (debug)
            {
                fprintf(debug, "skew_fac[%d] = %f\n", d, ddbox->skew_fac[d]);
                fprintf(debug, "normal[%d]  %.3f %.3f %.3f\n", d, normal[d][XX], normal[d][YY],
                        normal[d][ZZ]);
            }
        }
        else
        {
            ddbox->skew_fac[d] = 1;

            for (i = 0; i < DIM; i++)
            {
                clear_rvec(ddbox->v[d][i]);
                ddbox->v[d][i][i] = 1;
            }
            clear_rvec(normal[d]);
            normal[d][d] = 1;
        }
    }
}

/*! \brief This function calculates bounding box and pbc info and populates ddbox */
static void low_set_ddbox(int                            numPbcDimensions,
                          int                            numBoundedDimensions,
                          const ivec*                    dd_nc,
                          const matrix                   box,
                          bool                           calculateUnboundedSize,
                          gmx::ArrayRef<const gmx::RVec> x,
                          const MPI_Comm*                mpiCommunicator,
                          gmx_ddbox_t*                   ddbox)
{
    rvec av, stddev;
    real b0, b1;
    int  d;

    ddbox->npbcdim     = numPbcDimensions;
    ddbox->nboundeddim = numBoundedDimensions;

    for (d = 0; d < numBoundedDimensions; d++)
    {
        ddbox->box0[d]     = 0;
        ddbox->box_size[d] = box[d][d];
    }

    if (ddbox->nboundeddim < DIM && calculateUnboundedSize)
    {
        calc_pos_av_stddev(x, av, stddev, mpiCommunicator);

        /* GRID_STDDEV_FAC * stddev
         * gives a uniform load for a rectangular block of cg's.
         * For a sphere it is not a bad approximation for 4x1x1 up to 4x2x2.
         */
        for (d = ddbox->nboundeddim; d < DIM; d++)
        {
            b0 = av[d] - GRID_STDDEV_FAC * stddev[d];
            b1 = av[d] + GRID_STDDEV_FAC * stddev[d];
            if (debug)
            {
                fprintf(debug, "Setting global DD grid boundaries to %f - %f\n", b0, b1);
            }
            ddbox->box0[d]     = b0;
            ddbox->box_size[d] = b1 - b0;
        }
    }

    set_tric_dir(dd_nc, ddbox, box);
}

void set_ddbox(const gmx_domdec_t&            dd,
               bool                           masterRankHasTheSystemState,
               const matrix                   box,
               bool                           calculateUnboundedSize,
               gmx::ArrayRef<const gmx::RVec> x,
               gmx_ddbox_t*                   ddbox)
{
    if (!masterRankHasTheSystemState || DDMASTER(dd))
    {
        bool needToReduceCoordinateData     = (!masterRankHasTheSystemState && dd.nnodes > 1);
        gmx::ArrayRef<const gmx::RVec> xRef = constArrayRefFromArray(
                x.data(), masterRankHasTheSystemState ? x.size() : dd.comm->atomRanges.numHomeAtoms());

        low_set_ddbox(dd.unitCellInfo.npbcdim, dd.unitCellInfo.numBoundedDimensions, &dd.nc, box,
                      calculateUnboundedSize, xRef,
                      needToReduceCoordinateData ? &dd.mpi_comm_all : nullptr, ddbox);
    }

    if (masterRankHasTheSystemState)
    {
        dd_bcast(&dd, sizeof(gmx_ddbox_t), ddbox);
    }
}

void set_ddbox_cr(const t_commrec&               cr,
                  const ivec*                    dd_nc,
                  const t_inputrec&              ir,
                  const matrix                   box,
                  gmx::ArrayRef<const gmx::RVec> x,
                  gmx_ddbox_t*                   ddbox)
{
    if (MASTER(&cr))
    {
        low_set_ddbox(ePBC2npbcdim(ir.ePBC), inputrec2nboundeddim(&ir), dd_nc, box, true, x, nullptr, ddbox);
    }

    gmx_bcast(sizeof(gmx_ddbox_t), ddbox, &cr);
}
