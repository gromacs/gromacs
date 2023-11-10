/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements kernels and dispatch for position restraints
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_LISTEDFORCES_POSITIONRESTRAINTS_HPP
#define NBLIB_LISTEDFORCES_POSITIONRESTRAINTS_HPP

#include "gromacs/pbcutil/pbc.h"

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/traits.h"
#include "nblib/pbc.hpp"

namespace nblib
{

template<class BasicVector>
inline BasicVector posresCOM(const Box& box, RefCoordScaling refcoord_scaling, PbcType pbcType, BasicVector com)
{
    BasicVector com_sc{ 0, 0, 0 };
    auto        boxMatrix = box.legacyMatrix();

    int npbcdim = numPbcDimensions(pbcType);

    if (refcoord_scaling == RefCoordScaling::Com)
    {
        for (int m = 0; m < npbcdim; m++)
        {
            for (int d = m; d < npbcdim; d++)
            {
                com_sc[m] += com[d] * boxMatrix[d][m];
            }
        }
    }

    return com_sc;
}

/*! \brief returns dx, rdist, and dpdl for functions posres() and fbposres()
 */
template<class BasicVector, class Pbc>
util::tuple<BasicVector, BasicVector, BasicVector> posres_dx(BasicVector     x,
                                                             BasicVector     pos0A,
                                                             BasicVector     pos0B,
                                                             BasicVector     comA_sc,
                                                             BasicVector     comB_sc,
                                                             real            lambda,
                                                             const Pbc&      pbc,
                                                             const Box&      box,
                                                             RefCoordScaling refcoord_scaling,
                                                             int             npbcdim)
{
    using ValueType = VectorValueType_t<BasicVector>;

    ValueType   posA      = 0;
    ValueType   posB      = 0;
    ValueType   ref       = 0.;
    BasicVector pos       = { 0, 0, 0 };
    BasicVector dpdl      = { 0, 0, 0 };
    BasicVector rdist     = { 0, 0, 0 };
    auto        legacyBox = box.legacyMatrix();

    ValueType L1 = 1.0 - lambda;

    for (int m = 0; m < 3; m++)
    {
        posA = pos0A[m];
        posB = pos0B[m];

        if (m < npbcdim)
        {
            switch (refcoord_scaling)
            {
                case RefCoordScaling::No:
                    ref      = 0;
                    rdist[m] = L1 * posA + lambda * posB;
                    dpdl[m]  = posB - posA;
                    break;

                case RefCoordScaling::All:
                    /* Box relative coordinates are stored for dimensions with pbc */
                    posA *= legacyBox[m][m];
                    posB *= legacyBox[m][m];
                    assert(npbcdim <= DIM);

                    for (int d = m + 1; d < npbcdim && d < 3; d++)
                    {
                        posA += pos0A[d] * legacyBox[d][m];
                        posB += pos0B[d] * legacyBox[d][m];
                    }

                    ref      = L1 * posA + lambda * posB;
                    rdist[m] = 0;
                    dpdl[m]  = posB - posA;
                    break;

                case RefCoordScaling::Com:
                    ref      = L1 * comA_sc[m] + lambda * comB_sc[m];
                    rdist[m] = L1 * posA + lambda * posB;
                    dpdl[m]  = comB_sc[m] - comA_sc[m] + posB - posA;
                    break;

                default: throw InputException("No such scaling method implemented");
            }
        }
        else
        {
            ref      = L1 * posA + lambda * posB;
            rdist[m] = 0;
            dpdl[m]  = posB - posA;
        }

        /* We do pbc_dx with ref+rdist,
         * since with only ref we can be up to half a box vector wrong.
         */
        pos[m] = ref + rdist[m];
    }

    BasicVector dx = { 0, 0, 0 };
    pbc.dxAiuc(x, pos, dx);
    return util::make_tuple(dx, rdist, dpdl);
}

/*! \brief returns dx, rdist, for functions posres() and fbposres()
 */
template<class BasicVector, class Pbc>
util::tuple<BasicVector, BasicVector> posres_dx(BasicVector     x,
                                                BasicVector     pos0,
                                                BasicVector     com_sc,
                                                const Pbc&      pbc,
                                                const Box&      box,
                                                RefCoordScaling refcoord_scaling,
                                                int             npbcdim)
{
    using ValueType = VectorValueType_t<BasicVector>;
    ValueType   posA;
    ValueType   ref = 0.;
    BasicVector pos;
    BasicVector rdist;
    auto        legacyBox = box.legacyMatrix();

    for (int m = 0; m < 3; m++)
    {
        posA = pos0[m];
        if (m < npbcdim)
        {
            switch (refcoord_scaling)
            {
                case RefCoordScaling::No:
                    ref      = 0;
                    rdist[m] = posA;
                    break;

                case RefCoordScaling::All:
                    /* Box relative coordinates are stored for dimensions with pbc */
                    posA *= legacyBox[m][m];
                    assert(npbcdim <= DIM);

                    for (int d = m + 1; d < npbcdim && d < 3; d++)
                    {
                        posA += pos0[d] * legacyBox[d][m];
                    }

                    ref      = posA;
                    rdist[m] = 0;
                    break;

                case RefCoordScaling::Com:
                    ref      = com_sc[m];
                    rdist[m] = posA;
                    break;

                default: throw InputException("No such scaling method implemented");
            }
        }
        else
        {
            ref      = posA;
            rdist[m] = 0;
        }

        /* We do pbc_dx with ref+rdist,
         * since with only ref we can be up to half a box vector wrong.
         */
        pos[m] = ref + rdist[m];
    }

    BasicVector dx;
    pbc.dxAiuc(x, pos, dx);

    return util::make_tuple(dx, rdist);
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_POSITIONRESTRAINTS_HPP
