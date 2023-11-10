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
/*! \internal \file
 * \brief
 * Implements a RAII encapsulation objects for PbcAiuc and t_pbc that
 * can be passed to listed forces.
 * Note: this is a header implementation to make inlining easier if needed.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#ifndef NBLIB_PBC_HPP
#define NBLIB_PBC_HPP

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc.h"

#include "nblib/box.h"
#include "nblib/listed_forces/traits.h"
#include "nblib/util/annotation.hpp"

namespace nblib
{

/*! \brief life-time manager for pbcAiuc
 *
 */
class PbcHolderAiuc
{
    HOST_DEVICE_FUN static inline int xyzToShiftIndex(int x, int y, int z)
    {
        return (gmx::detail::c_nBoxX * (gmx::detail::c_nBoxY * ((z) + gmx::c_dBoxZ) + (y) + gmx::c_dBoxY)
                + (x) + gmx::c_dBoxX);
    }

public:
    PbcHolderAiuc() = default;

    //! \brief ctor
    PbcHolderAiuc(PbcType pbcType, const Box& box)
    {
        setPbcAiuc(numPbcDimensions(pbcType), box.legacyMatrix(), &pbcAiuc);
    }

    //! \brief calculate pbc-aware r1-r2
    template<class Vector1, class Vector2>
    HOST_DEVICE_FUN HOST_DEVICE_INLINE int dxAiuc(const Vector1 r1, const Vector1 r2, Vector2& dr) const
    {
        using ValueType = VectorValueType_t<Vector2>;

        dr[0] = r1[0] - r2[0];
        dr[1] = r1[1] - r2[1];
        dr[2] = r1[2] - r2[2];

        ValueType shz = std::rint(dr[2] * pbcAiuc.invBoxDiagZ);
        dr[0] -= shz * pbcAiuc.boxZX;
        dr[1] -= shz * pbcAiuc.boxZY;
        dr[2] -= shz * pbcAiuc.boxZZ;

        ValueType shy = std::rint(dr[1] * pbcAiuc.invBoxDiagY);
        dr[0] -= shy * pbcAiuc.boxYX;
        dr[1] -= shy * pbcAiuc.boxYY;

        ValueType shx = std::rint(dr[0] * pbcAiuc.invBoxDiagX);
        dr[0] -= shx * pbcAiuc.boxXX;

        return xyzToShiftIndex(-std::lrint(shx), -std::lrint(shy), -std::lrint(shz));
    }

private:
    PbcAiuc pbcAiuc;
};

/*! \brief life-time manager for t_pbc
 *
 */
class PbcHolder
{
public:
    //! \brief ctor
    explicit PbcHolder(PbcType pbcType, const Box& box)
    {
        set_pbc(&m_pbc, pbcType, box.legacyMatrix());
    }

    //! \brief calculate pbc-aware r1-r2, including the shift index
    template<class Vector1, class Vector2>
    inline int dxAiuc(const Vector1& r1, const Vector1& r2, Vector2& dr) const
    {
        return pbc_dx_aiuc(&m_pbc, r1.as_vec(), r2.as_vec(), dr.as_vec());
    }

private:
    t_pbc m_pbc;
};

/*! \brief dummy class used to turn off Pbc in listed forces
 *
 */
class NoPbc
{
public:
    //! \brief calculate r1-r2, ignore pbc
    template<class Vector1, class Vector2>
    HOST_DEVICE_FUN inline int dxAiuc(const Vector1& r1, const Vector1& r2, Vector2& dr) const
    {
        dr = r1 - r2;
        return gmx::c_centralShiftIndex;
    }
};

} // namespace nblib
#endif // NBLIB_PBC_HPP
