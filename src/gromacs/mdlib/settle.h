/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares interface to SETTLE code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_SETTLE_H
#define GMX_MDLIB_SETTLE_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct InteractionList;
struct t_inputrec;
struct t_pbc;

namespace gmx
{

template<typename>
class ArrayRef;
template<typename>
class ArrayRefWithPadding;
enum class ConstraintVariable : int;

/* \brief Parameters for SETTLE algorithm.
 *
 * \todo Remove duplicates, check if recomputing makes more sense in some cases.
 * \todo Move the projection parameters into separate structure.
 */
struct SettleParameters
{
    //! Mass of oxygen atom
    real mO;
    //! Mass of hydrogen atom
    real mH;
    //! Relative hydrogen mass (i.e. mH/(mO+2*mH))
    real wh;
    //! Target distance between oxygen and hydrogen
    real dOH;
    //! Target distance between hydrogen atoms
    real dHH;
    //! Double relative H mass times height of the H-O-H triangle.
    real ra;
    //! Height of the H-O-H triangle minus ra.
    real rb;
    //! Half the H-H distance.
    real rc;
    //! Reciprocal H-H distance
    real irc2;

    /* Parameters below are used for projection */
    //! Reciprocal oxygen mass
    real imO;
    //! Reciprocal hydrogen mass
    real imH;
    //! Reciprocal O-H distance
    real invdOH;
    //! Reciprocal H-H distance (again)
    real invdHH;
    //! Reciprocal projection matrix
    matrix invmat;
};

/*! \brief Computes and returns settle parameters.
 *
 * \param[in]  mO     Mass of oxygen atom
 * \param[in]  mH     Mass of hydrogen atom
 * \param[in]  invmO  Reciprocal mass of oxygen atom
 * \param[in]  invmH  Reciprocal mass of hydrogen atom
 * \param[in]  dOH    Target O-H bond length
 * \param[in]  dHH    Target H-H bond length
 */
SettleParameters settleParameters(real mO, real mH, real invmO, real invmH, real dOH, real dHH);

/*! \libinternal
 * \brief Data for executing SETTLE constraining
 */
class SettleData
{
public:
    //! Constructor
    SettleData(const gmx_mtop_t& mtop);

    //! Sets the constraints from the interaction list and the masses
    void setConstraints(const InteractionList&    il_settle,
                        int                       numHomeAtoms,
                        gmx::ArrayRef<const real> masses,
                        gmx::ArrayRef<const real> inverseMasses);

    //! Returns settle parameters for constraining coordinates and forces
    const SettleParameters& parametersMassWeighted() const { return parametersMassWeighted_; }

    //! Returns settle parameters for constraining forces: all masses are set to 1
    const SettleParameters& parametersAllMasses1() const { return parametersAllMasses1_; }

    //! Returns the number of SETTLEs
    int numSettles() const { return numSettles_; }

    //! Returns a pointer to the indices of oxygens for each SETTLE
    const int* ow1() const { return ow1_.data(); }
    //! Returns a pointer to the indices of the first hydrogen for each SETTLE
    const int* hw2() const { return hw2_.data(); }
    //! Returns a pointer to the indices of the second hydrogen for each SETTLE
    const int* hw3() const { return hw3_.data(); }
    //! Returns a pointer to the virial contribution factor (either 1 or 0) for each SETTLE
    const real* virfac() const { return virfac_.data(); }

    //! Returns whether we should use SIMD intrinsics code
    bool useSimd() const { return useSimd_; }

private:
    //! Parameters for SETTLE for coordinates
    SettleParameters parametersMassWeighted_;
    //! Parameters with all masses 1, for forces
    SettleParameters parametersAllMasses1_;

    //! The number of settles on our rank
    int numSettles_;

    //! Index to OW1 atoms, size numSettles_ + SIMD padding
    std::vector<int, gmx::AlignedAllocator<int>> ow1_;
    //! Index to HW2 atoms, size numSettles_ + SIMD padding
    std::vector<int, gmx::AlignedAllocator<int>> hw2_;
    //! Index to HW3 atoms, size numSettles_ + SIMD padding
    std::vector<int, gmx::AlignedAllocator<int>> hw3_;
    //! Virial factor 0 or 1, size numSettles_ + SIMD padding
    std::vector<real, gmx::AlignedAllocator<real>> virfac_;

    //! Tells whether we will use SIMD intrinsics code
    bool useSimd_;
};

/*! \brief Constrain coordinates using SETTLE.
 * Can be called on any number of threads.
 */
void csettle(const SettleData&               settled,     /* The SETTLE structure */
             int                             nthread,     /* The number of threads used */
             int                             thread,      /* Our thread index */
             const t_pbc*                    pbc,         /* PBC data pointer, can be NULL */
             ArrayRefWithPadding<const RVec> x,           /* Reference coordinates */
             ArrayRefWithPadding<RVec>       xprime,      /* New coords, to be settled */
             real                            invdt,       /* 1/delta_t */
             ArrayRefWithPadding<RVec>       v,           /* Also constrain v if v!=NULL */
             bool                            bCalcVirial, /* Calculate the virial contribution */
             tensor                          vir_r_m_dr,  /* sum r x m delta_r */
             bool*                           bErrorHasOccurred /* True if a settle error occurred */
);

/*! \brief Analytical algorithm to subtract the components of derivatives
 * of coordinates working on settle type constraint.
 */
void settle_proj(const SettleData&    settled,
                 ConstraintVariable   econq,
                 int                  nsettle,
                 const int            iatoms[],
                 const t_pbc*         pbc, /* PBC data pointer, can be NULL  */
                 ArrayRef<const RVec> x,
                 ArrayRef<RVec>       der,
                 ArrayRef<RVec>       derp,
                 int                  CalcVirAtomEnd,
                 tensor               vir_r_m_dder);

} // namespace gmx

#endif
