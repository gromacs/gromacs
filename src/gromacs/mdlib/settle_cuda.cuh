/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declares class for CUDA implementation of SETTLE
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_SETTLE_CUDA_CUH
#define GMX_MDLIB_SETTLE_CUDA_CUH

#include "gmxpre.h"

#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/* \brief Parameters for SETTLE algorithm.
 *
 * Following structure and subroutine are copy-pasted from the CPU version of SETTLE.
 * This is a temporary solution and they will be recycled in more appropriate way.
 * \todo Remove duplicates, check if recomputing makes more sense in some cases.
 * \todo Move the projection parameters into separate structure.
 */
struct SettleParameters
{
    //! Mass of oxygen atom
    float mO;
    //! Mass of hydrogen atom
    float mH;
    //! Relative hydrogen mass (i.e. mH/(mO+2*mH))
    float wh;
    //! Target distance between oxygen and hydrogen
    float dOH;
    //! Target distance between hydrogen atoms
    float dHH;
    //! Double relative H mass times height of the H-O-H triangle.
    float ra;
    //! Height of the H-O-H triangle minus ra.
    float rb;
    //! Half the H-H distance.
    float rc;
    //! Reciprocal H-H distance
    float irc2;

    /* Parameters below are used for projection */
    //! Reciprocal oxygen mass
    float imO;
    //! Reciprocal hydrogen mass
    float imH;
    //! Reciprocal O-H distance
    float invdOH;
    //! Reciprocal H-H distance (again)
    float invdHH;
    //! Reciprocal projection matrix
    matrix invmat;
};

/*! \brief Initializes a projection matrix.
 *
 * \todo This is identical to to CPU code. Unification is needed.
 *
 * \param[in]  invmO                  Reciprocal oxygen mass
 * \param[in]  invmH                  Reciprocal hydrogen mass
 * \param[in]  dOH                    Target O-H bond length
 * \param[in]  dHH                    Target H-H bond length
 * \param[out] inverseCouplingMatrix  Inverse bond coupling matrix for the projection version of SETTLE
 */
static void initializeProjectionMatrix(const real invmO,
                                       const real invmH,
                                       const real dOH,
                                       const real dHH,
                                       matrix     inverseCouplingMatrix)
{
    // We normalize the inverse masses with invmO for the matrix inversion.
    // so we can keep using masses of almost zero for frozen particles,
    // without running out of the float range in invertMatrix.
    double invmORelative = 1.0;
    double invmHRelative = invmH / static_cast<double>(invmO);
    double distanceRatio = dHH / static_cast<double>(dOH);

    /* Construct the constraint coupling matrix */
    matrix mat;
    mat[0][0] = invmORelative + invmHRelative;
    mat[0][1] = invmORelative * (1.0 - 0.5 * gmx::square(distanceRatio));
    mat[0][2] = invmHRelative * 0.5 * distanceRatio;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = invmHRelative + invmHRelative;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    gmx::invertMatrix(mat, inverseCouplingMatrix);

    msmul(inverseCouplingMatrix, 1 / invmO, inverseCouplingMatrix);
}

/*! \brief Initializes settle parameters.
 *
 * \todo This is identical to the CPU code. Unification is needed.
 *
 * \param[out] p    SettleParameters structure to initialize
 * \param[in]  mO   Mass of oxygen atom
 * \param[in]  mH   Mass of hydrogen atom
 * \param[in]  dOH  Target O-H bond length
 * \param[in]  dHH  Target H-H bond length
 */
gmx_unused // Temporary solution to keep clang happy
        static void
        initSettleParameters(SettleParameters* p, const real mO, const real mH, const real dOH, const real dHH)
{
    // We calculate parameters in double precision to minimize errors.
    // The velocity correction applied during SETTLE coordinate constraining
    // introduces a systematic error of approximately 1 bit per atom,
    // depending on what the compiler does with the code.
    double wohh;

    real invmO = 1.0 / mO;
    real invmH = 1.0 / mH;

    p->mO     = mO;
    p->mH     = mH;
    wohh      = mO + 2.0 * mH;
    p->wh     = mH / wohh;
    p->dOH    = dOH;
    p->dHH    = dHH;
    double rc = dHH / 2.0;
    double ra = 2.0 * mH * std::sqrt(dOH * dOH - rc * rc) / wohh;
    p->rb     = std::sqrt(dOH * dOH - rc * rc) - ra;
    p->rc     = rc;
    p->ra     = ra;
    p->irc2   = 1.0 / dHH;

    // For projection: inverse masses and coupling matrix inversion
    p->imO = invmO;
    p->imH = invmH;

    p->invdOH = 1.0 / dOH;
    p->invdHH = 1.0 / dHH;

    initializeProjectionMatrix(invmO, invmH, dOH, dHH, p->invmat);
}

/*! \internal \brief Class with interfaces and data for CUDA version of SETTLE. */
class SettleCuda
{

public:
    /*! \brief Create SETTLE object
     *
     *  Extracts masses for oxygen and hydrogen as well as the O-H and H-H target distances
     *  from the topology data (mtop), check their values for consistency and calls the
     *  following constructor.
     *
     * \param[in] mtop           Topology of the system to gen the masses for O and H atoms and
     *                           target O-H and H-H distances. These values are also checked for
     *                           consistency.
     * \param[in] commandStream  Device stream to use.
     */
    SettleCuda(const gmx_mtop_t& mtop, CommandStream commandStream);

    ~SettleCuda();

    /*! \brief Apply SETTLE.
     *
     * Applies SETTLE to coordinates and velocities, stored on GPU. Data at pointers d_xp and
     * d_v change in the GPU memory. The results are not automatically copied back to the CPU
     * memory. Method uses this class data structures which should be updated when needed using
     * update method.
     *
     * \param[in]     d_x               Coordinates before timestep (in GPU memory)
     * \param[in,out] d_xp              Coordinates after timestep (in GPU memory). The
     *                                  resulting constrained coordinates will be saved here.
     * \param[in]     updateVelocities  If the velocities should be updated.
     * \param[in,out] d_v               Velocities to update (in GPU memory, can be nullptr
     *                                  if not updated)
     * \param[in]     invdt             Reciprocal timestep (to scale Lagrange
     *                                  multipliers when velocities are updated)
     * \param[in]     computeVirial     If virial should be updated.
     * \param[in,out] virialScaled      Scaled virial tensor to be updated.
     */
    void apply(const float3* d_x,
               float3*       d_xp,
               const bool    updateVelocities,
               float3*       d_v,
               const real    invdt,
               const bool    computeVirial,
               tensor        virialScaled);

    /*! \brief
     * Update data-structures (e.g. after NB search step).
     *
     * Updates the constraints data and copies it to the GPU. Should be
     * called if the particles were sorted, redistributed between domains, etc.
     * Does not recycle the data preparation routines from the CPU version.
     * All three atoms from single water molecule should be handled by the same GPU.
     *
     * SETTLEs atom ID's is taken from idef.il[F_SETTLE].iatoms.
     *
     * \param[in] idef    System topology
     * \param[in] md      Atoms data. Can be used to update masses if needed (not used now).
     */
    void set(const t_idef& idef, const t_mdatoms& md);

    /*! \brief
     * Update PBC data.
     *
     * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \todo PBC should not be handled by constraints.
     *
     * \param[in] pbc The PBC data in t_pbc format.
     */
    void setPbc(const t_pbc* pbc);


private:
    //! CUDA stream
    CommandStream commandStream_;
    //! Periodic boundary data
    PbcAiuc pbcAiuc_;

    //! Scaled virial tensor (9 reals, GPU)
    std::vector<float> h_virialScaled_;
    //! Scaled virial tensor (9 reals, GPU)
    float* d_virialScaled_;

    //! Number of settles
    int numSettles_ = 0;

    //! Indexes of atoms (.x for oxygen, .y and.z for hydrogens, CPU)
    std::vector<int3> h_atomIds_;
    //! Indexes of atoms (.x for oxygen, .y and.z for hydrogens, GPU)
    int3* d_atomIds_;
    //! Current size of the array of atom IDs
    int numAtomIds_ = -1;
    //! Allocated size for the array of atom IDs
    int numAtomIdsAlloc_ = -1;

    //! Settle parameters
    SettleParameters settleParameters_;
};

} // namespace gmx

#endif
