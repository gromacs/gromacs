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
/*! \libinternal \file
 *
 * \brief Declares the class for CUDA implementation of LINCS.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LINCS_CUDA_CUH
#define GMX_MDLIB_LINCS_CUDA_CUH

#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc_aiuc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/* \brief LINCS parameters and GPU pointers
 *
 * This is used to accumulate all the parameters and pointers so they can be passed
 * to the GPU as a single structure.
 *
 */
struct LincsCudaKernelParameters
{
    //! Periodic boundary data
    PbcAiuc pbcAiuc;
    //! Order of expansion when inverting the matrix
    int expansionOrder;
    //! Number of iterations used to correct the projection
    int numIterations;
    //! 1/mass for all atoms (GPU)
    float* d_inverseMasses;
    //! Scaled virial tensor (6 floats: [XX, XY, XZ, YY, YZ, ZZ], GPU)
    float* d_virialScaled;
    /*! \brief Total number of threads.
     *
     *  This covers all constraints and gaps in the ends of the thread blocks
     *  that are necessary to avoid inter-block synchronizations.
     *  Should be a multiple of block size (the last block is filled with dummy to the end).
     */
    int numConstraintsThreads;
    //! List of constrained atoms (GPU memory)
    int2* d_constraints;
    //! Equilibrium distances for the constraints (GPU)
    float* d_constraintsTargetLengths;
    //! Number of constraints, coupled with the current one (GPU)
    int* d_coupledConstraintsCounts;
    //! List of coupled with the current one (GPU)
    int* d_coupledConstraintsIndices;
    //! Elements of the coupling matrix.
    float* d_matrixA;
    //! Mass factors (GPU)
    float* d_massFactors;
};

/*! \internal \brief Class with interfaces and data for CUDA version of LINCS. */
class LincsCuda
{

public:
    /*! \brief Constructor.
     *
     * \param[in] numIterations    Number of iteration for the correction of the projection.
     * \param[in] expansionOrder   Order of the matrix inversion algorithm.
     * \param[in] commandStream    Device command stream.
     */
    LincsCuda(int numIterations, int expansionOrder, CommandStream commandStream);
    /*! \brief Destructor.*/
    ~LincsCuda();

    /*! \brief Apply LINCS.
     *
     * Applies LINCS to coordinates and velocities, stored on GPU.
     * The results are not automatically copied back to the CPU memory.
     * Method uses this class data structures which should be updated
     * when needed using set() and setPbc() method.
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
     * This version uses common data formats so it can be called from anywhere
     * in the code. Does not recycle the data preparation routines from the CPU
     * version. Works only with simple case when all the constraints in idef are
     * are handled by a single GPU. Triangles are not handled as special case.
     *
     * Information about constraints is taken from:
     *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indexes (i1, i2)
     *     idef.iparams[T].constr.dA --- target length for constraint of type T
     * From t_mdatom, the code takes:
     *     md.invmass  --- array of inverse square root of masses for each atom in the system.
     *
     * \param[in] idef  Local topology data to get information on constraints from.
     * \param[in] md    Atoms data to get atom masses from.
     */
    void set(const t_idef& idef, const t_mdatoms& md);

    /*! \brief
     * Update PBC data.
     *
     * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \todo Remove this method. LINCS should not manage PBC.
     *
     * \param[in] pbc The PBC data in t_pbc format.
     */
    void setPbc(const t_pbc* pbc);


private:
    //! CUDA stream
    CommandStream commandStream_;

    //! Parameters and pointers, passed to the CUDA kernel
    LincsCudaKernelParameters kernelParams_;

    //! Scaled virial tensor (6 floats: [XX, XY, XZ, YY, YZ, ZZ])
    std::vector<float> h_virialScaled_;

    /*! \brief Maximum total number of constraints so far.
     *
     * If the new number of constraints is larger then previous maximum, the GPU data arrays are
     * reallocated.
     */
    int numConstraintsThreadsAlloc_;

    /*! \brief Maximum total number of atoms so far.
     *
     * If the new number of atoms is larger then previous maximum, the GPU array with masses is
     * reallocated.
     */
    int numAtomsAlloc_;
};

} // namespace gmx

#endif // GMX_MDLIB_LINCS_CUDA_CUH
