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
#ifndef GMX_MDLIB_SETTLE_CUDA_IMPL_H
#define GMX_MDLIB_SETTLE_CUDA_IMPL_H


#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/settle_cuda.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"

/*! \internal \brief Class with interfaces and data for CUDA version of SETTLE. */
class SettleCuda::Impl
{

    public:
        /*! \brief Create SETTLE object
         *
         *  Extracts masses for oxygen and hydrogen as well as the O-H and H-H target distances
         *  from the topology data (mtop), check their values for consistensy and calls the
         *  following constructor.
         *
         * \param [in] nAtom  Number of atoms that will be handles by SETTLE.
         *                    Used to compute the memory size in allocations and copy.
         * \param [in] mtop   Topology of the system to gen the masses for O and H atoms and
         *                    target O-H and H-H distances. These values are also checked for
         *                    consistency.
         */
        Impl(const int          nAtom,
             const gmx_mtop_t  &mtop);

        /*! \brief Create SETTLE object
         *
         * \param [in] nAtom  Number of atoms that will be handles by SETTLE.
         *                    Used to compute the memory size in allocations and copy.
         * \param [in] mO     Mass of the oxygen atom.
         * \param [in] mH     Mass of the hydrogen atom.
         * \param [in] dOH    Target distance for O-H bonds.
         * \param [in] dHH    Target for the distance between two hydrogen atoms.
         */
        Impl(const int nAtom,
             const real mO,  const real mH,
             const real dOH, const real dHH);

        /*! \brief Destructor.*/
        ~Impl();

        /*! \brief Apply SETTLE.
         *
         * Applies SETTLE to coordinates and velocities, stored on GPU.
         * Data at pointers xPrime and v (class fields) change in the GPU
         * memory. The results are not automatically copied back to the CPU
         * memory. Method uses this class data structures which should be
         * updated when needed using update method.
         *
         * \param[in] updateVelocities  If the velocities should be constrained.
         * \param[in] invdt             Inversed timestep (to scale Lagrange
         *                              multipliers when velocities are updated)
         * \param[in] computeVirial     If virial should be updated.
         * \param[in,out] virialScaled  Scaled virial tensor to be updated.
         */
        void apply(bool       updateVelocities,
                   real       invdt,
                   bool       computeVirial,
                   tensor     virialScaled);

        /*! \brief
         * Update data-structures (e.g. after NB search step).
         *
         * Updates the constraints data and copies it to the GPU. Should be
         * called if the particles were sorted, redistributed between domains, etc.
         * Does not recycle the data preparation routines from the CPU version.
         * All three atoms from the single water molecule should be handled by the same GPU.
         *
         * SETTLEs atom ID's is taken from idef.il[F_SETTLE].iatoms.
         *
         * \param [in] idef    System topology
         * \param [in] md      Atoms data. Can be used to update masses if needed (not used now).
         */
        void set(const t_idef         &idef,
                 const t_mdatoms      &md);

        /*! \brief
         * Update PBC data.
         *
         * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
         *
         * \param[in] *pbc The PBC data in t_pbc format.
         */
        void setPbc(t_pbc *pbc);

        /*! \brief
         * Copy coordinates from provided CPU location to GPU.
         *
         * Copies the coordinates before the integration step (x) and coordinates
         * after the integration step (xp) from the provided CPU location to GPU.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *x  CPU pointer where coordinates should be copied from.
         * \param[in] *xp CPU pointer where coordinates should be copied from.
         */
        void copyCoordinatesToGpu(const rvec * x, const rvec * xp);

        /*! \brief
         * Copy velocities from provided CPU location to GPU.
         *
         * Nothing is done if the argument provided is a nullptr.
         * The data are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *v  CPU pointer where velocities should be copied from.
         */
        void copyVelocitiesToGpu(const rvec * v);

        /*! \brief
         * Copy coordinates from GPU to provided CPU location.
         *
         * Copies the constrained coordinates to the provided location. The coordinates
         * are assumed to be in float3/fvec format (single precision).
         *
         * \param[out] *xp CPU pointer where coordinates should be copied to.
         */
        void copyCoordinatesFromGpu(rvec * xp);

        /*! \brief
         * Copy velocities from GPU to provided CPU location.
         *
         * The velocities are assumed to be in float3/fvec format (single precision).
         *
         * \param[in] *v  Pointer to velocities data.
         */
        void copyVelocitiesFromGpu(rvec * v);

        /*! \brief
         * Set the internal GPU-memory x, xprime and v pointers.
         *
         * Data is not copied. The data are assumed to be in float3/fvec format
         * (float3 is used internally, but the data layout should be identical).
         *
         * \param[in] *xDevice  Pointer to the coordinates before integrator update (on GPU)
         * \param[in] *xpDevice Pointer to the coordinates after integrator update, before update (on GPU)
         * \param[in] *vDevice  Pointer to the velocities before integrator update (on GPU)
         */
        void setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice);

    private:

        cudaStream_t       stream;                       //!< CUDA stream

        PbcAiuc            pbcAiuc;                      //!< Periodic boundary data

        int                nAtom;                        //!< Number of atoms

        float3            *xDevice;                      //!< Coordinates before the timestep (in the GPU memory)
        float3            *xpDevice;                     //!< Coordinates after the timestep, before constraining (on GPU).
                                                         //   These will be changed upon constraining.
        float3            *vDevice;                      //!< Velocities of atoms (on GPU)

        real              *virialScaledHost;             //!< Scaled virial tensor (9 reals, GPU)
        real              *virialScaledDevice;           //!< Scaled virial tensor (9 reals, GPU)

        int                nSettle;

        std::vector<int3>  atomIdsHost;                      //!< Indexes of atoms (.x for oxygen, .y and.z for hydrogens (CPU)
        int3              *atomIdsDevice;                    //!< Indexes of atoms (.x for oxygen, .y and.z for hydrogens (GPU)

        SettleParameters   settleParameters;
};


//! Initializes a projection matrix.
static void init_proj_matrix(real invmO, real invmH, real dOH, real dHH,
                             matrix inverseCouplingMatrix)
{
    /* We normalize the inverse masses with invmO for the matrix inversion.
     * so we can keep using masses of almost zero for frozen particles,
     * without running out of the float range in invertMatrix.
     */
    double invmORelative = 1.0;
    double invmHRelative = invmH/static_cast<double>(invmO);
    double distanceRatio = dHH/static_cast<double>(dOH);

    /* Construct the constraint coupling matrix */
    matrix mat;
    mat[0][0] = invmORelative + invmHRelative;
    mat[0][1] = invmORelative*(1.0 - 0.5*gmx::square(distanceRatio));
    mat[0][2] = invmHRelative*0.5*distanceRatio;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = invmHRelative + invmHRelative;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    gmx::invertMatrix(mat, inverseCouplingMatrix);

    msmul(inverseCouplingMatrix, 1/invmO, inverseCouplingMatrix);
}

//! Initializes settle parameters.
static void initSettleParameters(SettleCuda::SettleParameters *p,
                                 real mO,  real mH,
                                 real dOH, real dHH)
{
    /* We calculate parameters in double precision to minimize errors.
     * The velocity correction applied during SETTLE coordinate constraining
     * introduces a systematic error of approximately 1 bit per atom,
     * depending on what the compiler does with the code.
     */
    double wohh;

    real   invmO = 1.0/mO;
    real   invmH = 1.0/mH;

    p->mO     = mO;
    p->mH     = mH;
    wohh      = mO + 2.0*mH;
    p->wh     = mH/wohh;
    p->dOH    = dOH;
    p->dHH    = dHH;
    double rc = dHH/2.0;
    double ra = 2.0*mH*std::sqrt(dOH*dOH - rc*rc)/wohh;
    p->rb     = std::sqrt(dOH*dOH - rc*rc) - ra;
    p->rc     = rc;
    p->ra     = ra;
    p->irc2   = 1.0/dHH;

    /* For projection: inverse masses and coupling matrix inversion */
    p->imO    = invmO;
    p->imH    = invmH;

    p->invdOH = 1.0/dOH;
    p->invdHH = 1.0/dHH;

    init_proj_matrix(invmO, invmH, dOH, dHH, p->invmat);
}

#endif
