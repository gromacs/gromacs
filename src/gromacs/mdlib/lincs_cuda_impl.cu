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
#include "gmxpre.h"

#include "lincs_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs_cuda.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#include "gromacs/pbcutil/pbc.h"

#if defined(_MSVC)
#include <limits>
#endif

#define GMX_LINCS_CUDA_TPB 256


/*! \brief Main kernel for LINCS constraints.
 *
 * \todo Combine arguments
 * \todo Move everything to local/shared memory, try to get rid of atomics.
 * \todo Template updateVelocities and virial.
 *
 * \param[in]     ncons                     Total number of the constrain threads (empty spaces at the end of thread blocks included)
 * \param[in]     x                         Coordinates before the timestep
 * \param[in,out] xp                        Coordinates after the timestep. Will be updated to satisfy the constraints by this kernel.
 * \param[in]     nIter                     Number of iterations used to compute inverse matrix.
 * \param[in]     nOrder                    LINCS projection order for correcting the direction of constraint.
 * \param[in]     constraints               List of constraints.
 * \param[in]     constraintsR0             List of target distances for constraints.
 * \param[in]     coupledConstraintsCounts  Numbers of coupled constraints for each constraint
 * \param[in]     coupledConstraintsIdxes   Indexes of coupled constraints
 * \param[in]     massFactors               Mass factors: ( (+/-) * (1/sqrt(1/m1 + 1/m2)) * (1/m2) * 1/sqrt(1/m2 + 1/m3)),
 *                                          where m1 and m3 are coupled through m2 and sign + or - indicates the order,
 *                                          in which they are arranged in atoms array.)
 * \param[in]     matrixA                   Place to store constraints matrix.
 * \param[in]     pbcAiuc                   Periodic boundary information
 * \param[in]     mlambda                   Lagrange multiplier times mass, saved by this kernel for virial evaluation.
 * \param[in]     updateVelocities          If the velocities should be computed.
 * \param[in,out] v                         Velocities to update.
 * \param[in]     invdt                     Inverse timestep (needed to update velocities).
 */
__global__ void lincs_kernel(const int            ncons,
                             const float3        *x,
                             float3              *xp,
                             const int            nIter,
                             const int            nOrder,
                             const int2          *constraints,
                             const real          *constraintsR0,
                             const int           *coupledConstraintsCounts,
                             const int           *coupledConstraintsIdxes,
                             const real          *massFactors,
                             real                *matrixA,
                             const PbcAiuc        pbcAiuc,
                             real                *mlambda,
                             const real          *invmass,
                             const bool           updateVelocities,
                             float3              *v,
                             const real           invdt)
{

    int c           = blockIdx.x*blockDim.x+threadIdx.x;
    int cs          = threadIdx.x;
    int blockStarts = blockIdx.x*blockDim.x;

    extern __shared__ float3 r[];
    extern __shared__ float  rhs[];

    if (c < ncons)
    {

        float mvb;

        int2  pair = constraints[c];
        int   i    = pair.x;
        int   j    = pair.y;

        real  ml = 0.0;

        if (i != -1)
        {
            float rlen0 = constraintsR0[c];

            float im1      = invmass[i];
            float im2      = invmass[j];

            float sqrtmu = rsqrt(im1 + im2);


            float3 xi = x[i];
            float3 xj = x[j];

            float3 dx   = pbcDxAiucFloat3(pbcAiuc, xi, xj);
            float  rlen = rsqrtf(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);

            float3 rc = rlen*dx;
            r[cs] = rc;

            xi = xp[i];
            xj = xp[j];
            dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);

            mvb = sqrtmu*((rc.x*dx.x + rc.y*dx.y + rc.z*dx.z) - rlen0);

            float sol  = mvb;

            __syncthreads();

            int coupledConstraintsCount = coupledConstraintsCounts[c];

            for (int n = 0; n < coupledConstraintsCount; n++)
            {
                int    index = n*ncons + c;
                int    c1    = coupledConstraintsIdxes[index]; //Can be moved to local/shared memory

                float3 rc1 = r[c1-blockStarts];
                matrixA[index] = massFactors[index]*(rc.x*rc1.x + rc.y*rc1.y + rc.z*rc1.z);

            }

            rhs[cs] = mvb;
            __syncthreads();

            for (int rec = 0; rec < nOrder; rec++)
            {
                mvb = 0;

                for (int n = 0; n < coupledConstraintsCount; n++)
                {
                    int index = n*ncons + c;
                    int c1    = coupledConstraintsIdxes[index];

                    mvb = mvb + matrixA[index]*rhs[c1-blockStarts + blockDim.x*(rec % 2)];

                }
                rhs[cs + blockDim.x*((rec + 1) % 2)] = mvb;
                sol  = sol + mvb;
                __syncthreads();

            }

            ml = sqrtmu*sol;

            mvb      = ml;

            float3 tmp     = rc*mvb;

            atomicAdd(&xp[i], -tmp*im1);
            atomicAdd(&xp[j], tmp*im2);

            __syncthreads();

            for (int iter = 0; iter < nIter; iter++)
            {
                float len2, dlen2;

                xi = xp[i];
                xj = xp[j];


                dx = pbcDxAiucFloat3(pbcAiuc, xi, xj);

                len2  = rlen0*rlen0;
                dlen2 = 2.0f*len2 - norm2(dx);

                if (dlen2 > 0)
                {
                    mvb = sqrtmu*(rlen0 - dlen2*rsqrt(dlen2));
                }
                else
                {
                    mvb = sqrtmu*rlen0;
                }

                rhs[cs]  = mvb;
                sol      = mvb;
                __syncthreads();

                for (int rec = 0; rec < nOrder; rec++)
                {
                    mvb = 0;

                    for (int n = 0; n < coupledConstraintsCount; n++)
                    {
                        int index = n*ncons + c;
                        int c1    = coupledConstraintsIdxes[index];

                        mvb = mvb + matrixA[index]*rhs[c1-blockStarts + blockDim.x*(rec % 2)];

                    }
                    rhs[cs + blockDim.x*((rec + 1) % 2)] = mvb;
                    sol  = sol + mvb;
                    __syncthreads();

                }

                mvb         = sqrtmu*sol;
                float sqrtmu_sol  = mvb;
                ml += mvb;

                mvb      = sqrtmu_sol;

                float3 tmp = rc*mvb;

                atomicAdd(&xp[i], -tmp*im1);
                atomicAdd(&xp[j], tmp*im2);
                __syncthreads();
            }

            if (updateVelocities)
            {
                mvb      = invdt*ml;

                float3 tmp     = rc*mvb;

                atomicAdd(&v[i], -tmp*im1);
                atomicAdd(&v[j], tmp*im2);
            }

        }
        mlambda[c] = ml; // Needed for virial
    }

    return;
}

/*! \brief Kernel to compute virial.
 *
 * \todo Make it into a templated version of the main kernel?
 *
 * \param[in]     ncons                     Total number of the constrain threads (empty spaces at the end of thread blocks included)
 * \param[in]     x                         Coordinates before the timestep
 * \param[in,out] virialScaled              Scaled virial tensor to be updated.
 * \param[in]     mlambda                   Lagrange multiplier times mass, saved by this kernel for virial evaluation.
 * \param[in]     constraints               List of constraints.
 * \param[in]     constraintsR0             List of target distances for constraints.
 * \param[in]     pbcAiuc                   Periodic boundary information
 */
__global__ void lincs_virial_kernel(const int             ncons,
                                    const float3         *x,
                                    real                 *virialScaled,
                                    const real           *mlambda,
                                    const int2           *constraints,
                                    const real           *constraintsR0,
                                    const PbcAiuc         pbcAiuc)
{


    int c = blockIdx.x*blockDim.x+threadIdx.x;

    if (c < ncons)
    {
        float tmp0, tmp1;

        int2  pair = constraints[c];
        int   i    = pair.x;
        int   j    = pair.y;

        if (i != -1)
        {

            float3 xi = x[i];
            float3 xj = x[j];

            float3 dx   = pbcDxAiucFloat3(pbcAiuc, xi, xj);
            float  rlen = rsqrtf(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);

            rvec   rc;

            rc[0] = rlen*dx.x;
            rc[1] = rlen*dx.y;
            rc[2] = rlen*dx.z;

            tmp0 = -constraintsR0[c]*mlambda[c];
            for (int d1 = 0; d1 < DIM; d1++)
            {
                tmp1 = tmp0*rc[d1];
                for (int d2 = 0; d2 < DIM; d2++)
                {
                    atomicAdd(&(virialScaled[d1*DIM+d2]), -tmp1*rc[d2]);
                }
            }
        }
    }
    return;
}

/*
 * Temporary solution.
 */
#define cudaCheckError() {                                          \
        cudaError_t e = cudaGetLastError();                                 \
        if (e != cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));           \
            exit(0); \
        }                                                                 \
}

/*! \brief Create LINCS object
 *
 * \param [in] nAtom  Number of atoms that will be handles by LINCS.
 *                    Used to compute the memory size in allocations and copy.
 * \param [in] nIter  Number of iterations used to compute inverse matrix.
 * \param [in] nOrder LINCS projection order for correcting the direction of constraint.
 */
LincsCuda::Impl::Impl(int nAtom,
                      int nIter,
                      int nOrder)
    : nAtom(nAtom), nIter(nIter), nOrder(nOrder)

{
    GMX_ASSERT(sizeof(real) == sizeof(float), "Real numbers should be in single precision in GPU code.");
    cudaMalloc(&xDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&xpDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&vDevice, nAtom*DIM*sizeof(float));

    cudaMalloc(&virialScaledDevice, DIM*DIM*sizeof(float));
    maxConstraintsNumberSoFar = 0;
    cudaStreamCreate(&stream);
    cudaCheckError();
}

LincsCuda::Impl::~Impl()
{
}


/*! \brief Apply LINCS.
 *
 * Applies LINCS to coordinates and velocities, stored on GPU.
 * Data at pointers xPrime and v (class fields) change in the GPU
 * memory. The results are not automatically copied back to the CPU
 * memory. Method uses this class data structures which should be
 * updated when needed using update method.
 *
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] invdt             Inversed timestep (to scale Lagrange
 *                              multipliers when velocities are updated)
 * \param[in] bCalcVir          If virial should be updated.
 * \param[in] scaleLambda       If the Lagrange multipliers should be scaled
 *                              before virial is computed.
 * \param[in,out] virialScaled  Scaled virial tensor to be updated.
 */
void LincsCuda::Impl::apply(const bool       updateVelocities,
                            const real       invdt,
                            const gmx_bool   bCalcVir,
                            tensor           virialScaled)
{

    cudaCheckError();

    int blockSize  = GMX_LINCS_CUDA_TPB;
    int blockCount = (nConstraintsThreads + blockSize - 1)/blockSize;

    /*KernelLaunchConfig config;
       config.blockSize[0]     = blockSize;
       config.blockSize[1]     = 0;
       config.blockSize[2]     = 0;
       config.gridSize[0]      = (nConstraintsThreads + blockSize - 1)/blockSize;
       config.sharedMemorySize = blockSize*DIM*sizeof(float);
       config.stream           = stream;

       auto kernelPtr = lincs_kernel;
       const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
            &nConstraintsThreads, &xDevice, &xpDevice,
            &nIter, &nOrder,
            &constraintsDevice, &constraintsR0Device,
            &coupledConstraintsCountsDevice, &coupledConstraintsIdxesDevice,
            &massFactorsDevice, &matrixADevice,
            &pbcAiuc,
            &mlambdaDevice,
            &invmassDevice,
            &updateVelocities, &vDevice, &invdt);
       launchGpuKernel(lincs_kernel, config, nullptr, "lincs", kernelArgs);*/

    lincs_kernel
    <<< blockCount, blockSize, blockSize*DIM*sizeof(float), stream>>>
    (nConstraintsThreads, xDevice, xpDevice,
     nIter, nOrder,
     constraintsDevice, constraintsR0Device,
     coupledConstraintsCountsDevice, coupledConstraintsIdxesDevice,
     massFactorsDevice, matrixADevice,
     pbcAiuc,
     mlambdaDevice,
     invmassDevice,
     updateVelocities, vDevice, invdt);

    cudaCheckError();


    if (bCalcVir)
    {
        cudaMemcpy(virialScaledDevice, virialScaled, DIM*DIM*sizeof(real), cudaMemcpyHostToDevice);

        /*config.sharedMemorySize = 0;
           kernelPtr = lincs_virial_kernel;
           kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
            &nConstraintsThreads, &xDevice,
            &virialScaledDevice, &mlambdaDevice,
            &constraintsDevice, &constraintsR0Device,
            &pbcAiuc);
           launchGpuKernel(lincs_virial_kernel, config, nullptr, "lincs-virial", kernelArgs);*/

        lincs_virial_kernel
        <<< blockCount, blockSize, 0, stream>>>
        (nConstraintsThreads, xDevice, virialScaledDevice, mlambdaDevice, constraintsDevice, constraintsR0Device, pbcAiuc);

        cudaCheckError();

        cudaMemcpy(virialScaled, virialScaledDevice, DIM*DIM*sizeof(real), cudaMemcpyDeviceToHost);

        cudaCheckError();
    }

    cudaStreamSynchronize(stream);
    cudaCheckError();

    return;
}

/*! \brief Helper function to go through constraints recurrently
 *
 *  Counts the total number constraints, connected to an atom (including those, connected through other constraints).
 */
inline int countCoupled(int a, std::vector<int> *spaceNeeded,
                        std::vector<std::vector<std::tuple<int, int, int> > > *atomsAdjacencyList)

{
    int c2, a2, sign;
    int counted = 0;
    for (unsigned i = 0; i < atomsAdjacencyList->at(a).size(); i++)
    {
        std::tie(a2, c2, sign) = atomsAdjacencyList->at(a).at(i);
        if (spaceNeeded->at(c2) == -1)
        {
            spaceNeeded->at(c2) = 0; // To indicate we've been here
            counted            += 1 + countCoupled(a2, spaceNeeded, atomsAdjacencyList);
        }
    }
    return counted;
}

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
void LincsCuda::Impl::set(const t_idef    &idef,
                          const t_mdatoms &md)
{

    int blockSize = GMX_LINCS_CUDA_TPB;

    //t_idef idef = top.idef;
    t_iatom  *iatoms      = idef.il[F_CONSTR].iatoms;
    const int nConstraint = idef.il[F_CONSTR].nr/3;
    // Constructing adjacency list --- usefull intermediate structure
    std::vector<std::vector<std::tuple<int, int, int> > > atomsAdjacencyList(nAtom);
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];

        // Each constraint will be represented as a tuple, containing index of the second constrained atom,
        // index of the constraint and a sign that indicates the order of atoms in which they are listed.
        // Sign is needed to compute the mass factors.
        atomsAdjacencyList.at(a1).push_back(std::make_tuple(a2, c, +1));
        atomsAdjacencyList.at(a2).push_back(std::make_tuple(a1, c, -1));
    }

    // Compute, how many coupled constraints are in front of each constraint.
    // Needed to introduce splits in data so that all coupled constraints will be computed in a single GPU block.
    // The position 'c' of the vector spaceNeeded should have the number of constraints that are coupled to a constraint
    // 'c' and are after 'c' in the vector. Only first index of the connected group of the constraints is needed later in the
    // code, hence the spaceNeeded vetor is also used to keep track if the constrain was already counted.
    std::vector<int> spaceNeeded;
    spaceNeeded.resize(nConstraint, -1);
    std::fill(spaceNeeded.begin(), spaceNeeded.end(), -1);
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];
        if (spaceNeeded.at(c) == -1)
        {
            spaceNeeded.at(c) = countCoupled(a1, &spaceNeeded, &atomsAdjacencyList) +
                countCoupled(a2, &spaceNeeded, &atomsAdjacencyList);
        }

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        /*spaceNeeded.at(c) = atomsAdjacencyList.at(a1).size() + atomsAdjacencyList.at(a2).size() - 2;

           // Only space needed in front of the constraint should be counted.
           while (spaceNeeded.at(c) > 0 && c + 1 < nConstraint)
           {
            c++;
            spaceNeeded.at(c) = spaceNeeded.at(c-1) - 1;
           }*/
    }

    // Map of splits in the constraints data. For each 'old' constraint index gives 'new' which
    // takes into account the empty spaces which might be needed in the end of each thread block.
    std::vector<int> splitMap;
    splitMap.resize(nConstraint, -1);
    int              currentMapIndex = 0;
    for (int c = 0; c < nConstraint; c++)
    {
        if (currentMapIndex / blockSize != (currentMapIndex + spaceNeeded.at(c)) / blockSize)
        {
            currentMapIndex = ((currentMapIndex/blockSize) + 1) * blockSize;
        }
        splitMap.at(c) = currentMapIndex;
        currentMapIndex++;
    }
    nConstraintsThreads = currentMapIndex + blockSize - currentMapIndex % blockSize;


    // Initialize constraints and their target indexes taking into account the splits in the
    // data arrays.
    int2 pair;
    pair.x = -1;
    pair.y = -1;
    constraintsHost.resize(nConstraintsThreads, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsR0Host.resize(nConstraintsThreads, 0.0);
    std::fill(constraintsR0Host.begin(), constraintsR0Host.end(), 0.0);
    for (int c = 0; c < nConstraint; c++)
    {
        int  a1     = iatoms[3*c + 1];
        int  a2     = iatoms[3*c + 2];
        int  type   = iatoms[3*c];

        int2 pair;
        pair.x = a1;
        pair.y = a2;
        constraintsHost.at(splitMap.at(c))   = pair;
        constraintsR0Host.at(splitMap.at(c)) = idef.iparams[type].constr.dA;

    }

    // The adjacency list of constraints (i.e. the list of coupled constraints for each constraint).
    // We map a single thread to a single constraint, hence each thread 'c' will be using one element from
    // coupledConstraintsCountsHost array, which is the number of constraints coupled to the constraint 'c'.
    // The coupled constraints indexes are placed into the coupledConstraintsIdxesHost array. Latter is organized
    // as a one-dimensional array to ensure good memory alignment. It is addressed as [c + i*nConstraintsThreads],
    // where 'i' goes from zero to the number of constraints coupled to 'c'. 'nConstraintsThreads' is the width of
    // the array --- a number, greater then total number of constraints, taking into account the splits in the
    // constraints array due to the GPU block borders. This number can be adjusted to improve memory access pattern.
    // Mass factors are saved in a similar data structure.
    int              maxCoupledConstraints = 0;
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        int nCoupedConstraints = atomsAdjacencyList.at(a1).size() + atomsAdjacencyList.at(a2).size() - 2;

        if (nCoupedConstraints > maxCoupledConstraints)
        {
            maxCoupledConstraints = nCoupedConstraints;
        }
    }

    coupledConstraintsCountsHost.resize(nConstraintsThreads, 0);
    coupledConstraintsIdxesHost.resize(maxCoupledConstraints*nConstraintsThreads, -1);
    massFactorsHost.resize(maxCoupledConstraints*nConstraintsThreads, -1);

    for (int c1 = 0; c1 < nConstraint; c1++)
    {
        coupledConstraintsCountsHost.at(splitMap.at(c1))  = 0;
        int c1a1     = iatoms[3*c1 + 1];
        int c1a2     = iatoms[3*c1 + 2];
        int c2;
        int c2a1;
        int c2a2;

        int sign;

        c2a1 = c1a1;
        for (unsigned j = 0; j < atomsAdjacencyList.at(c1a1).size(); j++)
        {

            std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a1).at(j);

            if (c1 != c2)
            {
                int index = nConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIdxesHost.at(index) = splitMap.at(c2);

                int  center = c1a1;

                real sqrtmu1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                real sqrtmu2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = -sign*md.invmass[center]*sqrtmu1*sqrtmu2;

                coupledConstraintsCountsHost.at(splitMap.at(c1))++;

            }
        }

        c2a1 = c1a2;
        for (unsigned j = 0; j < atomsAdjacencyList.at(c1a2).size(); j++)
        {

            std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a2).at(j);

            if (c1 != c2)
            {
                int index = nConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIdxesHost.at(index) = splitMap.at(c2);

                int  center = c1a2;

                real sqrtmu1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                real sqrtmu2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = sign*md.invmass[center]*sqrtmu1*sqrtmu2;

                coupledConstraintsCountsHost.at(splitMap.at(c1))++;

            }
        }
    }

    if (nConstraint > maxConstraintsNumberSoFar)
    {

        if (maxConstraintsNumberSoFar > 0)
        {
            cudaFree(invmassDevice);
            cudaFree(mlambdaDevice);

            cudaFree(constraintsDevice);
            cudaFree(constraintsR0Device);

            cudaFree(coupledConstraintsCountsDevice);
            cudaFree(coupledConstraintsIdxesDevice);
            cudaFree(massFactorsDevice);
            cudaFree(matrixADevice);

        }
        maxConstraintsNumberSoFar = nConstraint;

        cudaMalloc(&invmassDevice, nAtom*sizeof(real));
        cudaMalloc(&mlambdaDevice, nConstraintsThreads*sizeof(float));

        cudaMalloc(&constraintsDevice, nConstraintsThreads*sizeof(int2));
        cudaMalloc(&constraintsR0Device, nConstraintsThreads*sizeof(float));

        cudaMalloc(&coupledConstraintsCountsDevice, nConstraintsThreads*sizeof(int));
        cudaMalloc(&coupledConstraintsIdxesDevice, maxCoupledConstraints*nConstraintsThreads*sizeof(int));
        cudaMalloc(&massFactorsDevice, maxCoupledConstraints*nConstraintsThreads*sizeof(float));
        cudaMalloc(&matrixADevice, maxCoupledConstraints*nConstraintsThreads*sizeof(float));

    }

    cudaMemcpy(constraintsDevice, constraintsHost.data(), nConstraintsThreads*sizeof(int2), cudaMemcpyHostToDevice);
    cudaMemcpy(constraintsR0Device, constraintsR0Host.data(), nConstraintsThreads*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(coupledConstraintsCountsDevice, coupledConstraintsCountsHost.data(),
               nConstraintsThreads*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(coupledConstraintsIdxesDevice, coupledConstraintsIdxesHost.data(),
               maxCoupledConstraints*nConstraintsThreads*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(massFactorsDevice, massFactorsHost.data(), maxCoupledConstraints*nConstraintsThreads*sizeof(float), cudaMemcpyHostToDevice);

    cudaCheckError();

    GMX_ASSERT(md.invmass != nullptr, "Masses of attoms should be specified.\n");
    cudaMemcpy(invmassDevice, md.invmass, nAtom*sizeof(real), cudaMemcpyHostToDevice);

    cudaCheckError();

}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] *pbc The PBC data in t_pbc format.
 */
void LincsCuda::Impl::setPbc(t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc);
}

/*! \brief
 * Copy coordinates and velocities from provided CPU location to GPU.
 *
 * Copies the coordinates before the integration step (x), coordinates
 * after the integration step (xp) and velocities (v) from the provided
 * CPU location to GPU. The data are assumed to be in float3/fvec format
 * (single precision).
 *
 * \param[in] *x  CPU pointer where coordinates should be copied from.
 * \param[in] *xp CPU pointer where coordinates should be copied from.
 * \param[in] *v  CPU pointer where velocities should be copied from.
 */
void LincsCuda::Impl::copyCoordinatesToGpu(const rvec * x, const rvec * xp, const rvec * v)
{
    cudaMemcpy(xDevice, x, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaMemcpy(xpDevice, xp, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    if (v != nullptr)
    {
        cudaMemcpy(vDevice, v, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    }

    cudaCheckError();
    cudaDeviceSynchronize();
}

/*! \brief
 * Copy coordinates from GPU to provided CPU location.
 *
 * Copies the constrained coordinates to the provided location. The coordinates
 * are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to.
 */
void LincsCuda::Impl::copyCoordinatesFromGpu(rvec * xp)
{
    cudaMemcpy(xp, xpDevice, nAtom*sizeof(float3), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaCheckError();
}

/*! \brief
 * Copy velocities from GPU to provided CPU location.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  Pointer to velocities data.
 */
void LincsCuda::Impl::copyVelocitiesFromGpu(rvec * v)
{
    cudaMemcpy(v, vDevice, nAtom*sizeof(float3), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaCheckError();
}

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
void LincsCuda::Impl::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    this->xDevice  = (float3*)xDevice;
    this->xpDevice = (float3*)xpDevice;
    this->vDevice  = (float3*)vDevice;
}


LincsCuda::LincsCuda(int nAtom,
                     int nIter,
                     int nOrder)
    : impl_(new Impl(nAtom, nIter, nOrder))
{
}

LincsCuda::~LincsCuda() = default;

void LincsCuda::apply(bool       updateVelocities,
                      real       invdt,
                      gmx_bool   bCalcVir,
                      tensor     virialScaled)
{
    impl_->apply(updateVelocities,
                 invdt,
                 bCalcVir,
                 virialScaled);
}

void LincsCuda::setPbc(t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void LincsCuda::set(const t_idef    &idef,
                    const t_mdatoms &md)
{
    impl_->set(idef, md);
}

void LincsCuda::copyCoordinatesToGpu(const rvec *x, const rvec *xp, const rvec *v)
{
    impl_->copyCoordinatesToGpu(x, xp, v);
}

void LincsCuda::copyCoordinatesFromGpu(rvec * xp)
{
    impl_->copyCoordinatesFromGpu(xp);
}

void LincsCuda::copyVelocitiesFromGpu(rvec * v)
{
    impl_->copyVelocitiesFromGpu(v);
}

void LincsCuda::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    impl_->setXVPointers(xDevice, xpDevice, vDevice);
}
