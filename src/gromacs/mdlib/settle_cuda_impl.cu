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

#include "settle_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/settle_cuda.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#include "gromacs/pbcutil/pbc.h"

#if defined(_MSVC)
#include <limits>
#endif

#define GMX_SETTLE_CUDA_TPB 256

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

// Settle Constraints Kernel - see original CPU version in settle.cpp
template <bool updateVelocities, bool computeVirial>
__global__ void settle_kernel(const int                            nSettle,
                              const int3                          *settles,
                              const SettleCuda::SettleParameters   pars,
                              const float3                        *x,
                              float3                              *xprime,
                              const PbcAiuc                        pbcAiuc,
                              float                                invdt,
                              float3                              *v,
                              float*                               virialScaled)
{
    /* ******************************************************************* */
    /*                                                                  ** */
    /*    Original code by Shuichi Miyamoto, last update Oct. 1, 1992   ** */
    /*                                                                  ** */
    /*    Algorithm changes by Berk Hess:                               ** */
    /*    2004-07-15 Convert COM to double precision to avoid drift     ** */
    /*    2006-10-16 Changed velocity update to use differences         ** */
    /*    2012-09-24 Use oxygen as reference instead of COM             ** */
    /*    2016-02    Complete rewrite of the code for SIMD              ** */
    /*                                                                  ** */
    /*    Reference for the SETTLE algorithm                            ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992).    ** */
    /*                                                                  ** */
    /* ******************************************************************* */

    real              almost_zero = real(1e-12);

    int               i = blockIdx.x*blockDim.x+threadIdx.x;

    if (i < nSettle)
    {
        int3            idxes = settles[i];

        float3          x_ow1, x_hw2, x_hw3;

        x_ow1 = x[idxes.x];
        x_hw2 = x[idxes.y];
        x_hw3 = x[idxes.z];

        float3 xprime_ow1, xprime_hw2, xprime_hw3;

        xprime_ow1 = xprime[idxes.x];
        xprime_hw2 = xprime[idxes.y];
        xprime_hw3 = xprime[idxes.z];


        float3 dist21, dist31;
        float3 doh2, doh3;
        float3 sh_hw2, sh_hw3;

        dist21 = pbcDxAiucFloat3(pbcAiuc, x_hw2, x_ow1);
        dist31 = pbcDxAiucFloat3(pbcAiuc, x_hw3, x_ow1);
        doh2   = pbcDxAiucFloat3(pbcAiuc, xprime_hw2, xprime_ow1);

        sh_hw2     = xprime_hw2 - (xprime_ow1 + doh2);
        xprime_hw2 = xprime_hw2 - sh_hw2;

        doh3 =   pbcDxAiucFloat3(pbcAiuc, xprime_hw3, xprime_ow1);

        sh_hw3     = xprime_hw3 - (xprime_ow1 + doh3);
        xprime_hw3 = xprime_hw3 - sh_hw3;

        float3 a1, com;

        a1  = (-doh2 - doh3) * pars.wh;
        com = xprime_ow1 - a1;

        float3 b1;
        b1  = xprime_hw2 - com;

        float3 c1;
        c1  = xprime_hw3 - com;

        float  xakszd = dist21.y * dist31.z - dist21.z * dist31.y;
        float  yakszd = dist21.z * dist31.x - dist21.x * dist31.z;
        float  zakszd = dist21.x * dist31.y - dist21.y * dist31.x;

        float  xaksxd = a1.y * zakszd - a1.z * yakszd;
        float  yaksxd = a1.z * xakszd - a1.x * zakszd;
        float  zaksxd = a1.x * yakszd - a1.y * xakszd;

        float  xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        float  yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        float  zaksyd = xakszd * yaksxd - yakszd * xaksxd;

        float  axlng = rsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        float  aylng = rsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        float  azlng = rsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        float3 trns1, trns2, trns3;

        trns1.x = xaksxd * axlng;
        trns2.x = yaksxd * axlng;
        trns3.x = zaksxd * axlng;
        trns1.y = xaksyd * aylng;
        trns2.y = yaksyd * aylng;
        trns3.y = zaksyd * aylng;
        trns1.z = xakszd * azlng;
        trns2.z = yakszd * azlng;
        trns3.z = zakszd * azlng;


        float2 b0d, c0d;

        b0d.x = trns1.x * dist21.x + trns2.x * dist21.y + trns3.x * dist21.z;
        b0d.y = trns1.y * dist21.x + trns2.y * dist21.y + trns3.y * dist21.z;

        c0d.x = trns1.x * dist31.x + trns2.x * dist31.y + trns3.x * dist31.z;
        c0d.y = trns1.y * dist31.x + trns2.y * dist31.y + trns3.y * dist31.z;

        float  a1d_z;
        float3 b1d, c1d;

        a1d_z = trns1.z * a1.x + trns2.z * a1.y + trns3.z * a1.z;

        b1d.x = trns1.x * b1.x + trns2.x * b1.y + trns3.x * b1.z;
        b1d.y = trns1.y * b1.x + trns2.y * b1.y + trns3.y * b1.z;
        b1d.z = trns1.z * b1.x + trns2.z * b1.y + trns3.z * b1.z;

        c1d.x = trns1.x * c1.x + trns2.x * c1.y + trns3.x * c1.z;
        c1d.y = trns1.y * c1.x + trns2.y * c1.y + trns3.y * c1.z;
        c1d.z = trns1.z * c1.x + trns2.z * c1.y + trns3.z * c1.z;

        float tmp, tmp2;

        float sinphi = a1d_z * rsqrt(pars.ra*pars.ra);
        tmp2     = 1.0 - sinphi * sinphi;

        if (almost_zero > tmp2)
        {
            tmp2 = almost_zero;
        }

        tmp      = rsqrt(tmp2);
        float cosphi = tmp2*tmp;
        float sinpsi = (b1d.z - c1d.z) * pars.irc2 * tmp;
        tmp2     = 1.0 - sinpsi * sinpsi;

        float cospsi = tmp2*rsqrt(tmp2);

        float a2d_y  =  pars.ra * cosphi;
        float b2d_x  = -pars.rc * cospsi;
        float t1     = -pars.rb * cosphi;
        float t2     =  pars.rc * sinpsi * sinphi;
        float b2d_y  =  t1 - t2;
        float c2d_y  =  t1 + t2;

        /*     --- Step3  al,be,ga            --- */
        float alpha  = b2d_x * (b0d.x - c0d.x) + b0d.y * b2d_y + c0d.y * c2d_y;
        float beta   = b2d_x * (c0d.y - b0d.y) + b0d.x * b2d_y + c0d.x * c2d_y;
        float gamma  = b0d.x * b1d.y - b1d.x * b0d.y + c0d.x * c1d.y - c1d.x * c0d.y;
        float al2be2 = alpha * alpha + beta * beta;
        tmp2     = (al2be2 - gamma * gamma);
        float sinthe = (alpha * gamma - beta * tmp2*rsqrt(tmp2)) * rsqrt(al2be2*al2be2);

        /*  --- Step4  A3' --- */
        tmp2     = 1.0 - sinthe * sinthe;
        float  costhe = tmp2*rsqrt(tmp2);

        float3 a3d, b3d, c3d;


        a3d.x  = -a2d_y * sinthe;
        a3d.y  =  a2d_y * costhe;
        a3d.z  =  a1d_z;
        b3d.x  =  b2d_x * costhe - b2d_y * sinthe;
        b3d.y  =  b2d_x * sinthe + b2d_y * costhe;
        b3d.z  =  b1d.z;
        c3d.x  = -b2d_x * costhe - c2d_y * sinthe;
        c3d.y  = -b2d_x * sinthe + c2d_y * costhe;
        c3d.z  =  c1d.z;

        /*    --- Step5  A3 --- */
        float3 a3, b3, c3;

        a3.x = trns1.x*a3d.x + trns1.y*a3d.y + trns1.z*a3d.z;
        a3.y = trns2.x*a3d.x + trns2.y*a3d.y + trns2.z*a3d.z;
        a3.z = trns3.x*a3d.x + trns3.y*a3d.y + trns3.z*a3d.z;

        b3.x = trns1.x*b3d.x + trns1.y*b3d.y + trns1.z*b3d.z;
        b3.y = trns2.x*b3d.x + trns2.y*b3d.y + trns2.z*b3d.z;
        b3.z = trns3.x*b3d.x + trns3.y*b3d.y + trns3.z*b3d.z;

        c3.x = trns1.x*c3d.x + trns1.y*c3d.y + trns1.z*c3d.z;
        c3.y = trns2.x*c3d.x + trns2.y*c3d.y + trns2.z*c3d.z;
        c3.z = trns3.x*c3d.x + trns3.y*c3d.y + trns3.z*c3d.z;


        /* Compute and store the corrected new coordinate */
        xprime_ow1 = com + a3;
        xprime_hw2 = com + b3 + sh_hw2;
        xprime_hw3 = com + c3 + sh_hw3;

        xprime[idxes.x]   = xprime_ow1;
        xprime[idxes.y]   = xprime_hw2;
        xprime[idxes.z]   = xprime_hw3;


        if (updateVelocities || computeVirial)
        {

            float3 da, db, dc;
            da = a3 - a1;
            db = b3 - b1;
            dc = c3 - c1;


            if (updateVelocities)
            {
                float3 v_ow1, v_hw2, v_hw3;

                v_ow1 = v[idxes.x];
                v_hw2 = v[idxes.y];
                v_hw3 = v[idxes.z];

                /* Add the position correction divided by dt to the velocity */
                v_ow1 = da*invdt + v_ow1;
                v_hw2 = db*invdt + v_hw2;
                v_hw3 = dc*invdt + v_hw3;

                v[idxes.x]   = v_ow1;
                v[idxes.y]   = v_hw2;
                v[idxes.z]   = v_hw3;

            }


            if (computeVirial)
            {
                float3 mdo, mdb, mdc;

                mdb = pars.mH*db;
                mdc = pars.mH*dc;
                mdo = pars.mO*da + mdb + mdc;

                atomicAdd(&virialScaled[XX*DIM+XX], -(x_ow1.x*mdo.x + dist21.x*mdb.x + dist31.x*mdc.x));
                atomicAdd(&virialScaled[XX*DIM+YY], -(x_ow1.x*mdo.y + dist21.x*mdb.y + dist31.x*mdc.y));
                atomicAdd(&virialScaled[XX*DIM+ZZ], -(x_ow1.x*mdo.z + dist21.x*mdb.z + dist31.x*mdc.z));

                atomicAdd(&virialScaled[YY*DIM+XX], -(x_ow1.y*mdo.x + dist21.y*mdb.x + dist31.y*mdc.x));
                atomicAdd(&virialScaled[YY*DIM+YY], -(x_ow1.y*mdo.y + dist21.y*mdb.y + dist31.y*mdc.y));
                atomicAdd(&virialScaled[YY*DIM+ZZ], -(x_ow1.y*mdo.z + dist21.y*mdb.z + dist31.y*mdc.z));

                atomicAdd(&virialScaled[ZZ*DIM+XX], -(x_ow1.z*mdo.x + dist21.z*mdb.x + dist31.z*mdc.x));
                atomicAdd(&virialScaled[ZZ*DIM+YY], -(x_ow1.z*mdo.y + dist21.z*mdb.y + dist31.z*mdc.y));
                atomicAdd(&virialScaled[ZZ*DIM+ZZ], -(x_ow1.z*mdo.z + dist21.z*mdb.z + dist31.z*mdc.z));
            }
        }
    }

    return;
}

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
void SettleCuda::Impl::apply(const bool       updateVelocities,
                             const real       invdt,
                             const gmx_bool   computeVirial,
                             tensor           virialScaled)
{

    cudaCheckError();

    int blockSize  = GMX_SETTLE_CUDA_TPB;
    int blockCount = (nSettle + blockSize - 1)/blockSize;

    if (computeVirial)
    {
        cudaMemcpy(virialScaledDevice, virialScaled, DIM*DIM*sizeof(real), cudaMemcpyHostToDevice);
    }
    // Is there a more elegant way of doing this?
    auto kernelPtr = settle_kernel<true, true>;
    if (updateVelocities && computeVirial)
    {
        kernelPtr = settle_kernel<true, true>;
    }
    else if (updateVelocities && !computeVirial)
    {
        kernelPtr = settle_kernel<true, false>;
    }
    else if (!updateVelocities && computeVirial)
    {
        kernelPtr = settle_kernel<false, true>;
    }
    else if (!updateVelocities && !computeVirial)
    {
        kernelPtr = settle_kernel<false, false>;
    }

    /*
       KernelLaunchConfig config;
       config.blockSize[0] = GMX_SETTLE_CUDA_TPB;
       config.blockSize[1] = 1;
       config.blockSize[2] = 1;
       config.gridSize[0]  = (nSettle + GMX_SETTLE_CUDA_TPB - 1)/GMX_SETTLE_CUDA_TPB;
       config.gridSize[1]  = 1;
       config.gridSize[2]  = 1;
       config.stream       = stream;

       const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
        &nSettle,
        &atomIdsDevice,
        &settleParameters,
        &xDevice,
        &xpDevice,
        &pbcAiuc,
        &invdt,
        &vDevice,
        &virialScaledDevice);

       launchGpuKernel(kernelPtr, config, nullptr,
        "settle_kernel<updateVelocities, computeVirial>", kernelArgs);*/

    kernelPtr<<< blockCount, blockSize, 0, stream >>>
    (nSettle,
     atomIdsDevice,
     settleParameters,
     xDevice,
     xpDevice,
     pbcAiuc,
     invdt,
     vDevice,
     virialScaledDevice);

    cudaStreamSynchronize(stream);
    cudaCheckError();

    if (computeVirial)
    {
        cudaMemcpy(virialScaled, virialScaledDevice, DIM*DIM*sizeof(real), cudaMemcpyDeviceToHost);
        cudaCheckError();
    }
    return;
}

/*! \brief Create SETTLE object
 *
 * \param [in] nAtom  Number of atoms that will be handles by SETTLE.
 *                    Used to compute the memory size in allocations and copy.
 * \param [in] mtop   Topology of the system to gen the masses for O and H atoms and
 *                    target O-H and H-H distances. These values are also checked for
 *                    consistency.
 */
SettleCuda::Impl::Impl(const int         nAtom,
                       const gmx_mtop_t &mtop)
    : nAtom(nAtom)
{
    GMX_RELEASE_ASSERT(sizeof(real) == sizeof(float), "Real numbers should be in single precision in GPU code.");
    /*
     * \todo This should be lifted to a separate subroutine that gets the values of Oxygen and Hydrogen
     * masses, checks if they are consistent across the topology and if there is no more than two values
     * for each mass if the free energy perturbation is enabled. In later case, masses may need to be
     * updated on a regular basis (i.e. in set(...) method).
     * \todo Do the checks for FEP
     */
    real mO = -1.0;
    real mH = -1.0;

    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype.at(mt).ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        for (unsigned i = 0; i < iatoms.size()/nral1; i++)
        {
            int3 settler;
            settler.x = iatoms[i*nral1 + 1]; // Oxygen index
            settler.y = iatoms[i*nral1 + 2]; // First hydrogen index
            settler.z = iatoms[i*nral1 + 3]; // Second hydrogen index
            t_atom ow1 = mtop.moltype.at(mt).atoms.atom[settler.x];
            t_atom hw2 = mtop.moltype.at(mt).atoms.atom[settler.y];
            t_atom hw3 = mtop.moltype.at(mt).atoms.atom[settler.z];

            if (mO < 0)
            {
                mO = ow1.m;
            }
            if (mH < 0)
            {
                mH = hw2.m;
            }
            GMX_RELEASE_ASSERT(mO == ow1.m, "Topology has different values for oxygen mass. Should be identical in order to use SETTLE.");
            GMX_RELEASE_ASSERT(hw2.m == hw3.m && hw2.m == mH, "Topology has different values for hydrogen mass. Should be identical in order to use SETTLE.");
        }
    }

    GMX_RELEASE_ASSERT(mO > 0, "Could not find oxygen mass in the topology. Needed in SETTLE.");
    GMX_RELEASE_ASSERT(mH > 0, "Could not find hydrogen mass in the topology. Needed in SETTLE.");

    /*
     * \todo Very similar to SETTLE initialization on CPU. Should be lifted to a separate method
     * (one that gets dOH and dHH values and checks them for consistency)
     */
    int                    settle_type = -1;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int       nral1           = 1 + NRAL(F_SETTLE);
        InteractionList interactionList = mtop.moltype.at(mt).ilist[F_SETTLE];
        for (int i = 0; i < interactionList.size(); i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = interactionList.iatoms[i];
            }
            else if (interactionList.iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }

    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    real dOH = mtop.ffparams.iparams[settle_type].settle.doh;
    real dHH = mtop.ffparams.iparams[settle_type].settle.dhh;

    //Impl(nAtom, mO, mH, dOH, dHH);

    initSettleParameters(&settleParameters, mO, mH, dOH, dHH);

    cudaMalloc(&xDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&xpDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&vDevice, nAtom*DIM*sizeof(float));

    cudaMalloc(&virialScaledDevice, DIM*DIM*sizeof(float));
    cudaStreamCreate(&stream);

    cudaCheckError();

    //Impl(nAtom, 15.99940, 1.00800, 0.09572, 0.15139);
}

/*! \brief Create SETTLE object
 *
 * \param [in] nAtom  Number of atoms that will be handles by SETTLE.
 *                    Used to compute the memory size in allocations and copy.
 * \param [in] mO     Mass of the oxygen atom.
 * \param [in] mH     Mass of the hydrogen atom.
 * \param [in] dOH    Target distance for O-H bonds.
 * \param [in] dHH    Target for the distance between two hydrogen atoms.
 */
SettleCuda::Impl::Impl(const int  nAtom,
                       const real mO,  const real mH,
                       const real dOH, const real dHH)
    : nAtom(nAtom)
{
    GMX_RELEASE_ASSERT(sizeof(real) == sizeof(float), "Real numbers should be in single precision in GPU code.");

    initSettleParameters(&settleParameters, mO, mH, dOH, dHH);

    cudaMalloc(&xDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&xpDevice, nAtom*DIM*sizeof(float));
    cudaMalloc(&vDevice, nAtom*DIM*sizeof(float));

    cudaMalloc(&virialScaledDevice, DIM*DIM*sizeof(float));
    cudaStreamCreate(&stream);

    cudaCheckError();
}

SettleCuda::Impl::~Impl()
{
}


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
void SettleCuda::Impl::set(const t_idef               &idef,
                           const t_mdatoms gmx_unused &md)
{
    const int  nral1     = 1 + NRAL(F_SETTLE);
    t_ilist    il_settle = idef.il[F_SETTLE];
    t_iatom   *iatoms    = il_settle.iatoms;
    nSettle   = il_settle.nr/nral1;
    if ((int)atomIdsHost.size() < nSettle)
    {
        if (atomIdsHost.size() != 0)
        {
            cudaFree(atomIdsDevice);
            cudaCheckError();
        }
        atomIdsHost.resize(nSettle);
        cudaMalloc(&atomIdsDevice, nSettle*sizeof(int3));
        cudaCheckError();
    }
    for (int i = 0; i < nSettle; i++)
    {
        int3 settler;
        settler.x         = iatoms[i*nral1 + 1]; // Oxygen index
        settler.y         = iatoms[i*nral1 + 2]; // First hydrogen index
        settler.z         = iatoms[i*nral1 + 3]; // Second hydrogen index
        atomIdsHost.at(i) = settler;
    }
    cudaMemcpy(atomIdsDevice, atomIdsHost.data(), nSettle*sizeof(int3), cudaMemcpyHostToDevice);
    cudaCheckError();
}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] *pbc The PBC data in t_pbc format.
 */
void SettleCuda::Impl::setPbc(t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc);
}

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
void SettleCuda::Impl::copyCoordinatesToGpu(const rvec * x, const rvec * xp)
{
    cudaMemcpy(xDevice, x, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaMemcpy(xpDevice, xp, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
    cudaCheckError();
    cudaDeviceSynchronize();
}

/*! \brief
 * Copy velocities from provided CPU location to GPU.
 *
 * Nothing is done if the argument provided is a nullptr.
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  CPU pointer where velocities should be copied from.
 */
void SettleCuda::Impl::copyVelocitiesToGpu(const rvec * v)
{
    if (v != nullptr)
    {
        cudaMemcpy(vDevice, v, nAtom*sizeof(float3), cudaMemcpyHostToDevice);
        cudaCheckError();
        cudaDeviceSynchronize();
    }
}

/*! \brief
 * Copy coordinates from GPU to provided CPU location.
 *
 * Copies the constrained coordinates to the provided location. The coordinates
 * are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to.
 */
void SettleCuda::Impl::copyCoordinatesFromGpu(rvec * xp)
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
void SettleCuda::Impl::copyVelocitiesFromGpu(rvec * v)
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
void SettleCuda::Impl::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    this->xDevice  = (float3*)xDevice;
    this->xpDevice = (float3*)xpDevice;
    this->vDevice  = (float3*)vDevice;
}


SettleCuda::SettleCuda(const int         nAtom,
                       const gmx_mtop_t &mtop)
    : impl_(new Impl(nAtom, mtop))
{
}

SettleCuda::SettleCuda(const int nAtom,
                       const real mO,  const real mH,
                       const real dOH, const real dHH)
    : impl_(new Impl(nAtom, mO, mH, dOH, dHH))
{
}

SettleCuda::~SettleCuda() = default;

void SettleCuda::apply(bool       updateVelocities,
                       real       invdt,
                       gmx_bool   bCalcVir,
                       tensor     virialScaled)
{
    impl_->apply(updateVelocities,
                 invdt,
                 bCalcVir,
                 virialScaled);
}

void SettleCuda::setPbc(t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void SettleCuda::set(const t_idef    &idef,
                     const t_mdatoms &md)
{
    impl_->set(idef, md);
}

void SettleCuda::copyCoordinatesToGpu(const rvec *x, const rvec *xp)
{
    impl_->copyCoordinatesToGpu(x, xp);
}

void SettleCuda::copyVelocitiesToGpu(const rvec *v)
{
    impl_->copyVelocitiesToGpu(v);
}

void SettleCuda::copyCoordinatesFromGpu(rvec * xp)
{
    impl_->copyCoordinatesFromGpu(xp);
}

void SettleCuda::copyVelocitiesFromGpu(rvec * v)
{
    impl_->copyVelocitiesFromGpu(v);
}

void SettleCuda::setXVPointers(rvec * xDevice, rvec * xpDevice, rvec * vDevice)
{
    impl_->setXVPointers(xDevice, xpDevice, vDevice);
}
