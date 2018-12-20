/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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


#ifndef GMX_MDLIB_UPDATE_CONSTRAINTS_H
#define GMX_MDLIB_UPDATE_CONSTRAINTS_H


#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/mdlib/nbnxn_cuda/gpuD2DConstraintsCUDA.h"
#include "gromacs/mdtypes/group.h"


struct settleparam_t
{
    real   mO;
    real   mH;
    real   wh;
    real   dOH;
    real   dHH;
    real   ra;
    real   rb;
    real   rc;
    real   irc2;
    /* For projection */
    real   imO;
    real   imH;
    real   invdOH;
    real   invdHH;
    matrix invmat;
};


struct settledata
{
    settleparam_t massw;    /* Parameters for SETTLE for coordinates */
    settleparam_t mass1;    /* Parameters with all masses 1, for forces */

    int           nsettle;  /* The number of settles on our rank */
    int          *ow1;      /* Index to OW1 atoms, size nsettle + SIMD padding */
    int          *hw2;      /* Index to HW2 atoms, size nsettle + SIMD padding */
    int          *hw3;      /* Index to HW3 atoms, size nsettle + SIMD padding */
    real         *virfac;   /* Virial factor 0 or 1, size nsettle + SIMD pad. */
    int           nalloc;   /* Allocation size of ow1, hw2, hw3, virfac */

    bool          bUseSimd; /* Use SIMD intrinsics code, if possible */
};


enum class NumTempScaleValues
{
    single,       //!< Single T-scaling value (either one group or all values =1)
    multiple      //!< Multiple T-scaling values, need to use T-group indices
};





//! Unit of work within LINCS.
struct Task
{
    //! First constraint for this task.
    int    b0 = 0;
    //! b1-1 is the last constraint for this task.
    int    b1 = 0;
    //! The number of constraints in triangles.
    int    ntriangle = 0;
    //! The list of triangle constraints.
    int   *triangle = nullptr;
    //! The bits tell if the matrix element should be used.
    int   *tri_bits = nullptr;
    //! Allocation size of triangle and tri_bits.
    int    tri_alloc = 0;
    //! Number of indices.
    int    nind = 0;
    //! Constraint index for updating atom data.
    int   *ind = nullptr;
    //! Number of indices.
    int    nind_r = 0;
    //! Constraint index for updating atom data.
    int   *ind_r = nullptr;
    //! Allocation size of ind and ind_r.
    int    ind_nalloc = 0;
    //! Temporary variable for virial calculation.
    tensor vir_r_m_dr = {{0}};
    //! Temporary variable for lambda derivative.
    real   dhdlambda;
};




class Lincs
{
    public:
        //! The global number of constraints.
        int             ncg = 0;
        //! The global number of flexible constraints.
        int             ncg_flex = 0;
        //! The global number of constraints in triangles.
        int             ncg_triangle = 0;
        //! The number of iterations.
        int             nIter = 0;
        //! The order of the matrix expansion.
        int             nOrder = 0;
        //! The maximum number of constraints connected to a single atom.
        int             max_connect = 0;

        //! The number of real constraints.
        int             nc_real = 0;
        //! The number of constraints including padding for SIMD.
        int             nc = 0;
        //! The number we allocated memory for.
        int             nc_alloc = 0;
        //! The number of constraint connections.
        int             ncc = 0;
        //! The number we allocated memory for.
        int             ncc_alloc = 0;
        //! The FE lambda value used for filling blc and blmf.
        real            matlam = 0;
        //! mapping from topology to LINCS constraints.
        int            *con_index = nullptr;
        //! The reference distance in topology A.
        real           *bllen0 = nullptr;
        //! The reference distance in top B - the r.d. in top A.
        real           *ddist = nullptr;
        //! The atom pairs involved in the constraints.
        int            *bla = nullptr;
        //! 1/sqrt(invmass1  invmass2).
        real           *blc = nullptr;
        //! As blc, but with all masses 1.
        real           *blc1 = nullptr;
        //! Index into blbnb and blmf.
        int            *blnr = nullptr;
        //! List of constraint connections.
        int            *blbnb = nullptr;
        //! The local number of constraints in triangles.
        int             ntriangle = 0;
        //! The number of constraint connections in triangles.
        int             ncc_triangle = 0;
        //! Communicate before each LINCS interation.
        bool            bCommIter = false;
        //! Matrix of mass factors for constraint connections.
        real           *blmf = nullptr;
        //! As blmf, but with all masses 1.
        real           *blmf1 = nullptr;
        //! The reference bond length.
        real           *bllen = nullptr;
        //! The local atom count per constraint, can be NULL.
        int            *nlocat = nullptr;

        /*! \brief The number of tasks used for LINCS work.
         *
         * \todo This is mostly used to loop over \c task, which would
         * be nicer to do with range-based for loops, but the thread
         * index is used for constructing bit masks and organizing the
         * virial output buffer, so other things need to change,
         * first. */
        int               ntask = 0;
        /*! \brief LINCS thread division */
        std::vector<Task> task;
        //! Atom flags for thread parallelization.
        gmx_bitmask_t    *atf = nullptr;
        //! Allocation size of atf
        int               atf_nalloc = 0;
        //! Are the LINCS tasks interdependent?
        bool              bTaskDep = false;
        //! Are there triangle constraints that cross task borders?
        bool              bTaskDepTri = false;
        //! Arrays for temporary storage in the LINCS algorithm.
        /*! @{ */
        rvec           *tmpv   = nullptr;
        real           *tmpncc = nullptr;
        real           *tmp1   = nullptr;
        real           *tmp2   = nullptr;
        real           *tmp3   = nullptr;
        real           *tmp4   = nullptr;
        /*! @} */
        //! The Lagrange multipliers times -1.
        real               *mlambda = nullptr;
        //! Storage for the constraint RMS relative deviation output.
        std::array<real, 2> rmsdData = {{0}};
};



GPU_FUNC_QUALIFIER
void do_lincs_gpu(void *lincsd_ptr,  t_pbc *pbc, rvec *x, rvec *xp, rvec* v,  tensor vir_r_m_dr, real *invmass, const real invdt, gmx_bool bCalcVir, const t_commrec *cr, matrix box) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void updateMDLeapfrogSimple_gpu(int                       nrend,
                                real                      dt,
                                real                      dtPressureCouple,
                                const rvec * gmx_restrict invMassPerDim,
                                const t_grp_tcstat      * tcstat,
                                const unsigned short    * cTC,
                                const rvec * gmx_restrict x,
                                rvec       * gmx_restrict xprime,
                                rvec       * gmx_restrict v,
                                const rvec * gmx_restrict f,
                                NumTempScaleValues        numTempScaleValues) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void settle_gpu(void *settledptr,
                const t_pbc* pbc,
                const real *x, real *xprime,
                real invdt, real * gmx_restrict v,
                tensor vir_r_m_dr,
                bool *bErrorHasOccurred, bool bCorrectVelocity, bool bCalcVirial) GPU_FUNC_TERM


GPU_FUNC_QUALIFIER
void gpuUpdateConstraintsSetSize(int xsize) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
int gpuUpdateConstraintsGetSize() GPU_FUNC_TERM_WITH_RETURN(0)

GPU_FUNC_QUALIFIER
void gpuUpdateConstraintsSetTimestepInfo(bool bNS, bool bNSNextStep,
                                         bool copybackVelocity) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuUpdateConstraintsSetGpuNB(gmx_nbnxn_gpu_t* gpu_nbv) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuUpdateConstraintsCopyXPToXOnDevice() GPU_FUNC_TERM

#if defined(__CUDACC__)
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"


enum {
    epbcdxRECTANGULAR = 1, epbcdxTRICLINIC,
    epbcdx2D_RECT,       epbcdx2D_TRIC,
    epbcdx1D_RECT,       epbcdx1D_TRIC,
    epbcdxSCREW_RECT,    epbcdxSCREW_TRIC,
    epbcdxNOPBC,         epbcdxUNSUPPORTED
};


//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
        cudaError_t e = cudaGetLastError();                                 \
        if (e != cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));           \
            exit(0); \
        }                                                                 \
}


/* maths operations */
/* imported from cpu versions in math/vec.h */
__forceinline__ __device__
void svmul_gpu(real a, const rvec v1, rvec v2)
{
    v2[XX] = a*v1[XX];
    v2[YY] = a*v1[YY];
    v2[ZZ] = a*v1[ZZ];
}


__forceinline__ __device__
void rvec_add_gpu(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
void ivec_add_gpu(const ivec a, const ivec b, ivec c)
{
    int x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
void rvec_inc_atomic(rvec a, const rvec b)
{
    atomicAdd(&a[XX], b[XX]);
    atomicAdd(&a[YY], b[YY]);
    atomicAdd(&a[ZZ], b[ZZ]);
}

__forceinline__ __device__
void rvec_inc_gpu(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__
void rvec_dec_atomic(rvec a, const rvec b)
{
/* real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;*/

    atomicAdd(&a[XX], -1.0f*b[XX]);
    atomicAdd(&a[YY], -1.0f*b[YY]);
    atomicAdd(&a[ZZ], -1.0f*b[ZZ]);
}

__forceinline__ __device__
void rvec_dec_gpu(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__
void cprod_gpu(const rvec a, const rvec b, rvec c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

__forceinline__ __device__
real iprod_gpu(const rvec a, const rvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

__forceinline__ __device__
real norm_gpu(const rvec a)
{
    return sqrt(iprod_gpu(a, a));
}

__forceinline__ __device__
real gmx_angle_gpu(const rvec a, const rvec b)
{
    rvec w;
    real wlen, s;

    cprod_gpu(a, b, w);

    wlen  = norm_gpu(w);
    s     = iprod_gpu(a, b);

    return atan2f(wlen, s); //requires float
}

__forceinline__ __device__
void clear_ivec_gpu(ivec a)
{
    a[XX] = 0;
    a[YY] = 0;
    a[ZZ] = 0;
}
__forceinline__ __device__
void rvec_sub_gpu(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
real norm2_gpu(const rvec a)
{
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

__forceinline__ __device__
void copy_rvec_gpu(const rvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__
void copy_ivec_gpu(const ivec a, ivec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__
real cos_angle_gpu(const rvec a, const rvec b)
{
    /*
     *                  ax*bx + ay*by + az*bz
     * cos-vec (a,b) =  ---------------------
     *                      ||a|| * ||b||
     */
    real   cosval;
    int    m;
    double aa, bb, ip, ipa, ipb, ipab; /* For accuracy these must be double! */

    ip = ipa = ipb = 0.0;
    for (m = 0; (m < DIM); m++) /* 18 */
    {
        aa   = a[m];
        bb   = b[m];
        ip  += aa*bb;
        ipa += aa*aa;
        ipb += bb*bb;
    }
    ipab = ipa*ipb;
    if (ipab > 0)
    {
        cosval = ip*rsqrt(ipab); /*  7 */  //double precision
    }
    else
    {
        cosval = 1;
    }
    /* 25 TOTAL */
    if (cosval > 1.0)
    {
        return 1.0;
    }
    if (cosval < -1.0)
    {
        return -1.0;
    }

    return cosval;
}



__device__ static inline float invsqrt(float x)
{
    return 1.0f/std::sqrt(x);
}


__device__ static inline void unitv_gpu(const rvec src, rvec dest)
{
    real linv;

    linv     = invsqrt(norm2_gpu(src));
    dest[XX] = linv*src[XX];
    dest[YY] = linv*src[YY];
    dest[ZZ] = linv*src[ZZ];
}




__device__ static inline void pbc_dx_aiuc_gpu(const t_pbc *pbc, const rvec x1, const rvec x2, rvec dx, const rvec pbc_hbox_diag, const matrix pbc_box, const rvec pbc_mhbox_diag, const rvec  pbc_fbox_diag,  const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift  )
{


    int  i, j;
    rvec dx_start, trial;
    real d2min, d2trial;
    ivec ishift, ishift_start;

    rvec_sub_gpu(x1, x2, dx);
    clear_ivec_gpu(ishift);


    switch (pbc->ePBCDX)
    {
        case epbcdxRECTANGULAR:
            for (i = 0; i < DIM; i++)
            {
                if (dx[i] > pbc_hbox_diag[i])
                {
                    dx[i] -=  pbc_fbox_diag[i];
                    ishift[i]--;
                }
                else if (dx[i] <= pbc_mhbox_diag[i])
                {
                    dx[i] +=  pbc_fbox_diag[i];
                    ishift[i]++;
                }
            }
            break;

        case epbcdxTRICLINIC:

            for (i = DIM-1; i >= 1; i--)
            {
                if (dx[i] > pbc_hbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] -= pbc_box[i][j];
                    }
                    ishift[i]--;
                }
                else if (dx[i] <= pbc_mhbox_diag[i])
                {
                    for (j = i; j >= 0; j--)
                    {
                        dx[j] += pbc_box[i][j];
                    }
                    ishift[i]++;
                }
            }

            /* Allow 2 shifts in x */
            if (dx[XX] > pbc_hbox_diag[XX])
            {
                dx[XX] -= pbc_fbox_diag[XX];
                ishift[XX]--;
                if (dx[XX] > pbc_hbox_diag[XX])
                {
                    dx[XX] -= pbc_fbox_diag[XX];
                    ishift[XX]--;
                }
            }
            else if (dx[XX] <= pbc_mhbox_diag[XX])
            {
                dx[XX] += pbc_fbox_diag[XX];
                ishift[XX]++;
                if (dx[XX] <= pbc_mhbox_diag[XX])
                {
                    dx[XX] += pbc_fbox_diag[XX];
                    ishift[XX]++;
                }
            }

            /* dx is the distance in a rectangular box */
            d2min = norm2_gpu(dx);
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec_gpu(dx, dx_start);
                copy_ivec_gpu(ishift, ishift_start);
                d2min = norm2_gpu(dx);
                /* Now try all possible shifts, when the distance is within max_cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add_gpu(dx_start, pbc_tric_vec[i], trial);
                    d2trial = norm2_gpu(trial);
                    if (d2trial < d2min)
                    {
                        copy_rvec_gpu(trial, dx);
                        ivec_add_gpu(ishift_start, pbc_tric_shift[i], ishift);
                        d2min = d2trial;
                    }
                    i++;
                }
            }

            break;
        case epbcdx2D_RECT:
            for (i = 0; i < DIM; i++)
            {
                if (i != pbc->dim)
                {
                    if (dx[i] > pbc_hbox_diag[i])
                    {
                        dx[i] -= pbc_fbox_diag[i];
                        ishift[i]--;
                    }
                    else if (dx[i] <= pbc_mhbox_diag[i])
                    {
                        dx[i] += pbc_fbox_diag[i];
                        ishift[i]++;
                    }
                }
            }
            break;
        case epbcdx2D_TRIC:
            d2min = 0;
            for (i = DIM-1; i >= 1; i--)
            {
                if (i != pbc->dim)
                {
                    if (dx[i] > pbc_hbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] -= pbc_box[i][j];
                        }
                        ishift[i]--;
                    }
                    else if (dx[i] <= pbc_mhbox_diag[i])
                    {
                        for (j = i; j >= 0; j--)
                        {
                            dx[j] += pbc_box[i][j];
                        }
                        ishift[i]++;
                    }
                    d2min += dx[i]*dx[i];
                }
            }
            if (pbc->dim != XX)
            {
                /* Allow 2 shifts in x */
                if (dx[XX] > pbc_hbox_diag[XX])
                {
                    dx[XX] -= pbc_fbox_diag[XX];
                    ishift[XX]--;
                    if (dx[XX] > pbc_hbox_diag[XX])
                    {
                        dx[XX] -= pbc_fbox_diag[XX];
                        ishift[XX]--;
                    }
                }
                else if (dx[XX] <= pbc_mhbox_diag[XX])
                {
                    dx[XX] += pbc_fbox_diag[XX];
                    ishift[XX]++;
                    if (dx[XX] <= pbc_mhbox_diag[XX])
                    {
                        dx[XX] += pbc_fbox_diag[XX];
                        ishift[XX]++;
                    }
                }
                d2min += dx[XX]*dx[XX];
            }
            if (d2min > pbc->max_cutoff2)
            {
                copy_rvec_gpu(dx, dx_start);
                copy_ivec_gpu(ishift, ishift_start);
                /* Now try all possible shifts, when the distance is within max\
                   _cutoff
                 * it must be the shortest possible distance.
                 */
                i = 0;
                while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
                {
                    rvec_add_gpu(dx_start, pbc_tric_vec[i], trial);
                    d2trial = 0;
                    for (j = 0; j < DIM; j++)
                    {
                        if (j != pbc->dim)
                        {
                            d2trial += trial[j]*trial[j];
                        }
                    }
                    if (d2trial < d2min)
                    {
                        copy_rvec_gpu(trial, dx);
                        ivec_add_gpu(ishift_start, pbc_tric_shift[i], ishift);
                        d2min = d2trial;
                    }
                    i++;
                }
            }
            break;
        case epbcdxNOPBC:
        case epbcdxUNSUPPORTED:
            break;


    }
    //  is = IVEC2IS(ishift);
    return;


}



#endif
#endif
