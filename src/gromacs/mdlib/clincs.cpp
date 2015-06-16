/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc-simd.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

/* MSVC 2010 produces buggy SIMD PBC code, disable SIMD for MSVC <= 2010 */
#if defined GMX_SIMD_HAVE_REAL && !(defined _MSC_VER && _MSC_VER < 1700) && !defined(__ICL)
#define LINCS_SIMD
#endif


#if defined(GMX_SIMD_X86_AVX_256) || defined(GMX_SIMD_X86_AVX2_256)

// This was originally work-in-progress for augmenting the SIMD module with
// masked load/store operations. Instead, that turned into and extended SIMD
// interface that supports gather/scatter in all platforms, which will be
// part of a future Gromacs version. However, since the code for bonded
// interactions and LINCS was already written it would be a pity not to get
// the performance gains in Gromacs-5.1. For this reason we have added it as
// a bit of a hack in the two files that use it. It will be replaced with the
// new generic functionality after version 5.1

#    ifdef GMX_DOUBLE
static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose4_r(gmx_simd_double_t *row0,
                           gmx_simd_double_t *row1,
                           gmx_simd_double_t *row2,
                           gmx_simd_double_t *row3)
{
    __m256d tmp0, tmp1, tmp2, tmp3;

    tmp0  = _mm256_unpacklo_pd(*row0, *row1);
    tmp2  = _mm256_unpacklo_pd(*row2, *row3);
    tmp1  = _mm256_unpackhi_pd(*row0, *row1);
    tmp3  = _mm256_unpackhi_pd(*row2, *row3);
    *row0 = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    *row1 = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    *row2 = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
    *row3 = _mm256_permute2f128_pd(tmp1, tmp3, 0x31);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd4_transpose_to_simd_r(const gmx_simd4_double_t *a,
                                   gmx_simd_double_t        *row0,
                                   gmx_simd_double_t        *row1,
                                   gmx_simd_double_t        *row2,
                                   gmx_simd_double_t        *row3)
{
    *row0 = a[0];
    *row1 = a[1];
    *row2 = a[2];
    *row3 = a[3];

    gmx_hack_simd_transpose4_r(row0, row1, row2, row3);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose_to_simd4_r(gmx_simd_double_t   row0,
                                   gmx_simd_double_t   row1,
                                   gmx_simd_double_t   row2,
                                   gmx_simd_double_t   row3,
                                   gmx_simd4_double_t *a)
{
    a[0] = row0;
    a[1] = row1;
    a[2] = row2;
    a[3] = row3;

    gmx_hack_simd_transpose4_r(&a[0], &a[1], &a[2], &a[3]);
}


#    ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#        define gmx_hack_simd4_load3_r(mem)      _mm256_maskload_pd((mem), _mm_castsi128_ps(_mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1)))
#        define gmx_hack_simd4_store3_r(mem, x)   _mm256_maskstore_pd((mem), _mm_castsi128_ps(_mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1)), (x))
#    else
#        define gmx_hack_simd4_load3_r(mem)      _mm256_maskload_pd((mem), _mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1))
#        define gmx_hack_simd4_store3_r(mem, x)   _mm256_maskstore_pd((mem), _mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1), (x))
#    endif

#    else /* single instead of double */
static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose4_r(gmx_simd_float_t *row0,
                           gmx_simd_float_t *row1,
                           gmx_simd_float_t *row2,
                           gmx_simd_float_t *row3)
{
    __m256 tmp0, tmp1, tmp2, tmp3;

    tmp0  = _mm256_unpacklo_ps(*row0, *row1);
    tmp2  = _mm256_unpacklo_ps(*row2, *row3);
    tmp1  = _mm256_unpackhi_ps(*row0, *row1);
    tmp3  = _mm256_unpackhi_ps(*row2, *row3);
    *row0 = _mm256_shuffle_ps(tmp0, tmp2, 0x44);
    *row1 = _mm256_shuffle_ps(tmp0, tmp2, 0xEE);
    *row2 = _mm256_shuffle_ps(tmp1, tmp3, 0x44);
    *row3 = _mm256_shuffle_ps(tmp1, tmp3, 0xEE);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd4_transpose_to_simd_r(const gmx_simd4_float_t *a,
                                   gmx_simd_float_t        *row0,
                                   gmx_simd_float_t        *row1,
                                   gmx_simd_float_t        *row2,
                                   gmx_simd_float_t        *row3)
{
    *row0 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[0]), a[4], 1);
    *row1 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[1]), a[5], 1);
    *row2 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[2]), a[6], 1);
    *row3 = _mm256_insertf128_ps(_mm256_castps128_ps256(a[3]), a[7], 1);

    gmx_hack_simd_transpose4_r(row0, row1, row2, row3);
}

static gmx_inline void gmx_simdcall
gmx_hack_simd_transpose_to_simd4_r(gmx_simd_float_t   row0,
                                   gmx_simd_float_t   row1,
                                   gmx_simd_float_t   row2,
                                   gmx_simd_float_t   row3,
                                   gmx_simd4_float_t *a)
{
    gmx_hack_simd_transpose4_r(&row0, &row1, &row2, &row3);

    a[0] = _mm256_extractf128_ps(row0, 0);
    a[1] = _mm256_extractf128_ps(row1, 0);
    a[2] = _mm256_extractf128_ps(row2, 0);
    a[3] = _mm256_extractf128_ps(row3, 0);
    a[4] = _mm256_extractf128_ps(row0, 1);
    a[5] = _mm256_extractf128_ps(row1, 1);
    a[6] = _mm256_extractf128_ps(row2, 1);
    a[7] = _mm256_extractf128_ps(row3, 1);
}
#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#        define gmx_hack_simd4_load3_r(mem)      _mm_maskload_ps((mem), _mm_castsi256_pd(_mm_set_epi32(0, -1, -1, -1)))
#        define gmx_hack_simd4_store3_r(mem, x)   _mm_maskstore_ps((mem), _mm_castsi256_pd(_mm_set_epi32(0, -1, -1, -1)), (x))
#else
#        define gmx_hack_simd4_load3_r(mem)      _mm_maskload_ps((mem), _mm_set_epi32(0, -1, -1, -1))
#        define gmx_hack_simd4_store3_r(mem, x)   _mm_maskstore_ps((mem), _mm_set_epi32(0, -1, -1, -1), (x))
#endif

#endif /* double */

#endif /* AVX */

#ifdef GMX_SIMD_HAVE_REAL
/*! \brief Store differences between indexed rvecs in SIMD registers.
 *
 * Returns SIMD register with the difference vectors:
 *     v[pair_index[i*2]] - v[pair_index[i*2 + 1]]
 *
 * \param[in]     v           Array of rvecs
 * \param[in]     pair_index  Index pairs for GMX_SIMD_REAL_WIDTH vector pairs
 * \param[in,out] buf         Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    dx          SIMD register with x difference
 * \param[out]    dy          SIMD register with y difference
 * \param[out]    dz          SIMD register with z difference
 */
static gmx_inline void gmx_simdcall
gmx_hack_simd_gather_rvec_dist_pair_index(const rvec      *v,
                                          const int       *pair_index,
                                          real gmx_unused *buf,
                                          gmx_simd_real_t *dx,
                                          gmx_simd_real_t *dy,
                                          gmx_simd_real_t *dz)
{
#if defined(GMX_SIMD_X86_AVX_256) || defined(GMX_SIMD_X86_AVX2_256)
    int              i;
    gmx_simd4_real_t d[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  tmp;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        d[i] = gmx_simd4_sub_r(gmx_hack_simd4_load3_r(&(v[pair_index[i*2 + 0]][0])),
                               gmx_hack_simd4_load3_r(&(v[pair_index[i*2 + 1]][0])));
    }

    gmx_hack_simd4_transpose_to_simd_r(d, dx, dy, dz, &tmp);
#else
#if GMX_ALIGNMENT
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) buf_aligned[3*GMX_SIMD_REAL_WIDTH];
#else
    real* buf_aligned = buf;
#endif

    int i, m;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        /* Store the distances packed and aligned */
        for (m = 0; m < DIM; m++)
        {
            buf_aligned[m*GMX_SIMD_REAL_WIDTH + i] =
                v[pair_index[i*2]][m] - v[pair_index[i*2 + 1]][m];
        }
    }
    *dx = gmx_simd_load_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH);
    *dy = gmx_simd_load_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH);
    *dz = gmx_simd_load_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH);
#endif
}


/*! \brief Stores SIMD vector into multiple rvecs.
 *
 * \param[in]     x           SIMD register with x-components of the vectors
 * \param[in]     y           SIMD register with y-components of the vectors
 * \param[in]     z           SIMD register with z-components of the vectors
 * \param[in,out] buf         Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    v           Array of GMX_SIMD_REAL_WIDTH rvecs
 */
static gmx_inline void gmx_simdcall
gmx_simd_store_vec_to_rvec(gmx_simd_real_t  x,
                           gmx_simd_real_t  y,
                           gmx_simd_real_t  z,
                           real gmx_unused *buf,
                           rvec            *v)
{
#if defined(GMX_SIMD_X86_AVX_256) || defined(GMX_SIMD_X86_AVX2_256)
    int              i;
    gmx_simd4_real_t s4[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  zero = gmx_simd_setzero_r();

    gmx_hack_simd_transpose_to_simd4_r(x, y, z, zero, s4);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        gmx_hack_simd4_store3_r(v[i], s4[i]);
    }
#else
#if GMX_ALIGNMENT
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) buf_aligned[3*GMX_SIMD_REAL_WIDTH];
#else
    real* buf_aligned = buf;
#endif

    int i, m;

    gmx_simd_store_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH, x);
    gmx_simd_store_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH, y);
    gmx_simd_store_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH, z);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        for (m = 0; m < DIM; m++)
        {
            v[i][m] = buf_aligned[m*GMX_SIMD_REAL_WIDTH + i];
        }
    }
#endif
}
#endif /* GMX_SIMD_HAVE_REAL */


typedef struct {
    int    b0;         /* first constraint for this task */
    int    b1;         /* b1-1 is the last constraint for this task */
    int    ntriangle;  /* the number of constraints in triangles */
    int   *triangle;   /* the list of triangle constraints */
    int   *tri_bits;   /* the bits tell if the matrix element should be used */
    int    tri_alloc;  /* allocation size of triangle and tri_bits */
    int    nind;       /* number of indices */
    int   *ind;        /* constraint index for updating atom data */
    int    nind_r;     /* number of indices */
    int   *ind_r;      /* constraint index for updating atom data */
    int    ind_nalloc; /* allocation size of ind and ind_r */
    tensor vir_r_m_dr; /* temporary variable for virial calculation */
    real   dhdlambda;  /* temporary variable for lambda derivative */
    real  *simd_buf;   /* Aligned work array for SIMD */
} lincs_task_t;

typedef struct gmx_lincsdata {
    int             ncg;          /* the global number of constraints */
    int             ncg_flex;     /* the global number of flexible constraints */
    int             ncg_triangle; /* the global number of constraints in triangles */
    int             nIter;        /* the number of iterations */
    int             nOrder;       /* the order of the matrix expansion */
    int             max_connect;  /* the maximum number of constrains connected to a single atom */

    int             nc_real;      /* the number of real constraints */
    int             nc;           /* the number of constraints including padding for SIMD */
    int             nc_alloc;     /* the number we allocated memory for */
    int             ncc;          /* the number of constraint connections */
    int             ncc_alloc;    /* the number we allocated memory for */
    real            matlam;       /* the FE lambda value used for filling blc and blmf */
    int            *con_index;    /* mapping from topology to LINCS constraints */
    real           *bllen0;       /* the reference distance in topology A */
    real           *ddist;        /* the reference distance in top B - the r.d. in top A */
    int            *bla;          /* the atom pairs involved in the constraints */
    real           *blc;          /* 1/sqrt(invmass1 + invmass2) */
    real           *blc1;         /* as blc, but with all masses 1 */
    int            *blnr;         /* index into blbnb and blmf */
    int            *blbnb;        /* list of constraint connections */
    int             ntriangle;    /* the local number of constraints in triangles */
    int             ncc_triangle; /* the number of constraint connections in triangles */
    gmx_bool        bCommIter;    /* communicate before each LINCS interation */
    real           *blmf;         /* matrix of mass factors for constraint connections */
    real           *blmf1;        /* as blmf, but with all masses 1 */
    real           *bllen;        /* the reference bond length */
    int            *nlocat;       /* the local atom count per constraint, can be NULL */

    int             ntask;        /* The number of tasks = #threads for LINCS */
    lincs_task_t   *task;         /* LINCS thread division */
    gmx_bitmask_t  *atf;          /* atom flags for thread parallelization */
    int             atf_nalloc;   /* allocation size of atf */
    gmx_bool        bTaskDep;     /* are the LINCS tasks interdependent? */
    /* arrays for temporary storage in the LINCS algorithm */
    rvec           *tmpv;
    real           *tmpncc;
    real           *tmp1;
    real           *tmp2;
    real           *tmp3;
    real           *tmp4;
    real           *mlambda; /* the Lagrange multipliers * -1 */
    /* storage for the constraint RMS relative deviation output */
    real            rmsd_data[3];
} t_gmx_lincsdata;

/* Define simd_width for memory allocation used for SIMD code */
#ifdef LINCS_SIMD
static const int simd_width = GMX_SIMD_REAL_WIDTH;
#else
static const int simd_width = 1;
#endif
/* We can't use small memory alignments on many systems, so use min 64 bytes*/
static const int align_bytes = std::max<int>(64, simd_width*sizeof(real));

real *lincs_rmsd_data(struct gmx_lincsdata *lincsd)
{
    return lincsd->rmsd_data;
}

real lincs_rmsd(struct gmx_lincsdata *lincsd, gmx_bool bSD2)
{
    if (lincsd->rmsd_data[0] > 0)
    {
        return sqrt(lincsd->rmsd_data[bSD2 ? 2 : 1]/lincsd->rmsd_data[0]);
    }
    else
    {
        return 0;
    }
}

/* Do a set of nrec LINCS matrix multiplications.
 * This function will return with up to date thread-local
 * constraint data, without an OpenMP barrier.
 */
static void lincs_matrix_expand(const struct gmx_lincsdata *lincsd,
                                const lincs_task_t *li_task,
                                const real *blcc,
                                real *rhs1, real *rhs2, real *sol)
{
    int        b0, b1, nrec, rec;
    const int *blnr  = lincsd->blnr;
    const int *blbnb = lincsd->blbnb;

    b0   = li_task->b0;
    b1   = li_task->b1;
    nrec = lincsd->nOrder;

    for (rec = 0; rec < nrec; rec++)
    {
        int b;

        if (lincsd->bTaskDep)
        {
#pragma omp barrier
        }
        for (b = b0; b < b1; b++)
        {
            real mvb;
            int  n;

            mvb = 0;
            for (n = blnr[b]; n < blnr[b+1]; n++)
            {
                mvb = mvb + blcc[n]*rhs1[blbnb[n]];
            }
            rhs2[b] = mvb;
            sol[b]  = sol[b] + mvb;
        }

        real *swap;

        swap = rhs1;
        rhs1 = rhs2;
        rhs2 = swap;
    } /* nrec*(ncons+2*nrtot) flops */

    if (lincsd->ntriangle > 0)
    {
        /* Perform an extra nrec recursions for only the constraints
         * involved in rigid triangles.
         * In this way their accuracy should come close to those of the other
         * constraints, since traingles of constraints can produce eigenvalues
         * around 0.7, while the effective eigenvalue for bond constraints
         * is around 0.4 (and 0.7*0.7=0.5).
         */

        if (lincsd->bTaskDep)
        {
            /* We need a barrier here, since other threads might still be
             * reading the contents of rhs1 and/o rhs2.
             * We could avoid this barrier by introducing two extra rhs
             * arrays for the triangle constraints only.
             */
#pragma omp barrier
        }

        /* Constraints involved in a triangle are ensured to be in the same
         * LINCS task. This means no barriers are required during the extra
         * iterations for the triangle constraints.
         */
        const int *triangle = li_task->triangle;
        const int *tri_bits = li_task->tri_bits;

        for (rec = 0; rec < nrec; rec++)
        {
            int tb;

            for (tb = 0; tb < li_task->ntriangle; tb++)
            {
                int  b, bits, nr0, nr1, n;
                real mvb;

                b    = triangle[tb];
                bits = tri_bits[tb];
                mvb  = 0;
                nr0  = blnr[b];
                nr1  = blnr[b+1];
                for (n = nr0; n < nr1; n++)
                {
                    if (bits & (1 << (n - nr0)))
                    {
                        mvb = mvb + blcc[n]*rhs1[blbnb[n]];
                    }
                }
                rhs2[b] = mvb;
                sol[b]  = sol[b] + mvb;
            }

            real *swap;

            swap = rhs1;
            rhs1 = rhs2;
            rhs2 = swap;
        } /* nrec*(ntriangle + ncc_triangle*2) flops */
    }
}

static void lincs_update_atoms_noind(int ncons, const int *bla,
                                     real prefac,
                                     const real *fac, rvec *r,
                                     const real *invmass,
                                     rvec *x)
{
    int  b, i, j;
    real mvb, im1, im2, tmp0, tmp1, tmp2;

    if (invmass != NULL)
    {
        for (b = 0; b < ncons; b++)
        {
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            im1      = invmass[i];
            im2      = invmass[j];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0*im1;
            x[i][1] -= tmp1*im1;
            x[i][2] -= tmp2*im1;
            x[j][0] += tmp0*im2;
            x[j][1] += tmp1*im2;
            x[j][2] += tmp2*im2;
        } /* 16 ncons flops */
    }
    else
    {
        for (b = 0; b < ncons; b++)
        {
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0;
            x[i][1] -= tmp1;
            x[i][2] -= tmp2;
            x[j][0] += tmp0;
            x[j][1] += tmp1;
            x[j][2] += tmp2;
        }
    }
}

static void lincs_update_atoms_ind(int ncons, const int *ind, const int *bla,
                                   real prefac,
                                   const real *fac, rvec *r,
                                   const real *invmass,
                                   rvec *x)
{
    int  bi, b, i, j;
    real mvb, im1, im2, tmp0, tmp1, tmp2;

    if (invmass != NULL)
    {
        for (bi = 0; bi < ncons; bi++)
        {
            b        = ind[bi];
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            im1      = invmass[i];
            im2      = invmass[j];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0*im1;
            x[i][1] -= tmp1*im1;
            x[i][2] -= tmp2*im1;
            x[j][0] += tmp0*im2;
            x[j][1] += tmp1*im2;
            x[j][2] += tmp2*im2;
        } /* 16 ncons flops */
    }
    else
    {
        for (bi = 0; bi < ncons; bi++)
        {
            b        = ind[bi];
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;
            x[i][0] -= tmp0;
            x[i][1] -= tmp1;
            x[i][2] -= tmp2;
            x[j][0] += tmp0;
            x[j][1] += tmp1;
            x[j][2] += tmp2;
        } /* 16 ncons flops */
    }
}

static void lincs_update_atoms(struct gmx_lincsdata *li, int th,
                               real prefac,
                               const real *fac, rvec *r,
                               const real *invmass,
                               rvec *x)
{
    if (li->ntask == 1)
    {
        /* Single thread, we simply update for all constraints */
        lincs_update_atoms_noind(li->nc_real,
                                 li->bla, prefac, fac, r, invmass, x);
    }
    else
    {
        /* Update the atom vector components for our thread local
         * constraints that only access our local atom range.
         * This can be done without a barrier.
         */
        lincs_update_atoms_ind(li->task[th].nind, li->task[th].ind,
                               li->bla, prefac, fac, r, invmass, x);

        if (li->task[li->ntask].nind > 0)
        {
            /* Update the constraints that operate on atoms
             * in multiple thread atom blocks on the master thread.
             */
#pragma omp barrier
#pragma omp master
            {
                lincs_update_atoms_ind(li->task[li->ntask].nind,
                                       li->task[li->ntask].ind,
                                       li->bla, prefac, fac, r, invmass, x);
            }
        }
    }
}

#ifdef LINCS_SIMD
/* Calculate the constraint distance vectors r to project on from x.
 * Determine the right-hand side of the matrix equation using quantity f.
 * This function only differs from calc_dr_x_xp_simd below in that
 * no constraint length is subtracted and no PBC is used for f.
 */
static void gmx_simdcall
calc_dr_x_f_simd(int                       b0,
                 int                       b1,
                 const int *               bla,
                 const rvec * gmx_restrict x,
                 const rvec * gmx_restrict f,
                 const real * gmx_restrict blc,
                 const pbc_simd_t *        pbc_simd,
                 real * gmx_restrict       vbuf1,
                 real * gmx_restrict       vbuf2,
                 rvec * gmx_restrict       r,
                 real * gmx_restrict       rhs,
                 real * gmx_restrict       sol)
{
    int bs;

    assert(b0 % GMX_SIMD_REAL_WIDTH == 0);

    for (bs = b0; bs < b1; bs += GMX_SIMD_REAL_WIDTH)
    {
        gmx_simd_real_t rx_S, ry_S, rz_S, n2_S, il_S;
        gmx_simd_real_t fx_S, fy_S, fz_S, ip_S, rhs_S;

        gmx_hack_simd_gather_rvec_dist_pair_index(x, bla + bs*2, vbuf1,
                                                  &rx_S, &ry_S, &rz_S);

        pbc_correct_dx_simd(&rx_S, &ry_S, &rz_S, pbc_simd);

        n2_S  = gmx_simd_norm2_r(rx_S, ry_S, rz_S);
        il_S  = gmx_simd_invsqrt_r(n2_S);

        rx_S  = gmx_simd_mul_r(rx_S, il_S);
        ry_S  = gmx_simd_mul_r(ry_S, il_S);
        rz_S  = gmx_simd_mul_r(rz_S, il_S);

        gmx_simd_store_vec_to_rvec(rx_S, ry_S, rz_S, vbuf1, r + bs);

        gmx_hack_simd_gather_rvec_dist_pair_index(f, bla + bs*2, vbuf2,
                                                  &fx_S, &fy_S, &fz_S);

        ip_S  = gmx_simd_iprod_r(rx_S, ry_S, rz_S,
                                 fx_S, fy_S, fz_S);

        rhs_S = gmx_simd_mul_r(gmx_simd_load_r(blc + bs), ip_S);

        gmx_simd_store_r(rhs + bs, rhs_S);
        gmx_simd_store_r(sol + bs, rhs_S);
    }
}
#endif /* LINCS_SIMD */

/* LINCS projection, works on derivatives of the coordinates */
static void do_lincsp(rvec *x, rvec *f, rvec *fp, t_pbc *pbc,
                      struct gmx_lincsdata *lincsd, int th,
                      real *invmass,
                      int econq, gmx_bool bCalcDHDL,
                      gmx_bool bCalcVir, tensor rmdf)
{
    int      b0, b1, b;
    int     *bla, *blnr, *blbnb;
    rvec    *r;
    real    *blc, *blmf, *blcc, *rhs1, *rhs2, *sol;

    b0 = lincsd->task[th].b0;
    b1 = lincsd->task[th].b1;

    bla    = lincsd->bla;
    r      = lincsd->tmpv;
    blnr   = lincsd->blnr;
    blbnb  = lincsd->blbnb;
    if (econq != econqForce)
    {
        /* Use mass-weighted parameters */
        blc  = lincsd->blc;
        blmf = lincsd->blmf;
    }
    else
    {
        /* Use non mass-weighted parameters */
        blc  = lincsd->blc1;
        blmf = lincsd->blmf1;
    }
    blcc   = lincsd->tmpncc;
    rhs1   = lincsd->tmp1;
    rhs2   = lincsd->tmp2;
    sol    = lincsd->tmp3;

#ifdef LINCS_SIMD

    /* This SIMD code does the same as the plain-C code after the #else.
     * The only difference is that we always call pbc code, as with SIMD
     * the overhead of pbc computation (when not needed) is small.
     */
    pbc_simd_t pbc_simd;

    /* Convert the pbc struct for SIMD */
    set_pbc_simd(pbc, &pbc_simd);

    /* Compute normalized x i-j vectors, store in r.
     * Compute the inner product of r and xp i-j and store in rhs1.
     */
    calc_dr_x_f_simd(b0, b1, bla, x, f, blc,
                     &pbc_simd,
                     lincsd->task[th].simd_buf,
                     lincsd->task[th].simd_buf + GMX_SIMD_REAL_WIDTH*DIM,
                     r, rhs1, sol);

#else /* LINCS_SIMD */

    /* Compute normalized i-j vectors */
    if (pbc)
    {
        for (b = b0; b < b1; b++)
        {
            rvec dx;

            pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);
        }
    }
    else
    {
        for (b = b0; b < b1; b++)
        {
            rvec dx;

            rvec_sub(x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);
        } /* 16 ncons flops */
    }

    for (b = b0; b < b1; b++)
    {
        int  i, j;
        real mvb;

        i       = bla[2*b];
        j       = bla[2*b+1];
        mvb     = blc[b]*(r[b][0]*(f[i][0] - f[j][0]) +
                          r[b][1]*(f[i][1] - f[j][1]) +
                          r[b][2]*(f[i][2] - f[j][2]));
        rhs1[b] = mvb;
        sol[b]  = mvb;
        /* 7 flops */
    }

#endif /* LINCS_SIMD */

    if (lincsd->bTaskDep)
    {
        /* We need a barrier, since the matrix construction below
         * can access entries in r of other threads.
         */
#pragma omp barrier
    }

    /* Construct the (sparse) LINCS matrix */
    for (b = b0; b < b1; b++)
    {
        int n;

        for (n = blnr[b]; n < blnr[b+1]; n++)
        {
            blcc[n] = blmf[n]*iprod(r[b], r[blbnb[n]]);
        } /* 6 nr flops */
    }
    /* Together: 23*ncons + 6*nrtot flops */

    lincs_matrix_expand(lincsd, &lincsd->task[th], blcc, rhs1, rhs2, sol);
    /* nrec*(ncons+2*nrtot) flops */

    if (econq == econqDeriv_FlexCon)
    {
        /* We only want to constraint the flexible constraints,
         * so we mask out the normal ones by setting sol to 0.
         */
        for (b = b0; b < b1; b++)
        {
            if (!(lincsd->bllen0[b] == 0 && lincsd->ddist[b] == 0))
            {
                sol[b] = 0;
            }
        }
    }

    /* We multiply sol by blc, so we can use lincs_update_atoms for OpenMP */
    for (b = b0; b < b1; b++)
    {
        sol[b] *= blc[b];
    }

    /* When constraining forces, we should not use mass weighting,
     * so we pass invmass=NULL, which results in the use of 1 for all atoms.
     */
    lincs_update_atoms(lincsd, th, 1.0, sol, r,
                       (econq != econqForce) ? invmass : NULL, fp);

    if (bCalcDHDL)
    {
        real dhdlambda;

        dhdlambda = 0;
        for (b = b0; b < b1; b++)
        {
            dhdlambda -= sol[b]*lincsd->ddist[b];
        }

        lincsd->task[th].dhdlambda = dhdlambda;
    }

    if (bCalcVir)
    {
        /* Constraint virial,
         * determines sum r_bond x delta f,
         * where delta f is the constraint correction
         * of the quantity that is being constrained.
         */
        for (b = b0; b < b1; b++)
        {
            real mvb, tmp1;
            int  i, j;

            mvb = lincsd->bllen[b]*sol[b];
            for (i = 0; i < DIM; i++)
            {
                tmp1 = mvb*r[b][i];
                for (j = 0; j < DIM; j++)
                {
                    rmdf[i][j] += tmp1*r[b][j];
                }
            }
        } /* 23 ncons flops */
    }
}

#ifdef LINCS_SIMD
/* Calculate the constraint distance vectors r to project on from x.
 * Determine the right-hand side of the matrix equation using coordinates xp.
 */
static void gmx_simdcall
calc_dr_x_xp_simd(int                       b0,
                  int                       b1,
                  const int *               bla,
                  const rvec * gmx_restrict x,
                  const rvec * gmx_restrict xp,
                  const real * gmx_restrict bllen,
                  const real * gmx_restrict blc,
                  const pbc_simd_t *        pbc_simd,
                  real * gmx_restrict       vbuf1,
                  real * gmx_restrict       vbuf2,
                  rvec * gmx_restrict       r,
                  real * gmx_restrict       rhs,
                  real * gmx_restrict       sol)
{
    int bs;

    assert(b0 % GMX_SIMD_REAL_WIDTH == 0);

    for (bs = b0; bs < b1; bs += GMX_SIMD_REAL_WIDTH)
    {
        gmx_simd_real_t rx_S, ry_S, rz_S, n2_S, il_S;
        gmx_simd_real_t rxp_S, ryp_S, rzp_S, ip_S, rhs_S;

        gmx_hack_simd_gather_rvec_dist_pair_index(x, bla + bs*2, vbuf1,
                                                  &rx_S, &ry_S, &rz_S);

        pbc_correct_dx_simd(&rx_S, &ry_S, &rz_S, pbc_simd);

        n2_S  = gmx_simd_norm2_r(rx_S, ry_S, rz_S);
        il_S  = gmx_simd_invsqrt_r(n2_S);

        rx_S  = gmx_simd_mul_r(rx_S, il_S);
        ry_S  = gmx_simd_mul_r(ry_S, il_S);
        rz_S  = gmx_simd_mul_r(rz_S, il_S);

        gmx_simd_store_vec_to_rvec(rx_S, ry_S, rz_S, vbuf1, r + bs);

        gmx_hack_simd_gather_rvec_dist_pair_index(xp, bla + bs*2, vbuf2,
                                                  &rxp_S, &ryp_S, &rzp_S);

        pbc_correct_dx_simd(&rxp_S, &ryp_S, &rzp_S, pbc_simd);

        ip_S  = gmx_simd_iprod_r(rx_S, ry_S, rz_S,
                                 rxp_S, ryp_S, rzp_S);

        rhs_S = gmx_simd_mul_r(gmx_simd_load_r(blc + bs),
                               gmx_simd_sub_r(ip_S, gmx_simd_load_r(bllen + bs)));

        gmx_simd_store_r(rhs + bs, rhs_S);
        gmx_simd_store_r(sol + bs, rhs_S);
    }
}
#endif /* LINCS_SIMD */

/* Determine the distances and right-hand side for the next iteration */
static void calc_dist_iter(int                       b0,
                           int                       b1,
                           const int *               bla,
                           const rvec * gmx_restrict xp,
                           const real * gmx_restrict bllen,
                           const real * gmx_restrict blc,
                           const t_pbc *             pbc,
                           real                      wfac,
                           real * gmx_restrict       rhs,
                           real * gmx_restrict       sol,
                           gmx_bool *                bWarn)
{
    int b;

    for (b = b0; b < b1; b++)
    {
        real len, len2, dlen2, mvb;
        rvec dx;

        len = bllen[b];
        if (pbc)
        {
            pbc_dx_aiuc(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx);
        }
        else
        {
            rvec_sub(xp[bla[2*b]], xp[bla[2*b+1]], dx);
        }
        len2  = len*len;
        dlen2 = 2*len2 - norm2(dx);
        if (dlen2 < wfac*len2)
        {
            /* not race free - see detailed comment in caller */
            *bWarn = TRUE;
        }
        if (dlen2 > 0)
        {
            mvb = blc[b]*(len - dlen2*gmx_invsqrt(dlen2));
        }
        else
        {
            mvb = blc[b]*len;
        }
        rhs[b]  = mvb;
        sol[b]  = mvb;
    } /* 20*ncons flops */
}

#ifdef LINCS_SIMD
/* As the function above, but using SIMD intrinsics */
static void gmx_simdcall
calc_dist_iter_simd(int                       b0,
                    int                       b1,
                    const int *               bla,
                    const rvec * gmx_restrict x,
                    const real * gmx_restrict bllen,
                    const real * gmx_restrict blc,
                    const pbc_simd_t *        pbc_simd,
                    real                      wfac,
                    real * gmx_restrict       vbuf,
                    real * gmx_restrict       rhs,
                    real * gmx_restrict       sol,
                    gmx_bool *                bWarn)
{
    gmx_simd_real_t min_S  = gmx_simd_set1_r(GMX_REAL_MIN);
    gmx_simd_real_t two_S  = gmx_simd_set1_r(2.0);
    gmx_simd_real_t wfac_S = gmx_simd_set1_r(wfac);
    gmx_simd_bool_t warn_S;

    int             bs;

    assert(b0 % GMX_SIMD_REAL_WIDTH == 0);

    /* Initialize all to FALSE */
    warn_S = gmx_simd_cmplt_r(two_S, gmx_simd_setzero_r());

    for (bs = b0; bs < b1; bs += GMX_SIMD_REAL_WIDTH)
    {
        gmx_simd_real_t rx_S, ry_S, rz_S, n2_S;
        gmx_simd_real_t len_S, len2_S, dlen2_S, lc_S, blc_S;

        gmx_hack_simd_gather_rvec_dist_pair_index(x, bla + bs*2, vbuf,
                                                  &rx_S, &ry_S, &rz_S);

        pbc_correct_dx_simd(&rx_S, &ry_S, &rz_S, pbc_simd);

        n2_S    = gmx_simd_norm2_r(rx_S, ry_S, rz_S);

        len_S   = gmx_simd_load_r(bllen + bs);
        len2_S  = gmx_simd_mul_r(len_S, len_S);

        dlen2_S = gmx_simd_fmsub_r(two_S, len2_S, n2_S);

        warn_S  = gmx_simd_or_b(warn_S,
                                gmx_simd_cmplt_r(dlen2_S,
                                                 gmx_simd_mul_r(wfac_S, len2_S)));

        /* Avoid 1/0 by taking the max with REAL_MIN.
         * Note: when dlen2 is close to zero (90 degree constraint rotation),
         * the accuracy of the algorithm is no longer relevant.
         */
        dlen2_S = gmx_simd_max_r(dlen2_S, min_S);

        lc_S    = gmx_simd_fnmadd_r(dlen2_S, gmx_simd_invsqrt_r(dlen2_S), len_S);

        blc_S   = gmx_simd_load_r(blc + bs);

        lc_S    = gmx_simd_mul_r(blc_S, lc_S);

        gmx_simd_store_r(rhs + bs, lc_S);
        gmx_simd_store_r(sol + bs, lc_S);
    }

    if (gmx_simd_anytrue_b(warn_S))
    {
        *bWarn = TRUE;
    }
}
#endif /* LINCS_SIMD */

static void do_lincs(rvec *x, rvec *xp, matrix box, t_pbc *pbc,
                     struct gmx_lincsdata *lincsd, int th,
                     const real *invmass,
                     t_commrec *cr,
                     gmx_bool bCalcDHDL,
                     real wangle, gmx_bool *bWarn,
                     real invdt, rvec * gmx_restrict v,
                     gmx_bool bCalcVir, tensor vir_r_m_dr)
{
    int      b0, b1, b, i, j, n, iter;
    int     *bla, *blnr, *blbnb;
    rvec    *r;
    real    *blc, *blmf, *bllen, *blcc, *rhs1, *rhs2, *sol, *blc_sol, *mlambda;
    int     *nlocat;

    b0 = lincsd->task[th].b0;
    b1 = lincsd->task[th].b1;

    bla     = lincsd->bla;
    r       = lincsd->tmpv;
    blnr    = lincsd->blnr;
    blbnb   = lincsd->blbnb;
    blc     = lincsd->blc;
    blmf    = lincsd->blmf;
    bllen   = lincsd->bllen;
    blcc    = lincsd->tmpncc;
    rhs1    = lincsd->tmp1;
    rhs2    = lincsd->tmp2;
    sol     = lincsd->tmp3;
    blc_sol = lincsd->tmp4;
    mlambda = lincsd->mlambda;
    nlocat  = lincsd->nlocat;

#ifdef LINCS_SIMD

    /* This SIMD code does the same as the plain-C code after the #else.
     * The only difference is that we always call pbc code, as with SIMD
     * the overhead of pbc computation (when not needed) is small.
     */
    pbc_simd_t pbc_simd;

    /* Convert the pbc struct for SIMD */
    set_pbc_simd(pbc, &pbc_simd);

    /* Compute normalized x i-j vectors, store in r.
     * Compute the inner product of r and xp i-j and store in rhs1.
     */
    calc_dr_x_xp_simd(b0, b1, bla, x, xp, bllen, blc,
                      &pbc_simd,
                      lincsd->task[th].simd_buf,
                      lincsd->task[th].simd_buf + GMX_SIMD_REAL_WIDTH*DIM,
                      r, rhs1, sol);

#else /* LINCS_SIMD */

    if (pbc)
    {
        /* Compute normalized i-j vectors */
        for (b = b0; b < b1; b++)
        {
            rvec dx;
            real mvb;

            pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
            unitv(dx, r[b]);

            pbc_dx_aiuc(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx);
            mvb     = blc[b]*(iprod(r[b], dx) - bllen[b]);
            rhs1[b] = mvb;
            sol[b]  = mvb;
        }
    }
    else
    {
        /* Compute normalized i-j vectors */
        for (b = b0; b < b1; b++)
        {
            real tmp0, tmp1, tmp2, rlen, mvb;

            i       = bla[2*b];
            j       = bla[2*b+1];
            tmp0    = x[i][0] - x[j][0];
            tmp1    = x[i][1] - x[j][1];
            tmp2    = x[i][2] - x[j][2];
            rlen    = gmx_invsqrt(tmp0*tmp0 + tmp1*tmp1 + tmp2*tmp2);
            r[b][0] = rlen*tmp0;
            r[b][1] = rlen*tmp1;
            r[b][2] = rlen*tmp2;
            /* 16 ncons flops */

            i       = bla[2*b];
            j       = bla[2*b+1];
            mvb     = blc[b]*(r[b][0]*(xp[i][0] - xp[j][0]) +
                              r[b][1]*(xp[i][1] - xp[j][1]) +
                              r[b][2]*(xp[i][2] - xp[j][2]) - bllen[b]);
            rhs1[b] = mvb;
            sol[b]  = mvb;
            /* 10 flops */
        }
        /* Together: 26*ncons + 6*nrtot flops */
    }

#endif /* LINCS_SIMD */

    if (lincsd->bTaskDep)
    {
        /* We need a barrier, since the matrix construction below
         * can access entries in r of other threads.
         */
#pragma omp barrier
    }

    /* Construct the (sparse) LINCS matrix */
    for (b = b0; b < b1; b++)
    {
        for (n = blnr[b]; n < blnr[b+1]; n++)
        {
            blcc[n] = blmf[n]*iprod(r[b], r[blbnb[n]]);
        }
    }
    /* Together: 26*ncons + 6*nrtot flops */

    lincs_matrix_expand(lincsd, &lincsd->task[th], blcc, rhs1, rhs2, sol);
    /* nrec*(ncons+2*nrtot) flops */

#ifndef LINCS_SIMD
    for (b = b0; b < b1; b++)
    {
        mlambda[b] = blc[b]*sol[b];
    }
#else
    for (b = b0; b < b1; b += GMX_SIMD_REAL_WIDTH)
    {
        gmx_simd_store_r(mlambda + b,
                         gmx_simd_mul_r(gmx_simd_load_r(blc + b),
                                        gmx_simd_load_r(sol + b)));
    }
#endif

    /* Update the coordinates */
    lincs_update_atoms(lincsd, th, 1.0, mlambda, r, invmass, xp);

    /*
     ********  Correction for centripetal effects  ********
     */

    real wfac;

    wfac = cos(DEG2RAD*wangle);
    wfac = wfac*wfac;

    for (iter = 0; iter < lincsd->nIter; iter++)
    {
        if ((lincsd->bCommIter && DOMAINDECOMP(cr) && cr->dd->constraints))
        {
#pragma omp barrier
#pragma omp master
            {
                /* Communicate the corrected non-local coordinates */
                if (DOMAINDECOMP(cr))
                {
                    dd_move_x_constraints(cr->dd, box, xp, NULL, FALSE);
                }
            }
#pragma omp barrier
        }
        else if (lincsd->bTaskDep)
        {
#pragma omp barrier
        }

#ifdef LINCS_SIMD
        calc_dist_iter_simd(b0, b1, bla, xp, bllen, blc, &pbc_simd, wfac,
                            lincsd->task[th].simd_buf, rhs1, sol, bWarn);
#else
        calc_dist_iter(b0, b1, bla, xp, bllen, blc, pbc, wfac,
                       rhs1, sol, bWarn);
        /* 20*ncons flops */
#endif

        lincs_matrix_expand(lincsd, &lincsd->task[th], blcc, rhs1, rhs2, sol);
        /* nrec*(ncons+2*nrtot) flops */

#ifndef LINCS_SIMD
        for (b = b0; b < b1; b++)
        {
            real mvb;

            mvb         = blc[b]*sol[b];
            blc_sol[b]  = mvb;
            mlambda[b] += mvb;
        }
#else
        for (b = b0; b < b1; b += GMX_SIMD_REAL_WIDTH)
        {
            gmx_simd_real_t mvb;

            mvb = gmx_simd_mul_r(gmx_simd_load_r(blc + b),
                                 gmx_simd_load_r(sol + b));
            gmx_simd_store_r(blc_sol + b, mvb);
            gmx_simd_store_r(mlambda + b,
                             gmx_simd_add_r(gmx_simd_load_r(mlambda + b), mvb));
        }
#endif

        /* Update the coordinates */
        lincs_update_atoms(lincsd, th, 1.0, blc_sol, r, invmass, xp);
    }
    /* nit*ncons*(37+9*nrec) flops */

    if (v != NULL)
    {
        /* Update the velocities */
        lincs_update_atoms(lincsd, th, invdt, mlambda, r, invmass, v);
        /* 16 ncons flops */
    }

    if (nlocat != NULL && (bCalcDHDL || bCalcVir))
    {
        if (lincsd->bTaskDep)
        {
            /* In lincs_update_atoms threads might cross-read mlambda */
#pragma omp barrier
        }

        /* Only account for local atoms */
        for (b = b0; b < b1; b++)
        {
            mlambda[b] *= 0.5*nlocat[b];
        }
    }

    if (bCalcDHDL)
    {
        real dhdl;

        dhdl = 0;
        for (b = b0; b < b1; b++)
        {
            /* Note that this this is dhdl*dt^2, the dt^2 factor is corrected
             * later after the contributions are reduced over the threads.
             */
            dhdl -= lincsd->mlambda[b]*lincsd->ddist[b];
        }
        lincsd->task[th].dhdlambda = dhdl;
    }

    if (bCalcVir)
    {
        /* Constraint virial */
        for (b = b0; b < b1; b++)
        {
            real tmp0, tmp1;

            tmp0 = -bllen[b]*mlambda[b];
            for (i = 0; i < DIM; i++)
            {
                tmp1 = tmp0*r[b][i];
                for (j = 0; j < DIM; j++)
                {
                    vir_r_m_dr[i][j] -= tmp1*r[b][j];
                }
            }
        } /* 22 ncons flops */
    }

    /* Total:
     * 26*ncons + 6*nrtot + nrec*(ncons+2*nrtot)
     * + nit * (20*ncons + nrec*(ncons+2*nrtot) + 17 ncons)
     *
     * (26+nrec)*ncons + (6+2*nrec)*nrtot
     * + nit * ((37+nrec)*ncons + 2*nrec*nrtot)
     * if nit=1
     * (63+nrec)*ncons + (6+4*nrec)*nrtot
     */
}

/* Sets the elements in the LINCS matrix for task li_task */
static void set_lincs_matrix_task(struct gmx_lincsdata *li,
                                  lincs_task_t         *li_task,
                                  const real           *invmass,
                                  int                  *ncc_triangle)
{
    int        i;

    /* Construct the coupling coefficient matrix blmf */
    li_task->ntriangle = 0;
    *ncc_triangle      = 0;
    for (i = li_task->b0; i < li_task->b1; i++)
    {
        int a1, a2, n;

        a1 = li->bla[2*i];
        a2 = li->bla[2*i+1];
        for (n = li->blnr[i]; (n < li->blnr[i+1]); n++)
        {
            int k, sign, center, end;

            k = li->blbnb[n];

            /* If we are using multiple, independent tasks for LINCS,
             * the calls to check_assign_connected should have
             * put all connected constraints in our task.
             */
            assert(li->bTaskDep || (k >= li_task->b0 && k < li_task->b1));

            if (a1 == li->bla[2*k] || a2 == li->bla[2*k+1])
            {
                sign = -1;
            }
            else
            {
                sign = 1;
            }
            if (a1 == li->bla[2*k] || a1 == li->bla[2*k+1])
            {
                center = a1;
                end    = a2;
            }
            else
            {
                center = a2;
                end    = a1;
            }
            li->blmf[n]  = sign*invmass[center]*li->blc[i]*li->blc[k];
            li->blmf1[n] = sign*0.5;
            if (li->ncg_triangle > 0)
            {
                int nk, kk;

                /* Look for constraint triangles */
                for (nk = li->blnr[k]; (nk < li->blnr[k+1]); nk++)
                {
                    kk = li->blbnb[nk];
                    if (kk != i && kk != k &&
                        (li->bla[2*kk] == end || li->bla[2*kk+1] == end))
                    {
                        /* If we are using multiple tasks for LINCS,
                         * the calls to check_assign_triangle should have
                         * put all constraints in the triangle in our task.
                         */
                        assert(k  >= li_task->b0 && k  < li_task->b1);
                        assert(kk >= li_task->b0 && kk < li_task->b1);

                        if (li_task->ntriangle == 0 ||
                            li_task->triangle[li_task->ntriangle - 1] < i)
                        {
                            /* Add this constraint to the triangle list */
                            li_task->triangle[li_task->ntriangle] = i;
                            li_task->tri_bits[li_task->ntriangle] = 0;
                            li_task->ntriangle++;
                            if (li->blnr[i+1] - li->blnr[i] > static_cast<int>(sizeof(li_task->tri_bits[0])*8 - 1))
                            {
                                gmx_fatal(FARGS, "A constraint is connected to %d constraints, this is more than the %d allowed for constraints participating in triangles",
                                          li->blnr[i+1] - li->blnr[i],
                                          sizeof(li_task->tri_bits[0])*8-1);
                            }
                        }
                        li_task->tri_bits[li_task->ntriangle-1] |= (1 << (n - li->blnr[i]));
                        (*ncc_triangle)++;
                    }
                }
            }
        }
    }
}

/* Sets the elements in the LINCS matrix */
void set_lincs_matrix(struct gmx_lincsdata *li, real *invmass, real lambda)
{
    int        i;
    const real invsqrt2 = 0.7071067811865475244;

    for (i = 0; (i < li->nc); i++)
    {
        int a1, a2;

        a1          = li->bla[2*i];
        a2          = li->bla[2*i+1];
        li->blc[i]  = gmx_invsqrt(invmass[a1] + invmass[a2]);
        li->blc1[i] = invsqrt2;
    }

    /* Construct the coupling coefficient matrix blmf */
    int th, ntriangle = 0, ncc_triangle = 0;
#pragma omp parallel for reduction(+: ntriangle, ncc_triangle) num_threads(li->ntask) schedule(static)
    for (th = 0; th < li->ntask; th++)
    {
        set_lincs_matrix_task(li, &li->task[th], invmass, &ncc_triangle);
        ntriangle = li->task[th].ntriangle;
    }
    li->ntriangle    = ntriangle;
    li->ncc_triangle = ncc_triangle;

    if (debug)
    {
        fprintf(debug, "Of the %d constraints %d participate in triangles\n",
                li->nc, li->ntriangle);
        fprintf(debug, "There are %d couplings of which %d in triangles\n",
                li->ncc, li->ncc_triangle);
    }

    /* Set matlam,
     * so we know with which lambda value the masses have been set.
     */
    li->matlam = lambda;
}

static int count_triangle_constraints(t_ilist *ilist, t_blocka *at2con)
{
    int      ncon1, ncon_tot;
    int      c0, a00, a01, n1, c1, a10, a11, ac1, n2, c2, a20, a21;
    int      ncon_triangle;
    gmx_bool bTriangle;
    t_iatom *ia1, *ia2, *iap;

    ncon1    = ilist[F_CONSTR].nr/3;
    ncon_tot = ncon1 + ilist[F_CONSTRNC].nr/3;

    ia1 = ilist[F_CONSTR].iatoms;
    ia2 = ilist[F_CONSTRNC].iatoms;

    ncon_triangle = 0;
    for (c0 = 0; c0 < ncon_tot; c0++)
    {
        bTriangle = FALSE;
        iap       = constr_iatomptr(ncon1, ia1, ia2, c0);
        a00       = iap[1];
        a01       = iap[2];
        for (n1 = at2con->index[a01]; n1 < at2con->index[a01+1]; n1++)
        {
            c1 = at2con->a[n1];
            if (c1 != c0)
            {
                iap = constr_iatomptr(ncon1, ia1, ia2, c1);
                a10 = iap[1];
                a11 = iap[2];
                if (a10 == a01)
                {
                    ac1 = a11;
                }
                else
                {
                    ac1 = a10;
                }
                for (n2 = at2con->index[ac1]; n2 < at2con->index[ac1+1]; n2++)
                {
                    c2 = at2con->a[n2];
                    if (c2 != c0 && c2 != c1)
                    {
                        iap = constr_iatomptr(ncon1, ia1, ia2, c2);
                        a20 = iap[1];
                        a21 = iap[2];
                        if (a20 == a00 || a21 == a00)
                        {
                            bTriangle = TRUE;
                        }
                    }
                }
            }
        }
        if (bTriangle)
        {
            ncon_triangle++;
        }
    }

    return ncon_triangle;
}

static gmx_bool more_than_two_sequential_constraints(const t_ilist  *ilist,
                                                     const t_blocka *at2con)
{
    t_iatom  *ia1, *ia2, *iap;
    int       ncon1, ncon_tot, c;
    int       a1, a2;
    gmx_bool  bMoreThanTwoSequentialConstraints;

    ncon1    = ilist[F_CONSTR].nr/3;
    ncon_tot = ncon1 + ilist[F_CONSTRNC].nr/3;

    ia1 = ilist[F_CONSTR].iatoms;
    ia2 = ilist[F_CONSTRNC].iatoms;

    bMoreThanTwoSequentialConstraints = FALSE;
    for (c = 0; c < ncon_tot && !bMoreThanTwoSequentialConstraints; c++)
    {
        iap = constr_iatomptr(ncon1, ia1, ia2, c);
        a1  = iap[1];
        a2  = iap[2];
        /* Check if this constraint has constraints connected at both atoms */
        if (at2con->index[a1+1] - at2con->index[a1] > 1 &&
            at2con->index[a2+1] - at2con->index[a2] > 1)
        {
            bMoreThanTwoSequentialConstraints = TRUE;
        }
    }

    return bMoreThanTwoSequentialConstraints;
}

static int int_comp(const void *a, const void *b)
{
    return (*(int *)a) - (*(int *)b);
}

gmx_lincsdata_t init_lincs(FILE *fplog, gmx_mtop_t *mtop,
                           int nflexcon_global, t_blocka *at2con,
                           gmx_bool bPLINCS, int nIter, int nProjOrder)
{
    struct gmx_lincsdata *li;
    int                   mt, mb;
    gmx_moltype_t        *molt;
    gmx_bool              bMoreThanTwoSeq;

    if (fplog)
    {
        fprintf(fplog, "\nInitializing%s LINear Constraint Solver\n",
                bPLINCS ? " Parallel" : "");
    }

    snew(li, 1);

    li->ncg      =
        gmx_mtop_ftype_count(mtop, F_CONSTR) +
        gmx_mtop_ftype_count(mtop, F_CONSTRNC);
    li->ncg_flex = nflexcon_global;

    li->nIter  = nIter;
    li->nOrder = nProjOrder;

    li->max_connect = 0;
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        int a;

        molt = &mtop->moltype[mt];
        for (a = 0; a < molt->atoms.nr; a++)
        {
            li->max_connect = std::max(li->max_connect,
                                       at2con[mt].index[a + 1] - at2con[mt].index[a]);
        }
    }

    li->ncg_triangle = 0;
    bMoreThanTwoSeq  = FALSE;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt              = &mtop->moltype[mtop->molblock[mb].type];

        li->ncg_triangle +=
            mtop->molblock[mb].nmol*
            count_triangle_constraints(molt->ilist,
                                       &at2con[mtop->molblock[mb].type]);

        if (!bMoreThanTwoSeq &&
            more_than_two_sequential_constraints(molt->ilist, &at2con[mtop->molblock[mb].type]))
        {
            bMoreThanTwoSeq = TRUE;
        }
    }

    /* Check if we need to communicate not only before LINCS,
     * but also before each iteration.
     * The check for only two sequential constraints is only
     * useful for the common case of H-bond only constraints.
     * With more effort we could also make it useful for small
     * molecules with nr. sequential constraints <= nOrder-1.
     */
    li->bCommIter = (bPLINCS && (li->nOrder < 1 || bMoreThanTwoSeq));

    if (debug && bPLINCS)
    {
        fprintf(debug, "PLINCS communication before each iteration: %d\n",
                li->bCommIter);
    }

    /* LINCS can run on any number of threads.
     * Currently the number is fixed for the whole simulation,
     * but it could be set in set_lincs().
     * The current constraint to task assignment code can create independent
     * tasks only when not more than two constraints are connected sequentially.
     */
    li->ntask    = gmx_omp_nthreads_get(emntLINCS);
    li->bTaskDep = (li->ntask > 1 && bMoreThanTwoSeq);
    if (debug)
    {
        fprintf(debug, "LINCS: using %d threads, tasks are %sdependent\n",
                li->ntask, li->bTaskDep ? "" : "in");
    }
    if (li->ntask == 1)
    {
        snew(li->task, 1);
    }
    else
    {
        /* Allocate an extra elements for "task-overlap" constraints */
        snew(li->task, li->ntask + 1);
    }
    int th;
#pragma omp parallel for num_threads(li->ntask)
    for (th = 0; th < li->ntask; th++)
    {
        /* Per thread SIMD load buffer for loading 2 simd_width rvecs */
        snew_aligned(li->task[th].simd_buf, 2*simd_width*DIM,
                     align_bytes);
    }

    if (bPLINCS || li->ncg_triangle > 0)
    {
        please_cite(fplog, "Hess2008a");
    }
    else
    {
        please_cite(fplog, "Hess97a");
    }

    if (fplog)
    {
        fprintf(fplog, "The number of constraints is %d\n", li->ncg);
        if (bPLINCS)
        {
            fprintf(fplog, "There are inter charge-group constraints,\n"
                    "will communicate selected coordinates each lincs iteration\n");
        }
        if (li->ncg_triangle > 0)
        {
            fprintf(fplog,
                    "%d constraints are involved in constraint triangles,\n"
                    "will apply an additional matrix expansion of order %d for couplings\n"
                    "between constraints inside triangles\n",
                    li->ncg_triangle, li->nOrder);
        }
    }

    return li;
}

/* Sets up the work division over the threads */
static void lincs_thread_setup(struct gmx_lincsdata *li, int natoms)
{
    lincs_task_t   *li_m;
    int             th;
    gmx_bitmask_t  *atf;
    int             a;

    if (natoms > li->atf_nalloc)
    {
        li->atf_nalloc = over_alloc_large(natoms);
        srenew(li->atf, li->atf_nalloc);
    }

    atf = li->atf;
    /* Clear the atom flags */
    for (a = 0; a < natoms; a++)
    {
        bitmask_clear(&atf[a]);
    }

    if (li->ntask > BITMASK_SIZE)
    {
        gmx_fatal(FARGS, "More than %d threads is not supported for LINCS.", BITMASK_SIZE);
    }

    for (th = 0; th < li->ntask; th++)
    {
        lincs_task_t *li_task;
        int           b;

        li_task = &li->task[th];

        /* For each atom set a flag for constraints from each */
        for (b = li_task->b0; b < li_task->b1; b++)
        {
            bitmask_set_bit(&atf[li->bla[b*2    ]], th);
            bitmask_set_bit(&atf[li->bla[b*2 + 1]], th);
        }
    }

#pragma omp parallel for num_threads(li->ntask) schedule(static)
    for (th = 0; th < li->ntask; th++)
    {
        lincs_task_t  *li_task;
        gmx_bitmask_t  mask;
        int            b;

        li_task = &li->task[th];

        if (li_task->b1 - li_task->b0 > li_task->ind_nalloc)
        {
            li_task->ind_nalloc = over_alloc_large(li_task->b1-li_task->b0);
            srenew(li_task->ind, li_task->ind_nalloc);
            srenew(li_task->ind_r, li_task->ind_nalloc);
        }

        bitmask_init_low_bits(&mask, th);

        li_task->nind   = 0;
        li_task->nind_r = 0;
        for (b = li_task->b0; b < li_task->b1; b++)
        {
            /* We let the constraint with the lowest thread index
             * operate on atoms with constraints from multiple threads.
             */
            if (bitmask_is_disjoint(atf[li->bla[b*2]], mask) &&
                bitmask_is_disjoint(atf[li->bla[b*2+1]], mask))
            {
                /* Add the constraint to the local atom update index */
                li_task->ind[li_task->nind++] = b;
            }
            else
            {
                /* Add the constraint to the rest block */
                li_task->ind_r[li_task->nind_r++] = b;
            }
        }
    }

    /* We need to copy all constraints which have not be assigned
     * to a thread to a separate list which will be handled by one thread.
     */
    li_m = &li->task[li->ntask];

    li_m->nind = 0;
    for (th = 0; th < li->ntask; th++)
    {
        lincs_task_t *li_task;
        int           b;

        li_task   = &li->task[th];

        if (li_m->nind + li_task->nind_r > li_m->ind_nalloc)
        {
            li_m->ind_nalloc = over_alloc_large(li_m->nind+li_task->nind_r);
            srenew(li_m->ind, li_m->ind_nalloc);
        }

        for (b = 0; b < li_task->nind_r; b++)
        {
            li_m->ind[li_m->nind++] = li_task->ind_r[b];
        }

        if (debug)
        {
            fprintf(debug, "LINCS thread %d: %d constraints\n",
                    th, li_task->nind);
        }
    }

    if (debug)
    {
        fprintf(debug, "LINCS thread r: %d constraints\n",
                li_m->nind);
    }
}

/* There is no realloc with alignment, so here we make one for reals.
 * Note that this function does not preserve the contents of the memory.
 */
static void resize_real_aligned(real **ptr, int nelem)
{
    sfree_aligned(*ptr);
    snew_aligned(*ptr, nelem, align_bytes);
}

static void assign_constraint(struct gmx_lincsdata *li,
                              int constraint_index,
                              int a1, int a2,
                              real lenA, real lenB,
                              const t_blocka *at2con)
{
    int con;

    con = li->nc;

    /* Make an mapping of local topology constraint index to LINCS index */
    li->con_index[constraint_index] = con;

    li->bllen0[con]  = lenA;
    li->ddist[con]   = lenB - lenA;
    /* Set the length to the topology A length */
    li->bllen[con]   = lenA;
    li->bla[2*con]   = a1;
    li->bla[2*con+1] = a2;

    /* Make space in the constraint connection matrix for constraints
     * connected to both end of the current constraint.
     */
    li->ncc +=
        at2con->index[a1 + 1] - at2con->index[a1] - 1 +
        at2con->index[a2 + 1] - at2con->index[a2] - 1;

    li->blnr[con + 1] = li->ncc;

    /* Increase the constraint count */
    li->nc++;
}

/* Check if constraint with topology index constraint_index is connected
 * to other constraints, and if so add those connected constraints to our task.
 */
static void check_assign_connected(struct gmx_lincsdata *li,
                                   const t_iatom *iatom,
                                   const t_idef *idef,
                                   int bDynamics,
                                   int a1, int a2,
                                   const t_blocka *at2con)
{
    /* Currently this function only supports constraint groups
     * in which all constraints share at least one atom
     * (e.g. H-bond constraints).
     * Check both ends of the current constraint for
     * connected constraints. We need to assign those
     * to the same task.
     */
    int end;

    for (end = 0; end < 2; end++)
    {
        int a, k;

        a = (end == 0 ? a1 : a2);

        for (k = at2con->index[a]; k < at2con->index[a + 1]; k++)
        {
            int cc;

            cc = at2con->a[k];
            /* Check if constraint cc has not yet been assigned */
            if (li->con_index[cc] == -1)
            {
                int  type;
                real lenA, lenB;

                type = iatom[cc*3];
                lenA = idef->iparams[type].constr.dA;
                lenB = idef->iparams[type].constr.dB;

                if (bDynamics || lenA != 0 || lenB != 0)
                {
                    assign_constraint(li, cc, iatom[3*cc + 1], iatom[3*cc + 2], lenA, lenB, at2con);
                }
            }
        }
    }
}

/* Check if constraint with topology index constraint_index is involved
 * in a constraint triangle, and if so add the other two constraints
 * in the triangle to our task.
 */
static void check_assign_triangle(struct gmx_lincsdata *li,
                                  const t_iatom *iatom,
                                  const t_idef *idef,
                                  int bDynamics,
                                  int constraint_index,
                                  int a1, int a2,
                                  const t_blocka *at2con)
{
    int nca, cc[32], ca[32], k;
    int c_triangle[2] = { -1, -1 };

    nca = 0;
    for (k = at2con->index[a1]; k < at2con->index[a1 + 1]; k++)
    {
        int c;

        c = at2con->a[k];
        if (c != constraint_index)
        {
            int aa1, aa2;

            aa1 = iatom[c*3 + 1];
            aa2 = iatom[c*3 + 2];
            if (aa1 != a1)
            {
                cc[nca] = c;
                ca[nca] = aa1;
                nca++;
            }
            if (aa2 != a1)
            {
                cc[nca] = c;
                ca[nca] = aa2;
                nca++;
            }
        }
    }

    for (k = at2con->index[a2]; k < at2con->index[a2 + 1]; k++)
    {
        int c;

        c = at2con->a[k];
        if (c != constraint_index)
        {
            int aa1, aa2, i;

            aa1 = iatom[c*3 + 1];
            aa2 = iatom[c*3 + 2];
            if (aa1 != a2)
            {
                for (i = 0; i < nca; i++)
                {
                    if (aa1 == ca[i])
                    {
                        c_triangle[0] = cc[i];
                        c_triangle[1] = c;
                    }
                }
            }
            if (aa2 != a2)
            {
                for (i = 0; i < nca; i++)
                {
                    if (aa2 == ca[i])
                    {
                        c_triangle[0] = cc[i];
                        c_triangle[1] = c;
                    }
                }
            }
        }
    }

    if (c_triangle[0] >= 0)
    {
        int end;

        for (end = 0; end < 2; end++)
        {
            /* Check if constraint c_triangle[end] has not yet been assigned */
            if (li->con_index[c_triangle[end]] == -1)
            {
                int  i, type;
                real lenA, lenB;

                i    = c_triangle[end]*3;
                type = iatom[i];
                lenA = idef->iparams[type].constr.dA;
                lenB = idef->iparams[type].constr.dB;

                if (bDynamics || lenA != 0 || lenB != 0)
                {
                    assign_constraint(li, c_triangle[end], iatom[i + 1], iatom[i + 2], lenA, lenB, at2con);
                }
            }
        }
    }
}

static void set_matrix_indices(struct gmx_lincsdata *li,
                               const lincs_task_t   *li_task,
                               const t_blocka       *at2con,
                               gmx_bool              bSortMatrix)
{
    int b;

    for (b = li_task->b0; b < li_task->b1; b++)
    {
        int a1, a2, i, k;

        a1 = li->bla[b*2];
        a2 = li->bla[b*2 + 1];

        i = li->blnr[b];
        for (k = at2con->index[a1]; k < at2con->index[a1 + 1]; k++)
        {
            int concon;

            concon = li->con_index[at2con->a[k]];
            if (concon != b)
            {
                li->blbnb[i++] = concon;
            }
        }
        for (k = at2con->index[a2]; k < at2con->index[a2 + 1]; k++)
        {
            int concon;

            concon = li->con_index[at2con->a[k]];
            if (concon != b)
            {
                li->blbnb[i++] = concon;
            }
        }

        if (bSortMatrix)
        {
            /* Order the blbnb matrix to optimize memory access */
            qsort(&(li->blbnb[li->blnr[b]]), li->blnr[b + 1] - li->blnr[b],
                  sizeof(li->blbnb[0]), int_comp);
        }
    }
}

void set_lincs(t_idef *idef, t_mdatoms *md,
               gmx_bool bDynamics, t_commrec *cr,
               struct gmx_lincsdata *li)
{
    int          natoms, nflexcon;
    t_blocka     at2con;
    t_iatom     *iatom;
    int          i, ncc_alloc_old, ncon_tot;

    li->nc_real = 0;
    li->nc      = 0;
    li->ncc     = 0;
    /* Zero the thread index ranges.
     * Otherwise without local constraints we could return with old ranges.
     */
    for (i = 0; i < li->ntask; i++)
    {
        li->task[i].b0   = 0;
        li->task[i].b1   = 0;
        li->task[i].nind = 0;
    }
    if (li->ntask > 1)
    {
        li->task[li->ntask].nind = 0;
    }

    /* This is the local topology, so there are only F_CONSTR constraints */
    if (idef->il[F_CONSTR].nr == 0)
    {
        /* There are no constraints,
         * we do not need to fill any data structures.
         */
        return;
    }

    if (debug)
    {
        fprintf(debug, "Building the LINCS connectivity\n");
    }

    if (DOMAINDECOMP(cr))
    {
        if (cr->dd->constraints)
        {
            int start;

            dd_get_constraint_range(cr->dd, &start, &natoms);
        }
        else
        {
            natoms = cr->dd->nat_home;
        }
    }
    else
    {
        natoms = md->homenr;
    }
    at2con = make_at2con(0, natoms, idef->il, idef->iparams, bDynamics,
                         &nflexcon);

    ncon_tot = idef->il[F_CONSTR].nr/3;

    /* Ensure we have enough padding for aligned loads for each thread */
    if (ncon_tot + li->ntask*simd_width > li->nc_alloc || li->nc_alloc == 0)
    {
        li->nc_alloc = over_alloc_dd(ncon_tot + li->ntask*simd_width);
        srenew(li->con_index, li->nc_alloc);
        resize_real_aligned(&li->bllen0, li->nc_alloc);
        resize_real_aligned(&li->ddist, li->nc_alloc);
        srenew(li->bla, 2*li->nc_alloc);
        resize_real_aligned(&li->blc, li->nc_alloc);
        resize_real_aligned(&li->blc1, li->nc_alloc);
        srenew(li->blnr, li->nc_alloc + 1);
        resize_real_aligned(&li->bllen, li->nc_alloc);
        srenew(li->tmpv, li->nc_alloc);
        if (DOMAINDECOMP(cr))
        {
            srenew(li->nlocat, li->nc_alloc);
        }
        resize_real_aligned(&li->tmp1, li->nc_alloc);
        resize_real_aligned(&li->tmp2, li->nc_alloc);
        resize_real_aligned(&li->tmp3, li->nc_alloc);
        resize_real_aligned(&li->tmp4, li->nc_alloc);
        resize_real_aligned(&li->mlambda, li->nc_alloc);
    }

    iatom = idef->il[F_CONSTR].iatoms;

    ncc_alloc_old = li->ncc_alloc;
    li->blnr[0]   = li->ncc;

    /* Assign the constraints for li->ntask LINCS tasks.
     * We target a uniform distribution of constraints over the tasks.
     * Note that when flexible constraints are present, but are removed here
     * (e.g. because we are doing EM) we get imbalance, but since that doesn't
     * happen during normal MD, that's ok.
     */
    int ncon_assign, ncon_target, con, th;

    /* Determine the number of constraints we need to assign here */
    ncon_assign      = ncon_tot;
    if (!bDynamics)
    {
        /* With energy minimization, flexible constraints are ignored
         * (and thus minimized, as they should be).
         */
        ncon_assign -= nflexcon;
    }

    /* Set the target constraint count per task to exactly uniform,
     * this might be overridden below.
     */
    ncon_target = (ncon_assign + li->ntask - 1)/li->ntask;

    /* Mark all constraints as unassigned by setting their index to -1 */
    for (con = 0; con < ncon_tot; con++)
    {
        li->con_index[con] = -1;
    }

    con = 0;
    for (th = 0; th < li->ntask; th++)
    {
        lincs_task_t *li_task;

        li_task = &li->task[th];

#ifdef LINCS_SIMD
        /* With indepedent tasks we likely have H-bond constraints or constraint
         * pairs. The connected constraints will be pulled into the task, so the
         * constraints per task will often exceed ncon_target.
         * Triangle constraints can also increase the count, but there are
         * relatively few of those, so we usually expect to get ncon_target.
         */
        if (li->bTaskDep)
        {
            /* We round ncon_target to a multiple of GMX_SIMD_WIDTH,
             * since otherwise a lot of operations can be wasted.
             * There are several ways to round here, we choose the one
             * that alternates block sizes, which helps with Intel HT.
             */
            ncon_target = ((ncon_assign*(th + 1))/li->ntask - li->nc_real + GMX_SIMD_REAL_WIDTH - 1) & ~(GMX_SIMD_REAL_WIDTH - 1);
        }
#endif

        /* Continue filling the arrays where we left off with the previous task,
         * including padding for SIMD.
         */
        li_task->b0 = li->nc;

        while (con < ncon_tot && li->nc - li_task->b0 < ncon_target)
        {
            if (li->con_index[con] == -1)
            {
                int  type, a1, a2;
                real lenA, lenB;

                type   = iatom[3*con];
                a1     = iatom[3*con + 1];
                a2     = iatom[3*con + 2];
                lenA   = idef->iparams[type].constr.dA;
                lenB   = idef->iparams[type].constr.dB;
                /* Skip the flexible constraints when not doing dynamics */
                if (bDynamics || lenA != 0 || lenB != 0)
                {
                    assign_constraint(li, con, a1, a2, lenA, lenB, &at2con);

                    if (li->ntask > 1 && !li->bTaskDep)
                    {
                        /* We can generate independent tasks. Check if we
                         * need to assign connected constraints to our task.
                         */
                        check_assign_connected(li, iatom, idef, bDynamics,
                                               a1, a2, &at2con);
                    }
                    if (li->ntask > 1 && li->ncg_triangle > 0)
                    {
                        /* Ensure constraints in one triangle are assigned
                         * to the same task.
                         */
                        check_assign_triangle(li, iatom, idef, bDynamics,
                                              con, a1, a2, &at2con);
                    }
                }
            }

            con++;
        }

        li_task->b1 = li->nc;

        if (simd_width > 1)
        {
            /* Copy the last atom pair indices and lengths for constraints
             * up to a multiple of simd_width, such that we can do all
             * SIMD operations without having to worry about end effects.
             */
            int i, last;

            li->nc = ((li_task->b1 + simd_width - 1)/simd_width)*simd_width;
            last   = li_task->b1 - 1;
            for (i = li_task->b1; i < li->nc; i++)
            {
                li->bla[i*2    ] = li->bla[last*2    ];
                li->bla[i*2 + 1] = li->bla[last*2 + 1];
                li->bllen0[i]    = li->bllen0[last];
                li->ddist[i]     = li->ddist[last];
                li->bllen[i]     = li->bllen[last];
                li->blnr[i + 1]  = li->blnr[last + 1];
            }
        }

        /* Keep track of how many constraints we assigned */
        li->nc_real += li_task->b1 - li_task->b0;

        if (debug)
        {
            fprintf(debug, "LINCS task %d constraints %d - %d\n",
                    th, li_task->b0, li_task->b1);
        }
    }

    assert(li->nc_real == ncon_assign);

    gmx_bool bSortMatrix;

    /* Without DD we order the blbnb matrix to optimize memory access.
     * With DD the overhead of sorting is more than the gain during access.
     */
    bSortMatrix = !DOMAINDECOMP(cr);

    if (li->ncc > li->ncc_alloc)
    {
        li->ncc_alloc = over_alloc_small(li->ncc);
        srenew(li->blbnb, li->ncc_alloc);
    }

#pragma omp parallel for num_threads(li->ntask) schedule(static)
    for (th = 0; th < li->ntask; th++)
    {
        lincs_task_t *li_task;

        li_task = &li->task[th];

        if (li->ncg_triangle > 0 &&
            li_task->b1 - li_task->b0 > li_task->tri_alloc)
        {
            /* This is allocating too much, but it is difficult to improve */
            li_task->tri_alloc = over_alloc_dd(li_task->b1 - li_task->b0);
            srenew(li_task->triangle, li_task->tri_alloc);
            srenew(li_task->tri_bits, li_task->tri_alloc);
        }

        set_matrix_indices(li, li_task, &at2con, bSortMatrix);
    }

    done_blocka(&at2con);

    if (cr->dd == NULL)
    {
        /* Since the matrix is static, we should free some memory */
        li->ncc_alloc = li->ncc;
        srenew(li->blbnb, li->ncc_alloc);
    }

    if (li->ncc_alloc > ncc_alloc_old)
    {
        srenew(li->blmf, li->ncc_alloc);
        srenew(li->blmf1, li->ncc_alloc);
        srenew(li->tmpncc, li->ncc_alloc);
    }

    if (DOMAINDECOMP(cr) && dd_constraints_nlocalatoms(cr->dd) != NULL)
    {
        int *nlocat_dd;

        nlocat_dd = dd_constraints_nlocalatoms(cr->dd);

        /* Convert nlocat from local topology to LINCS constraint indexing */
        for (con = 0; con < ncon_tot; con++)
        {
            li->nlocat[li->con_index[con]] = nlocat_dd[con];
        }
    }
    else
    {
        li->nlocat = NULL;
    }

    if (debug)
    {
        fprintf(debug, "Number of constraints is %d, padded %d, couplings %d\n",
                li->nc_real, li->nc, li->ncc);
    }

    if (li->ntask > 1)
    {
        lincs_thread_setup(li, md->nr);
    }

    set_lincs_matrix(li, md->invmass, md->lambda);
}

static void lincs_warning(FILE *fplog,
                          gmx_domdec_t *dd, rvec *x, rvec *xprime, t_pbc *pbc,
                          int ncons, int *bla, real *bllen, real wangle,
                          int maxwarn, int *warncount)
{
    int  b, i, j;
    rvec v0, v1;
    real wfac, d0, d1, cosine;
    char buf[STRLEN];

    wfac = cos(DEG2RAD*wangle);

    sprintf(buf, "bonds that rotated more than %g degrees:\n"
            " atom 1 atom 2  angle  previous, current, constraint length\n",
            wangle);
    fprintf(stderr, "%s", buf);
    if (fplog)
    {
        fprintf(fplog, "%s", buf);
    }

    for (b = 0; b < ncons; b++)
    {
        i = bla[2*b];
        j = bla[2*b+1];
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[i], x[j], v0);
            pbc_dx_aiuc(pbc, xprime[i], xprime[j], v1);
        }
        else
        {
            rvec_sub(x[i], x[j], v0);
            rvec_sub(xprime[i], xprime[j], v1);
        }
        d0     = norm(v0);
        d1     = norm(v1);
        cosine = iprod(v0, v1)/(d0*d1);
        if (cosine < wfac)
        {
            sprintf(buf, " %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
                    ddglatnr(dd, i), ddglatnr(dd, j),
                    RAD2DEG*acos(cosine), d0, d1, bllen[b]);
            fprintf(stderr, "%s", buf);
            if (fplog)
            {
                fprintf(fplog, "%s", buf);
            }
            if (!gmx_isfinite(d1))
            {
                gmx_fatal(FARGS, "Bond length not finite.");
            }

            (*warncount)++;
        }
    }
    if (*warncount > maxwarn)
    {
        too_many_constraint_warnings(econtLINCS, *warncount);
    }
}

static void cconerr(const struct gmx_lincsdata *lincsd,
                    rvec *x, t_pbc *pbc,
                    real *ncons_loc, real *ssd, real *max, int *imax)
{
    const int  *bla, *nlocat;
    const real *bllen;
    real        ma, ssd2;
    int         count, im, task;

    bla    = lincsd->bla;
    bllen  = lincsd->bllen;
    nlocat = lincsd->nlocat;

    ma    = 0;
    ssd2  = 0;
    im    = 0;
    count = 0;
    for (task = 0; task < lincsd->ntask; task++)
    {
        int b;

        for (b = lincsd->task[task].b0; b < lincsd->task[task].b1; b++)
        {
            real len, d, r2;
            rvec dx;

            if (pbc)
            {
                pbc_dx_aiuc(pbc, x[bla[2*b]], x[bla[2*b+1]], dx);
            }
            else
            {
                rvec_sub(x[bla[2*b]], x[bla[2*b+1]], dx);
            }
            r2  = norm2(dx);
            len = r2*gmx_invsqrt(r2);
            d   = fabs(len/bllen[b]-1);
            if (d > ma && (nlocat == NULL || nlocat[b]))
            {
                ma = d;
                im = b;
            }
            if (nlocat == NULL)
            {
                ssd2 += d*d;
                count++;
            }
            else
            {
                ssd2  += nlocat[b]*d*d;
                count += nlocat[b];
            }
        }
    }

    *ncons_loc = (nlocat ? 0.5 : 1)*count;
    *ssd       = (nlocat ? 0.5 : 1)*ssd2;
    *max       = ma;
    *imax      = im;
}

gmx_bool constrain_lincs(FILE *fplog, gmx_bool bLog, gmx_bool bEner,
                         t_inputrec *ir,
                         gmx_int64_t step,
                         struct gmx_lincsdata *lincsd, t_mdatoms *md,
                         t_commrec *cr,
                         rvec *x, rvec *xprime, rvec *min_proj,
                         matrix box, t_pbc *pbc,
                         real lambda, real *dvdlambda,
                         real invdt, rvec *v,
                         gmx_bool bCalcVir, tensor vir_r_m_dr,
                         int econq,
                         t_nrnb *nrnb,
                         int maxwarn, int *warncount)
{
    gmx_bool  bCalcDHDL;
    char      buf[STRLEN], buf2[22], buf3[STRLEN];
    int       i, p_imax;
    real      ncons_loc, p_ssd, p_max = 0;
    rvec      dx;
    gmx_bool  bOK, bWarn;

    bOK = TRUE;

    /* This boolean should be set by a flag passed to this routine.
     * We can also easily check if any constraint length is changed,
     * if not dH/dlambda=0 and we can also set the boolean to FALSE.
     */
    bCalcDHDL = (ir->efep != efepNO && dvdlambda != NULL);

    if (lincsd->nc == 0 && cr->dd == NULL)
    {
        if (bLog || bEner)
        {
            lincsd->rmsd_data[0] = 0;
            if (ir->eI == eiSD2 && v == NULL)
            {
                i = 2;
            }
            else
            {
                i = 1;
            }
            lincsd->rmsd_data[i] = 0;
        }

        return bOK;
    }

    if (econq == econqCoord)
    {
        /* We can't use bCalcDHDL here, since NULL can be passed for dvdlambda
         * also with efep!=fepNO.
         */
        if (ir->efep != efepNO)
        {
            if (md->nMassPerturbed && lincsd->matlam != md->lambda)
            {
                set_lincs_matrix(lincsd, md->invmass, md->lambda);
            }

            for (i = 0; i < lincsd->nc; i++)
            {
                lincsd->bllen[i] = lincsd->bllen0[i] + lambda*lincsd->ddist[i];
            }
        }

        if (lincsd->ncg_flex)
        {
            /* Set the flexible constraint lengths to the old lengths */
            if (pbc != NULL)
            {
                for (i = 0; i < lincsd->nc; i++)
                {
                    if (lincsd->bllen[i] == 0)
                    {
                        pbc_dx_aiuc(pbc, x[lincsd->bla[2*i]], x[lincsd->bla[2*i+1]], dx);
                        lincsd->bllen[i] = norm(dx);
                    }
                }
            }
            else
            {
                for (i = 0; i < lincsd->nc; i++)
                {
                    if (lincsd->bllen[i] == 0)
                    {
                        lincsd->bllen[i] =
                            sqrt(distance2(x[lincsd->bla[2*i]],
                                           x[lincsd->bla[2*i+1]]));
                    }
                }
            }
        }

        if (bLog && fplog)
        {
            cconerr(lincsd, xprime, pbc,
                    &ncons_loc, &p_ssd, &p_max, &p_imax);
        }

        /* This bWarn var can be updated by multiple threads
         * at the same time. But as we only need to detect
         * if a warning occured or not, this is not an issue.
         */
        bWarn = FALSE;

        /* The OpenMP parallel region of constrain_lincs for coords */
#pragma omp parallel num_threads(lincsd->ntask)
        {
            int th = gmx_omp_get_thread_num();

            clear_mat(lincsd->task[th].vir_r_m_dr);

            do_lincs(x, xprime, box, pbc, lincsd, th,
                     md->invmass, cr,
                     bCalcDHDL,
                     ir->LincsWarnAngle, &bWarn,
                     invdt, v, bCalcVir,
                     th == 0 ? vir_r_m_dr : lincsd->task[th].vir_r_m_dr);
        }

        if (bLog && fplog && lincsd->nc > 0)
        {
            fprintf(fplog, "   Rel. Constraint Deviation:  RMS         MAX     between atoms\n");
            fprintf(fplog, "       Before LINCS          %.6f    %.6f %6d %6d\n",
                    sqrt(p_ssd/ncons_loc), p_max,
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
        }
        if (bLog || bEner)
        {
            cconerr(lincsd, xprime, pbc,
                    &ncons_loc, &p_ssd, &p_max, &p_imax);
            /* Check if we are doing the second part of SD */
            if (ir->eI == eiSD2 && v == NULL)
            {
                i = 2;
            }
            else
            {
                i = 1;
            }
            lincsd->rmsd_data[0] = ncons_loc;
            lincsd->rmsd_data[i] = p_ssd;
        }
        else
        {
            lincsd->rmsd_data[0] = 0;
            lincsd->rmsd_data[1] = 0;
            lincsd->rmsd_data[2] = 0;
        }
        if (bLog && fplog && lincsd->nc > 0)
        {
            fprintf(fplog,
                    "        After LINCS          %.6f    %.6f %6d %6d\n\n",
                    sqrt(p_ssd/ncons_loc), p_max,
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                    ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
        }

        if (bWarn)
        {
            if (maxwarn >= 0)
            {
                cconerr(lincsd, xprime, pbc,
                        &ncons_loc, &p_ssd, &p_max, &p_imax);
                if (MULTISIM(cr))
                {
                    sprintf(buf3, " in simulation %d", cr->ms->sim);
                }
                else
                {
                    buf3[0] = 0;
                }
                sprintf(buf, "\nStep %s, time %g (ps)  LINCS WARNING%s\n"
                        "relative constraint deviation after LINCS:\n"
                        "rms %.6f, max %.6f (between atoms %d and %d)\n",
                        gmx_step_str(step, buf2), ir->init_t+step*ir->delta_t,
                        buf3,
                        sqrt(p_ssd/ncons_loc), p_max,
                        ddglatnr(cr->dd, lincsd->bla[2*p_imax]),
                        ddglatnr(cr->dd, lincsd->bla[2*p_imax+1]));
                if (fplog)
                {
                    fprintf(fplog, "%s", buf);
                }
                fprintf(stderr, "%s", buf);
                lincs_warning(fplog, cr->dd, x, xprime, pbc,
                              lincsd->nc, lincsd->bla, lincsd->bllen,
                              ir->LincsWarnAngle, maxwarn, warncount);
            }
            bOK = (p_max < 0.5);
        }

        if (lincsd->ncg_flex)
        {
            for (i = 0; (i < lincsd->nc); i++)
            {
                if (lincsd->bllen0[i] == 0 && lincsd->ddist[i] == 0)
                {
                    lincsd->bllen[i] = 0;
                }
            }
        }
    }
    else
    {
        /* The OpenMP parallel region of constrain_lincs for derivatives */
#pragma omp parallel num_threads(lincsd->ntask)
        {
            int th = gmx_omp_get_thread_num();

            do_lincsp(x, xprime, min_proj, pbc, lincsd, th,
                      md->invmass, econq, bCalcDHDL,
                      bCalcVir, th == 0 ? vir_r_m_dr : lincsd->task[th].vir_r_m_dr);
        }
    }

    if (bCalcDHDL)
    {
        /* Reduce the dH/dlambda contributions over the threads */
        real dhdlambda;
        int  th;

        dhdlambda = 0;
        for (th = 0; th < lincsd->ntask; th++)
        {
            dhdlambda += lincsd->task[th].dhdlambda;
        }
        if (econqCoord)
        {
            /* dhdlambda contains dH/dlambda*dt^2, correct for this */
            /* TODO This should probably use invdt, so that sd integrator scaling works properly */
            dhdlambda /= ir->delta_t*ir->delta_t;
        }
        *dvdlambda += dhdlambda;
    }

    if (bCalcVir && lincsd->ntask > 1)
    {
        for (i = 1; i < lincsd->ntask; i++)
        {
            m_add(vir_r_m_dr, lincsd->task[i].vir_r_m_dr, vir_r_m_dr);
        }
    }

    /* count assuming nit=1 */
    inc_nrnb(nrnb, eNR_LINCS, lincsd->nc_real);
    inc_nrnb(nrnb, eNR_LINCSMAT, (2+lincsd->nOrder)*lincsd->ncc);
    if (lincsd->ntriangle > 0)
    {
        inc_nrnb(nrnb, eNR_LINCSMAT, lincsd->nOrder*lincsd->ncc_triangle);
    }
    if (v)
    {
        inc_nrnb(nrnb, eNR_CONSTR_V, lincsd->nc_real*2);
    }
    if (bCalcVir)
    {
        inc_nrnb(nrnb, eNR_CONSTR_VIR, lincsd->nc_real);
    }

    return bOK;
}
