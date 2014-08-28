/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "typedefs.h"
#include "txtdump.h"
#include "vec.h"
#include "gromacs/utility/smalloc.h"
#include "coulomb.h"
#include "gmx_fatal.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "macros.h"

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"

/* Include the SIMD macro file and then check for support */
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#ifdef GMX_SIMD_HAVE_REAL
/* Turn on arbitrary width SIMD intrinsics for PME solve */
#    define PME_SIMD_SOLVE
#endif

#define PME_GRID_QA    0 /* Gridindex for A-state for Q */
#define PME_GRID_C6A   2 /* Gridindex for A-state for LJ */
#define DO_Q           2 /* Electrostatic grids have index q<2 */
#define DO_Q_AND_LJ    4 /* non-LB LJ grids have index 2 <= q < 4 */
#define DO_Q_AND_LJ_LB 9 /* With LB rules we need a total of 2+7 grids */

/* Pascal triangle coefficients scaled with (1/2)^6 for LJ-PME with LB-rules */
const real lb_scale_factor[] = {
    1.0/64, 6.0/64, 15.0/64, 20.0/64,
    15.0/64, 6.0/64, 1.0/64
};

/* Pascal triangle coefficients used in solve_pme_lj_yzx, only need to do 4 calculations due to symmetry */
const real lb_scale_factor_symm[] = { 2.0/64, 12.0/64, 30.0/64, 20.0/64 };

/* Check if we have 4-wide SIMD macro support */
#if (defined GMX_SIMD4_HAVE_REAL)
/* Do PME spread and gather with 4-wide SIMD.
 * NOTE: SIMD is only used with PME order 4 and 5 (which are the most common).
 */
#    define PME_SIMD4_SPREAD_GATHER

#    if (defined GMX_SIMD_HAVE_LOADU) && (defined GMX_SIMD_HAVE_STOREU)
/* With PME-order=4 on x86, unaligned load+store is slightly faster
 * than doubling all SIMD operations when using aligned load+store.
 */
#        define PME_SIMD4_UNALIGNED
#    endif
#endif

#define DFT_TOL 1e-7
/* #define PRT_FORCE */
/* conditions for on the fly time-measurement */
/* #define TAKETIME (step > 1 && timesteps < 10) */
#define TAKETIME FALSE

/* #define PME_TIME_THREADS */

#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
#    define SIMD4_ALIGNMENT  (GMX_SIMD4_WIDTH*sizeof(real))
#else
/* We can use any alignment, apart from 0, so we use 4 reals */
#    define SIMD4_ALIGNMENT  (4*sizeof(real))
#endif

/* GMX_CACHE_SEP should be a multiple of the SIMD and SIMD4 register size
 * to preserve alignment.
 */
#define GMX_CACHE_SEP 64

/* We only define a maximum to be able to use local arrays without allocation.
 * An order larger than 12 should never be needed, even for test cases.
 * If needed it can be changed here.
 */
#define PME_ORDER_MAX 12

/* Internal datastructures */
typedef struct {
    int send_index0;
    int send_nindex;
    int recv_index0;
    int recv_nindex;
    int recv_size;   /* Receive buffer width, used with OpenMP */
} pme_grid_comm_t;

typedef struct {
#ifdef GMX_MPI
    MPI_Comm         mpi_comm;
#endif
    int              nnodes, nodeid;
    int             *s2g0;
    int             *s2g1;
    int              noverlap_nodes;
    int             *send_id, *recv_id;
    int              send_size; /* Send buffer width, used with OpenMP */
    pme_grid_comm_t *comm_data;
    real            *sendbuf;
    real            *recvbuf;
} pme_overlap_t;

typedef struct {
    int *n;      /* Cumulative counts of the number of particles per thread */
    int  nalloc; /* Allocation size of i */
    int *i;      /* Particle indices ordered on thread index (n) */
} thread_plist_t;

typedef struct {
    int      *thread_one;
    int       n;
    int      *ind;
    splinevec theta;
    real     *ptr_theta_z;
    splinevec dtheta;
    real     *ptr_dtheta_z;
} splinedata_t;

typedef struct {
    int      dimind;        /* The index of the dimension, 0=x, 1=y */
    int      nslab;
    int      nodeid;
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif

    int     *node_dest;     /* The nodes to send x and q to with DD */
    int     *node_src;      /* The nodes to receive x and q from with DD */
    int     *buf_index;     /* Index for commnode into the buffers */

    int      maxshift;

    int      npd;
    int      pd_nalloc;
    int     *pd;
    int     *count;         /* The number of atoms to send to each node */
    int    **count_thread;
    int     *rcount;        /* The number of atoms to receive */

    int      n;
    int      nalloc;
    rvec    *x;
    real    *coefficient;
    rvec    *f;
    gmx_bool bSpread;       /* These coordinates are used for spreading */
    int      pme_order;
    ivec    *idx;
    rvec    *fractx;            /* Fractional coordinate relative to
                                 * the lower cell boundary
                                 */
    int             nthread;
    int            *thread_idx; /* Which thread should spread which coefficient */
    thread_plist_t *thread_plist;
    splinedata_t   *spline;
} pme_atomcomm_t;

#define FLBS  3
#define FLBSZ 4

typedef struct {
    ivec  ci;     /* The spatial location of this grid         */
    ivec  n;      /* The used size of *grid, including order-1 */
    ivec  offset; /* The grid offset from the full node grid   */
    int   order;  /* PME spreading order                       */
    ivec  s;      /* The allocated size of *grid, s >= n       */
    real *grid;   /* The grid local thread, size n             */
} pmegrid_t;

typedef struct {
    pmegrid_t  grid;         /* The full node grid (non thread-local)            */
    int        nthread;      /* The number of threads operating on this grid     */
    ivec       nc;           /* The local spatial decomposition over the threads */
    pmegrid_t *grid_th;      /* Array of grids for each thread                   */
    real      *grid_all;     /* Allocated array for the grids in *grid_th        */
    int      **g2t;          /* The grid to thread index                         */
    ivec       nthread_comm; /* The number of threads to communicate with        */
} pmegrids_t;

typedef struct {
#ifdef PME_SIMD4_SPREAD_GATHER
    /* Masks for 4-wide SIMD aligned spreading and gathering */
    gmx_simd4_bool_t mask_S0[6], mask_S1[6];
#else
    int              dummy; /* C89 requires that struct has at least one member */
#endif
} pme_spline_work_t;

typedef struct {
    /* work data for solve_pme */
    int      nalloc;
    real *   mhx;
    real *   mhy;
    real *   mhz;
    real *   m2;
    real *   denom;
    real *   tmp1_alloc;
    real *   tmp1;
    real *   tmp2;
    real *   eterm;
    real *   m2inv;

    real     energy_q;
    matrix   vir_q;
    real     energy_lj;
    matrix   vir_lj;
} pme_work_t;

typedef struct gmx_pme {
    int           ndecompdim; /* The number of decomposition dimensions */
    int           nodeid;     /* Our nodeid in mpi->mpi_comm */
    int           nodeid_major;
    int           nodeid_minor;
    int           nnodes;    /* The number of nodes doing PME */
    int           nnodes_major;
    int           nnodes_minor;

    MPI_Comm      mpi_comm;
    MPI_Comm      mpi_comm_d[2]; /* Indexed on dimension, 0=x, 1=y */
#ifdef GMX_MPI
    MPI_Datatype  rvec_mpi;      /* the pme vector's MPI type */
#endif

    gmx_bool   bUseThreads;   /* Does any of the PME ranks have nthread>1 ?  */
    int        nthread;       /* The number of threads doing PME on our rank */

    gmx_bool   bPPnode;       /* Node also does particle-particle forces */
    gmx_bool   bFEP;          /* Compute Free energy contribution */
    gmx_bool   bFEP_q;
    gmx_bool   bFEP_lj;
    int        nkx, nky, nkz; /* Grid dimensions */
    gmx_bool   bP3M;          /* Do P3M: optimize the influence function */
    int        pme_order;
    real       epsilon_r;

    int        ljpme_combination_rule;  /* Type of combination rule in LJ-PME */

    int        ngrids;                  /* number of grids we maintain for pmegrid, (c)fftgrid and pfft_setups*/

    pmegrids_t pmegrid[DO_Q_AND_LJ_LB]; /* Grids on which we do spreading/interpolation,
                                         * includes overlap Grid indices are ordered as
                                         * follows:
                                         * 0: Coloumb PME, state A
                                         * 1: Coloumb PME, state B
                                         * 2-8: LJ-PME
                                         * This can probably be done in a better way
                                         * but this simple hack works for now
                                         */
    /* The PME coefficient spreading grid sizes/strides, includes pme_order-1 */
    int        pmegrid_nx, pmegrid_ny, pmegrid_nz;
    /* pmegrid_nz might be larger than strictly necessary to ensure
     * memory alignment, pmegrid_nz_base gives the real base size.
     */
    int     pmegrid_nz_base;
    /* The local PME grid starting indices */
    int     pmegrid_start_ix, pmegrid_start_iy, pmegrid_start_iz;

    /* Work data for spreading and gathering */
    pme_spline_work_t     *spline_work;

    real                 **fftgrid; /* Grids for FFT. With 1D FFT decomposition this can be a pointer */
    /* inside the interpolation grid, but separate for 2D PME decomp. */
    int                    fftgrid_nx, fftgrid_ny, fftgrid_nz;

    t_complex            **cfftgrid;  /* Grids for complex FFT data */

    int                    cfftgrid_nx, cfftgrid_ny, cfftgrid_nz;

    gmx_parallel_3dfft_t  *pfft_setup;

    int                   *nnx, *nny, *nnz;
    real                  *fshx, *fshy, *fshz;

    pme_atomcomm_t         atc[2]; /* Indexed on decomposition index */
    matrix                 recipbox;
    splinevec              bsp_mod;
    /* Buffers to store data for local atoms for L-B combination rule
     * calculations in LJ-PME. lb_buf1 stores either the coefficients
     * for spreading/gathering (in serial), or the C6 coefficient for
     * local atoms (in parallel).  lb_buf2 is only used in parallel,
     * and stores the sigma values for local atoms. */
    real                 *lb_buf1, *lb_buf2;
    int                   lb_buf_nalloc; /* Allocation size for the above buffers. */

    pme_overlap_t         overlap[2];    /* Indexed on dimension, 0=x, 1=y */

    pme_atomcomm_t        atc_energy;    /* Only for gmx_pme_calc_energy */

    rvec                 *bufv;          /* Communication buffer */
    real                 *bufr;          /* Communication buffer */
    int                   buf_nalloc;    /* The communication buffer size */

    /* thread local work data for solve_pme */
    pme_work_t *work;

    /* Work data for sum_qgrid */
    real *   sum_qgrid_tmp;
    real *   sum_qgrid_dd_tmp;
} t_gmx_pme;

static void calc_interpolation_idx(gmx_pme_t pme, pme_atomcomm_t *atc,
                                   int start, int grid_index, int end, int thread)
{
    int             i;
    int            *idxptr, tix, tiy, tiz;
    real           *xptr, *fptr, tx, ty, tz;
    real            rxx, ryx, ryy, rzx, rzy, rzz;
    int             nx, ny, nz;
    int             start_ix, start_iy, start_iz;
    int            *g2tx, *g2ty, *g2tz;
    gmx_bool        bThreads;
    int            *thread_idx = NULL;
    thread_plist_t *tpl        = NULL;
    int            *tpl_n      = NULL;
    int             thread_i;

    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;

    start_ix = pme->pmegrid_start_ix;
    start_iy = pme->pmegrid_start_iy;
    start_iz = pme->pmegrid_start_iz;

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    g2tx = pme->pmegrid[grid_index].g2t[XX];
    g2ty = pme->pmegrid[grid_index].g2t[YY];
    g2tz = pme->pmegrid[grid_index].g2t[ZZ];

    bThreads = (atc->nthread > 1);
    if (bThreads)
    {
        thread_idx = atc->thread_idx;

        tpl   = &atc->thread_plist[thread];
        tpl_n = tpl->n;
        for (i = 0; i < atc->nthread; i++)
        {
            tpl_n[i] = 0;
        }
    }

    for (i = start; i < end; i++)
    {
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];

        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );

        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fptr[XX] = tx - tix + pme->fshx[tix];
        fptr[YY] = ty - tiy + pme->fshy[tiy];
        fptr[ZZ] = tz - tiz;

        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];

#ifdef DEBUG
        range_check(idxptr[XX], 0, pme->pmegrid_nx);
        range_check(idxptr[YY], 0, pme->pmegrid_ny);
        range_check(idxptr[ZZ], 0, pme->pmegrid_nz);
#endif

        if (bThreads)
        {
            thread_i      = g2tx[idxptr[XX]] + g2ty[idxptr[YY]] + g2tz[idxptr[ZZ]];
            thread_idx[i] = thread_i;
            tpl_n[thread_i]++;
        }
    }

    if (bThreads)
    {
        /* Make a list of particle indices sorted on thread */

        /* Get the cumulative count */
        for (i = 1; i < atc->nthread; i++)
        {
            tpl_n[i] += tpl_n[i-1];
        }
        /* The current implementation distributes particles equally
         * over the threads, so we could actually allocate for that
         * in pme_realloc_atomcomm_things.
         */
        if (tpl_n[atc->nthread-1] > tpl->nalloc)
        {
            tpl->nalloc = over_alloc_large(tpl_n[atc->nthread-1]);
            srenew(tpl->i, tpl->nalloc);
        }
        /* Set tpl_n to the cumulative start */
        for (i = atc->nthread-1; i >= 1; i--)
        {
            tpl_n[i] = tpl_n[i-1];
        }
        tpl_n[0] = 0;

        /* Fill our thread local array with indices sorted on thread */
        for (i = start; i < end; i++)
        {
            tpl->i[tpl_n[atc->thread_idx[i]]++] = i;
        }
        /* Now tpl_n contains the cummulative count again */
    }
}

static void make_thread_local_ind(pme_atomcomm_t *atc,
                                  int thread, splinedata_t *spline)
{
    int             n, t, i, start, end;
    thread_plist_t *tpl;

    /* Combine the indices made by each thread into one index */

    n     = 0;
    start = 0;
    for (t = 0; t < atc->nthread; t++)
    {
        tpl = &atc->thread_plist[t];
        /* Copy our part (start - end) from the list of thread t */
        if (thread > 0)
        {
            start = tpl->n[thread-1];
        }
        end = tpl->n[thread];
        for (i = start; i < end; i++)
        {
            spline->ind[n++] = tpl->i[i];
        }
    }

    spline->n = n;
}


static void pme_calc_pidx(int start, int end,
                          matrix recipbox, rvec x[],
                          pme_atomcomm_t *atc, int *count)
{
    int   nslab, i;
    int   si;
    real *xptr, s;
    real  rxx, ryx, rzx, ryy, rzy;
    int  *pd;

    /* Calculate PME task index (pidx) for each grid index.
     * Here we always assign equally sized slabs to each node
     * for load balancing reasons (the PME grid spacing is not used).
     */

    nslab = atc->nslab;
    pd    = atc->pd;

    /* Reset the count */
    for (i = 0; i < nslab; i++)
    {
        count[i] = 0;
    }

    if (atc->dimind == 0)
    {
        rxx = recipbox[XX][XX];
        ryx = recipbox[YY][XX];
        rzx = recipbox[ZZ][XX];
        /* Calculate the node index in x-dimension */
        for (i = start; i < end; i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab*(xptr[XX]*rxx + xptr[YY]*ryx + xptr[ZZ]*rzx);
            si    = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
    else
    {
        ryy = recipbox[YY][YY];
        rzy = recipbox[ZZ][YY];
        /* Calculate the node index in y-dimension */
        for (i = start; i < end; i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab*(xptr[YY]*ryy + xptr[ZZ]*rzy);
            si    = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
}

static void pme_calc_pidx_wrapper(int natoms, matrix recipbox, rvec x[],
                                  pme_atomcomm_t *atc)
{
    int nthread, thread, slab;

    nthread = atc->nthread;

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        pme_calc_pidx(natoms* thread   /nthread,
                      natoms*(thread+1)/nthread,
                      recipbox, x, atc, atc->count_thread[thread]);
    }
    /* Non-parallel reduction, since nslab is small */

    for (thread = 1; thread < nthread; thread++)
    {
        for (slab = 0; slab < atc->nslab; slab++)
        {
            atc->count_thread[0][slab] += atc->count_thread[thread][slab];
        }
    }
}

static void realloc_splinevec(splinevec th, real **ptr_z, int nalloc)
{
    const int padding = 4;
    int       i;

    srenew(th[XX], nalloc);
    srenew(th[YY], nalloc);
    /* In z we add padding, this is only required for the aligned SIMD code */
    sfree_aligned(*ptr_z);
    snew_aligned(*ptr_z, nalloc+2*padding, SIMD4_ALIGNMENT);
    th[ZZ] = *ptr_z + padding;

    for (i = 0; i < padding; i++)
    {
        (*ptr_z)[               i] = 0;
        (*ptr_z)[padding+nalloc+i] = 0;
    }
}

static void pme_realloc_splinedata(splinedata_t *spline, pme_atomcomm_t *atc)
{
    int i, d;

    srenew(spline->ind, atc->nalloc);
    /* Initialize the index to identity so it works without threads */
    for (i = 0; i < atc->nalloc; i++)
    {
        spline->ind[i] = i;
    }

    realloc_splinevec(spline->theta, &spline->ptr_theta_z,
                      atc->pme_order*atc->nalloc);
    realloc_splinevec(spline->dtheta, &spline->ptr_dtheta_z,
                      atc->pme_order*atc->nalloc);
}

static void pme_realloc_atomcomm_things(pme_atomcomm_t *atc)
{
    int nalloc_old, i, j, nalloc_tpl;

    /* We have to avoid a NULL pointer for atc->x to avoid
     * possible fatal errors in MPI routines.
     */
    if (atc->n > atc->nalloc || atc->nalloc == 0)
    {
        nalloc_old  = atc->nalloc;
        atc->nalloc = over_alloc_dd(max(atc->n, 1));

        if (atc->nslab > 1)
        {
            srenew(atc->x, atc->nalloc);
            srenew(atc->coefficient, atc->nalloc);
            srenew(atc->f, atc->nalloc);
            for (i = nalloc_old; i < atc->nalloc; i++)
            {
                clear_rvec(atc->f[i]);
            }
        }
        if (atc->bSpread)
        {
            srenew(atc->fractx, atc->nalloc);
            srenew(atc->idx, atc->nalloc);

            if (atc->nthread > 1)
            {
                srenew(atc->thread_idx, atc->nalloc);
            }

            for (i = 0; i < atc->nthread; i++)
            {
                pme_realloc_splinedata(&atc->spline[i], atc);
            }
        }
    }
}

static void pme_dd_sendrecv(pme_atomcomm_t gmx_unused *atc,
                            gmx_bool gmx_unused bBackward, int gmx_unused shift,
                            void gmx_unused *buf_s, int gmx_unused nbyte_s,
                            void gmx_unused *buf_r, int gmx_unused nbyte_r)
{
#ifdef GMX_MPI
    int        dest, src;
    MPI_Status stat;

    if (bBackward == FALSE)
    {
        dest = atc->node_dest[shift];
        src  = atc->node_src[shift];
    }
    else
    {
        dest = atc->node_src[shift];
        src  = atc->node_dest[shift];
    }

    if (nbyte_s > 0 && nbyte_r > 0)
    {
        MPI_Sendrecv(buf_s, nbyte_s, MPI_BYTE,
                     dest, shift,
                     buf_r, nbyte_r, MPI_BYTE,
                     src, shift,
                     atc->mpi_comm, &stat);
    }
    else if (nbyte_s > 0)
    {
        MPI_Send(buf_s, nbyte_s, MPI_BYTE,
                 dest, shift,
                 atc->mpi_comm);
    }
    else if (nbyte_r > 0)
    {
        MPI_Recv(buf_r, nbyte_r, MPI_BYTE,
                 src, shift,
                 atc->mpi_comm, &stat);
    }
#endif
}

static void dd_pmeredist_pos_coeffs(gmx_pme_t pme,
                                    int n, gmx_bool bX, rvec *x, real *data,
                                    pme_atomcomm_t *atc)
{
    int *commnode, *buf_index;
    int  nnodes_comm, i, nsend, local_pos, buf_pos, node, scount, rcount;

    commnode  = atc->node_dest;
    buf_index = atc->buf_index;

    nnodes_comm = min(2*atc->maxshift, atc->nslab-1);

    nsend = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        buf_index[commnode[i]] = nsend;
        nsend                 += atc->count[commnode[i]];
    }
    if (bX)
    {
        if (atc->count[atc->nodeid] + nsend != n)
        {
            gmx_fatal(FARGS, "%d particles communicated to PME rank %d are more than 2/3 times the cut-off out of the domain decomposition cell of their charge group in dimension %c.\n"
                      "This usually means that your system is not well equilibrated.",
                      n - (atc->count[atc->nodeid] + nsend),
                      pme->nodeid, 'x'+atc->dimind);
        }

        if (nsend > pme->buf_nalloc)
        {
            pme->buf_nalloc = over_alloc_dd(nsend);
            srenew(pme->bufv, pme->buf_nalloc);
            srenew(pme->bufr, pme->buf_nalloc);
        }

        atc->n = atc->count[atc->nodeid];
        for (i = 0; i < nnodes_comm; i++)
        {
            scount = atc->count[commnode[i]];
            /* Communicate the count */
            if (debug)
            {
                fprintf(debug, "dimind %d PME rank %d send to rank %d: %d\n",
                        atc->dimind, atc->nodeid, commnode[i], scount);
            }
            pme_dd_sendrecv(atc, FALSE, i,
                            &scount, sizeof(int),
                            &atc->rcount[i], sizeof(int));
            atc->n += atc->rcount[i];
        }

        pme_realloc_atomcomm_things(atc);
    }

    local_pos = 0;
    for (i = 0; i < n; i++)
    {
        node = atc->pd[i];
        if (node == atc->nodeid)
        {
            /* Copy direct to the receive buffer */
            if (bX)
            {
                copy_rvec(x[i], atc->x[local_pos]);
            }
            atc->coefficient[local_pos] = data[i];
            local_pos++;
        }
        else
        {
            /* Copy to the send buffer */
            if (bX)
            {
                copy_rvec(x[i], pme->bufv[buf_index[node]]);
            }
            pme->bufr[buf_index[node]] = data[i];
            buf_index[node]++;
        }
    }

    buf_pos = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        scount = atc->count[commnode[i]];
        rcount = atc->rcount[i];
        if (scount > 0 || rcount > 0)
        {
            if (bX)
            {
                /* Communicate the coordinates */
                pme_dd_sendrecv(atc, FALSE, i,
                                pme->bufv[buf_pos], scount*sizeof(rvec),
                                atc->x[local_pos], rcount*sizeof(rvec));
            }
            /* Communicate the coefficients */
            pme_dd_sendrecv(atc, FALSE, i,
                            pme->bufr+buf_pos, scount*sizeof(real),
                            atc->coefficient+local_pos, rcount*sizeof(real));
            buf_pos   += scount;
            local_pos += atc->rcount[i];
        }
    }
}

static void dd_pmeredist_f(gmx_pme_t pme, pme_atomcomm_t *atc,
                           int n, rvec *f,
                           gmx_bool bAddF)
{
    int *commnode, *buf_index;
    int  nnodes_comm, local_pos, buf_pos, i, scount, rcount, node;

    commnode  = atc->node_dest;
    buf_index = atc->buf_index;

    nnodes_comm = min(2*atc->maxshift, atc->nslab-1);

    local_pos = atc->count[atc->nodeid];
    buf_pos   = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        scount = atc->rcount[i];
        rcount = atc->count[commnode[i]];
        if (scount > 0 || rcount > 0)
        {
            /* Communicate the forces */
            pme_dd_sendrecv(atc, TRUE, i,
                            atc->f[local_pos], scount*sizeof(rvec),
                            pme->bufv[buf_pos], rcount*sizeof(rvec));
            local_pos += scount;
        }
        buf_index[commnode[i]] = buf_pos;
        buf_pos               += rcount;
    }

    local_pos = 0;
    if (bAddF)
    {
        for (i = 0; i < n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Add from the local force array */
                rvec_inc(f[i], atc->f[local_pos]);
                local_pos++;
            }
            else
            {
                /* Add from the receive buffer */
                rvec_inc(f[i], pme->bufv[buf_index[node]]);
                buf_index[node]++;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Copy from the local force array */
                copy_rvec(atc->f[local_pos], f[i]);
                local_pos++;
            }
            else
            {
                /* Copy from the receive buffer */
                copy_rvec(pme->bufv[buf_index[node]], f[i]);
                buf_index[node]++;
            }
        }
    }
}

#ifdef GMX_MPI
static void gmx_sum_qgrid_dd(gmx_pme_t pme, real *grid, int direction)
{
    pme_overlap_t *overlap;
    int            send_index0, send_nindex;
    int            recv_index0, recv_nindex;
    MPI_Status     stat;
    int            i, j, k, ix, iy, iz, icnt;
    int            ipulse, send_id, recv_id, datasize;
    real          *p;
    real          *sendptr, *recvptr;

    /* Start with minor-rank communication. This is a bit of a pain since it is not contiguous */
    overlap = &pme->overlap[1];

    for (ipulse = 0; ipulse < overlap->noverlap_nodes; ipulse++)
    {
        /* Since we have already (un)wrapped the overlap in the z-dimension,
         * we only have to communicate 0 to nkz (not pmegrid_nz).
         */
        if (direction == GMX_SUM_GRID_FORWARD)
        {
            send_id       = overlap->send_id[ipulse];
            recv_id       = overlap->recv_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            recv_index0   = overlap->comm_data[ipulse].recv_index0;
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
        }
        else
        {
            send_id       = overlap->recv_id[ipulse];
            recv_id       = overlap->send_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].recv_index0;
            send_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recv_index0   = overlap->comm_data[ipulse].send_index0;
            recv_nindex   = overlap->comm_data[ipulse].send_nindex;
        }

        /* Copy data to contiguous send buffer */
        if (debug)
        {
            fprintf(debug, "PME send rank %d %d -> %d grid start %d Communicating %d to %d\n",
                    pme->nodeid, overlap->nodeid, send_id,
                    pme->pmegrid_start_iy,
                    send_index0-pme->pmegrid_start_iy,
                    send_index0-pme->pmegrid_start_iy+send_nindex);
        }
        icnt = 0;
        for (i = 0; i < pme->pmegrid_nx; i++)
        {
            ix = i;
            for (j = 0; j < send_nindex; j++)
            {
                iy = j + send_index0 - pme->pmegrid_start_iy;
                for (k = 0; k < pme->nkz; k++)
                {
                    iz = k;
                    overlap->sendbuf[icnt++] = grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz];
                }
            }
        }

        datasize      = pme->pmegrid_nx * pme->nkz;

        MPI_Sendrecv(overlap->sendbuf, send_nindex*datasize, GMX_MPI_REAL,
                     send_id, ipulse,
                     overlap->recvbuf, recv_nindex*datasize, GMX_MPI_REAL,
                     recv_id, ipulse,
                     overlap->mpi_comm, &stat);

        /* Get data from contiguous recv buffer */
        if (debug)
        {
            fprintf(debug, "PME recv rank %d %d <- %d grid start %d Communicating %d to %d\n",
                    pme->nodeid, overlap->nodeid, recv_id,
                    pme->pmegrid_start_iy,
                    recv_index0-pme->pmegrid_start_iy,
                    recv_index0-pme->pmegrid_start_iy+recv_nindex);
        }
        icnt = 0;
        for (i = 0; i < pme->pmegrid_nx; i++)
        {
            ix = i;
            for (j = 0; j < recv_nindex; j++)
            {
                iy = j + recv_index0 - pme->pmegrid_start_iy;
                for (k = 0; k < pme->nkz; k++)
                {
                    iz = k;
                    if (direction == GMX_SUM_GRID_FORWARD)
                    {
                        grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz] += overlap->recvbuf[icnt++];
                    }
                    else
                    {
                        grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz]  = overlap->recvbuf[icnt++];
                    }
                }
            }
        }
    }

    /* Major dimension is easier, no copying required,
     * but we might have to sum to separate array.
     * Since we don't copy, we have to communicate up to pmegrid_nz,
     * not nkz as for the minor direction.
     */
    overlap = &pme->overlap[0];

    for (ipulse = 0; ipulse < overlap->noverlap_nodes; ipulse++)
    {
        if (direction == GMX_SUM_GRID_FORWARD)
        {
            send_id       = overlap->send_id[ipulse];
            recv_id       = overlap->recv_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            recv_index0   = overlap->comm_data[ipulse].recv_index0;
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recvptr       = overlap->recvbuf;
        }
        else
        {
            send_id       = overlap->recv_id[ipulse];
            recv_id       = overlap->send_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].recv_index0;
            send_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recv_index0   = overlap->comm_data[ipulse].send_index0;
            recv_nindex   = overlap->comm_data[ipulse].send_nindex;
            recvptr       = grid + (recv_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
        }

        sendptr       = grid + (send_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
        datasize      = pme->pmegrid_ny * pme->pmegrid_nz;

        if (debug)
        {
            fprintf(debug, "PME send rank %d %d -> %d grid start %d Communicating %d to %d\n",
                    pme->nodeid, overlap->nodeid, send_id,
                    pme->pmegrid_start_ix,
                    send_index0-pme->pmegrid_start_ix,
                    send_index0-pme->pmegrid_start_ix+send_nindex);
            fprintf(debug, "PME recv rank %d %d <- %d grid start %d Communicating %d to %d\n",
                    pme->nodeid, overlap->nodeid, recv_id,
                    pme->pmegrid_start_ix,
                    recv_index0-pme->pmegrid_start_ix,
                    recv_index0-pme->pmegrid_start_ix+recv_nindex);
        }

        MPI_Sendrecv(sendptr, send_nindex*datasize, GMX_MPI_REAL,
                     send_id, ipulse,
                     recvptr, recv_nindex*datasize, GMX_MPI_REAL,
                     recv_id, ipulse,
                     overlap->mpi_comm, &stat);

        /* ADD data from contiguous recv buffer */
        if (direction == GMX_SUM_GRID_FORWARD)
        {
            p = grid + (recv_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
            for (i = 0; i < recv_nindex*datasize; i++)
            {
                p[i] += overlap->recvbuf[i];
            }
        }
    }
}
#endif


static int copy_pmegrid_to_fftgrid(gmx_pme_t pme, real *pmegrid, real *fftgrid, int grid_index)
{
    ivec    local_fft_ndata, local_fft_offset, local_fft_size;
    ivec    local_pme_size;
    int     i, ix, iy, iz;
    int     pmeidx, fftidx;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    local_pme_size[0] = pme->pmegrid_nx;
    local_pme_size[1] = pme->pmegrid_ny;
    local_pme_size[2] = pme->pmegrid_nz;

    /* The fftgrid is always 'justified' to the lower-left corner of the PME grid,
       the offset is identical, and the PME grid always has more data (due to overlap)
     */
    {
#ifdef DEBUG_PME
        FILE *fp, *fp2;
        char  fn[STRLEN];
        real  val;
        sprintf(fn, "pmegrid%d.pdb", pme->nodeid);
        fp = gmx_ffopen(fn, "w");
        sprintf(fn, "pmegrid%d.txt", pme->nodeid);
        fp2 = gmx_ffopen(fn, "w");
#endif

        for (ix = 0; ix < local_fft_ndata[XX]; ix++)
        {
            for (iy = 0; iy < local_fft_ndata[YY]; iy++)
            {
                for (iz = 0; iz < local_fft_ndata[ZZ]; iz++)
                {
                    pmeidx          = ix*(local_pme_size[YY]*local_pme_size[ZZ])+iy*(local_pme_size[ZZ])+iz;
                    fftidx          = ix*(local_fft_size[YY]*local_fft_size[ZZ])+iy*(local_fft_size[ZZ])+iz;
                    fftgrid[fftidx] = pmegrid[pmeidx];
#ifdef DEBUG_PME
                    val = 100*pmegrid[pmeidx];
                    if (pmegrid[pmeidx] != 0)
                    {
                        gmx_fprintf_pdb_atomline(fp, epdbATOM, pmeidx, "CA", ' ', "GLY", ' ', pmeidx, ' ',
                                                 5.0*ix, 5.0*iy, 5.0*iz, 1.0, val, "");
                    }
                    if (pmegrid[pmeidx] != 0)
                    {
                        fprintf(fp2, "%-12s  %5d  %5d  %5d  %12.5e\n",
                                "qgrid",
                                pme->pmegrid_start_ix + ix,
                                pme->pmegrid_start_iy + iy,
                                pme->pmegrid_start_iz + iz,
                                pmegrid[pmeidx]);
                    }
#endif
                }
            }
        }
#ifdef DEBUG_PME
        gmx_ffclose(fp);
        gmx_ffclose(fp2);
#endif
    }
    return 0;
}


static gmx_cycles_t omp_cyc_start()
{
    return gmx_cycles_read();
}

static gmx_cycles_t omp_cyc_end(gmx_cycles_t c)
{
    return gmx_cycles_read() - c;
}


static int copy_fftgrid_to_pmegrid(gmx_pme_t pme, const real *fftgrid, real *pmegrid, int grid_index,
                                   int nthread, int thread)
{
    ivec          local_fft_ndata, local_fft_offset, local_fft_size;
    ivec          local_pme_size;
    int           ixy0, ixy1, ixy, ix, iy, iz;
    int           pmeidx, fftidx;
#ifdef PME_TIME_THREADS
    gmx_cycles_t  c1;
    static double cs1 = 0;
    static int    cnt = 0;
#endif

#ifdef PME_TIME_THREADS
    c1 = omp_cyc_start();
#endif
    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    local_pme_size[0] = pme->pmegrid_nx;
    local_pme_size[1] = pme->pmegrid_ny;
    local_pme_size[2] = pme->pmegrid_nz;

    /* The fftgrid is always 'justified' to the lower-left corner of the PME grid,
       the offset is identical, and the PME grid always has more data (due to overlap)
     */
    ixy0 = ((thread  )*local_fft_ndata[XX]*local_fft_ndata[YY])/nthread;
    ixy1 = ((thread+1)*local_fft_ndata[XX]*local_fft_ndata[YY])/nthread;

    for (ixy = ixy0; ixy < ixy1; ixy++)
    {
        ix = ixy/local_fft_ndata[YY];
        iy = ixy - ix*local_fft_ndata[YY];

        pmeidx = (ix*local_pme_size[YY] + iy)*local_pme_size[ZZ];
        fftidx = (ix*local_fft_size[YY] + iy)*local_fft_size[ZZ];
        for (iz = 0; iz < local_fft_ndata[ZZ]; iz++)
        {
            pmegrid[pmeidx+iz] = fftgrid[fftidx+iz];
        }
    }

#ifdef PME_TIME_THREADS
    c1   = omp_cyc_end(c1);
    cs1 += (double)c1;
    cnt++;
    if (cnt % 20 == 0)
    {
        printf("copy %.2f\n", cs1*1e-9);
    }
#endif

    return 0;
}


static void wrap_periodic_pmegrid(gmx_pme_t pme, real *pmegrid)
{
    int     nx, ny, nz, pnx, pny, pnz, ny_x, overlap, ix, iy, iz;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    pnx = pme->pmegrid_nx;
    pny = pme->pmegrid_ny;
    pnz = pme->pmegrid_nz;

    overlap = pme->pme_order - 1;

    /* Add periodic overlap in z */
    for (ix = 0; ix < pme->pmegrid_nx; ix++)
    {
        for (iy = 0; iy < pme->pmegrid_ny; iy++)
        {
            for (iz = 0; iz < overlap; iz++)
            {
                pmegrid[(ix*pny+iy)*pnz+iz] +=
                    pmegrid[(ix*pny+iy)*pnz+nz+iz];
            }
        }
    }

    if (pme->nnodes_minor == 1)
    {
        for (ix = 0; ix < pme->pmegrid_nx; ix++)
        {
            for (iy = 0; iy < overlap; iy++)
            {
                for (iz = 0; iz < nz; iz++)
                {
                    pmegrid[(ix*pny+iy)*pnz+iz] +=
                        pmegrid[(ix*pny+ny+iy)*pnz+iz];
                }
            }
        }
    }

    if (pme->nnodes_major == 1)
    {
        ny_x = (pme->nnodes_minor == 1 ? ny : pme->pmegrid_ny);

        for (ix = 0; ix < overlap; ix++)
        {
            for (iy = 0; iy < ny_x; iy++)
            {
                for (iz = 0; iz < nz; iz++)
                {
                    pmegrid[(ix*pny+iy)*pnz+iz] +=
                        pmegrid[((nx+ix)*pny+iy)*pnz+iz];
                }
            }
        }
    }
}


static void unwrap_periodic_pmegrid(gmx_pme_t pme, real *pmegrid)
{
    int     nx, ny, nz, pnx, pny, pnz, ny_x, overlap, ix;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    pnx = pme->pmegrid_nx;
    pny = pme->pmegrid_ny;
    pnz = pme->pmegrid_nz;

    overlap = pme->pme_order - 1;

    if (pme->nnodes_major == 1)
    {
        ny_x = (pme->nnodes_minor == 1 ? ny : pme->pmegrid_ny);

        for (ix = 0; ix < overlap; ix++)
        {
            int iy, iz;

            for (iy = 0; iy < ny_x; iy++)
            {
                for (iz = 0; iz < nz; iz++)
                {
                    pmegrid[((nx+ix)*pny+iy)*pnz+iz] =
                        pmegrid[(ix*pny+iy)*pnz+iz];
                }
            }
        }
    }

    if (pme->nnodes_minor == 1)
    {
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
        for (ix = 0; ix < pme->pmegrid_nx; ix++)
        {
            int iy, iz;

            for (iy = 0; iy < overlap; iy++)
            {
                for (iz = 0; iz < nz; iz++)
                {
                    pmegrid[(ix*pny+ny+iy)*pnz+iz] =
                        pmegrid[(ix*pny+iy)*pnz+iz];
                }
            }
        }
    }

    /* Copy periodic overlap in z */
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
    for (ix = 0; ix < pme->pmegrid_nx; ix++)
    {
        int iy, iz;

        for (iy = 0; iy < pme->pmegrid_ny; iy++)
        {
            for (iz = 0; iz < overlap; iz++)
            {
                pmegrid[(ix*pny+iy)*pnz+nz+iz] =
                    pmegrid[(ix*pny+iy)*pnz+iz];
            }
        }
    }
}


/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                            \
    for (ithx = 0; (ithx < order); ithx++)                    \
    {                                                    \
        index_x = (i0+ithx)*pny*pnz;                     \
        valx    = coefficient*thx[ithx];                          \
                                                     \
        for (ithy = 0; (ithy < order); ithy++)                \
        {                                                \
            valxy    = valx*thy[ithy];                   \
            index_xy = index_x+(j0+ithy)*pnz;            \
                                                     \
            for (ithz = 0; (ithz < order); ithz++)            \
            {                                            \
                index_xyz        = index_xy+(k0+ithz);   \
                grid[index_xyz] += valxy*thz[ithz];      \
            }                                            \
        }                                                \
    }


static void spread_coefficients_bsplines_thread(pmegrid_t                    *pmegrid,
                                                pme_atomcomm_t               *atc,
                                                splinedata_t                 *spline,
                                                pme_spline_work_t gmx_unused *work)
{

    /* spread coefficients from home atoms to local grid */
    real          *grid;
    pme_overlap_t *ol;
    int            b, i, nn, n, ithx, ithy, ithz, i0, j0, k0;
    int       *    idxptr;
    int            order, norder, index_x, index_xy, index_xyz;
    real           valx, valxy, coefficient;
    real          *thx, *thy, *thz;
    int            localsize, bndsize;
    int            pnx, pny, pnz, ndatatot;
    int            offx, offy, offz;

#if defined PME_SIMD4_SPREAD_GATHER && !defined PME_SIMD4_UNALIGNED
    real           thz_buffer[GMX_SIMD4_WIDTH*3], *thz_aligned;

    thz_aligned = gmx_simd4_align_r(thz_buffer);
#endif

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx*pny*pnz;
    grid     = pmegrid->grid;
    for (i = 0; i < ndatatot; i++)
    {
        grid[i] = 0;
    }

    order = pmegrid->order;

    for (nn = 0; nn < spline->n; nn++)
    {
        n           = spline->ind[nn];
        coefficient = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = nn*order;

            i0   = idxptr[XX] - offx;
            j0   = idxptr[YY] - offy;
            k0   = idxptr[ZZ] - offz;

            thx = spline->theta[XX] + norder;
            thy = spline->theta[YY] + norder;
            thz = spline->theta[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_SPREAD_SIMD4_ORDER4
#else
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme_simd4.h"
#else
                    DO_BSPLINE(4);
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 5
#include "pme_simd4.h"
#else
                    DO_BSPLINE(5);
#endif
                    break;
                default:
                    DO_BSPLINE(order);
                    break;
            }
        }
    }
}

static void set_grid_alignment(int gmx_unused *pmegrid_nz, int gmx_unused pme_order)
{
#ifdef PME_SIMD4_SPREAD_GATHER
    if (pme_order == 5
#ifndef PME_SIMD4_UNALIGNED
        || pme_order == 4
#endif
        )
    {
        /* Round nz up to a multiple of 4 to ensure alignment */
        *pmegrid_nz = ((*pmegrid_nz + 3) & ~3);
    }
#endif
}

static void set_gridsize_alignment(int gmx_unused *gridsize, int gmx_unused pme_order)
{
#ifdef PME_SIMD4_SPREAD_GATHER
#ifndef PME_SIMD4_UNALIGNED
    if (pme_order == 4)
    {
        /* Add extra elements to ensured aligned operations do not go
         * beyond the allocated grid size.
         * Note that for pme_order=5, the pme grid z-size alignment
         * ensures that we will not go beyond the grid size.
         */
        *gridsize += 4;
    }
#endif
#endif
}

static void pmegrid_init(pmegrid_t *grid,
                         int cx, int cy, int cz,
                         int x0, int y0, int z0,
                         int x1, int y1, int z1,
                         gmx_bool set_alignment,
                         int pme_order,
                         real *ptr)
{
    int nz, gridsize;

    grid->ci[XX]     = cx;
    grid->ci[YY]     = cy;
    grid->ci[ZZ]     = cz;
    grid->offset[XX] = x0;
    grid->offset[YY] = y0;
    grid->offset[ZZ] = z0;
    grid->n[XX]      = x1 - x0 + pme_order - 1;
    grid->n[YY]      = y1 - y0 + pme_order - 1;
    grid->n[ZZ]      = z1 - z0 + pme_order - 1;
    copy_ivec(grid->n, grid->s);

    nz = grid->s[ZZ];
    set_grid_alignment(&nz, pme_order);
    if (set_alignment)
    {
        grid->s[ZZ] = nz;
    }
    else if (nz != grid->s[ZZ])
    {
        gmx_incons("pmegrid_init call with an unaligned z size");
    }

    grid->order = pme_order;
    if (ptr == NULL)
    {
        gridsize = grid->s[XX]*grid->s[YY]*grid->s[ZZ];
        set_gridsize_alignment(&gridsize, pme_order);
        snew_aligned(grid->grid, gridsize, SIMD4_ALIGNMENT);
    }
    else
    {
        grid->grid = ptr;
    }
}

static int div_round_up(int enumerator, int denominator)
{
    return (enumerator + denominator - 1)/denominator;
}

static void make_subgrid_division(const ivec n, int ovl, int nthread,
                                  ivec nsub)
{
    int gsize_opt, gsize;
    int nsx, nsy, nsz;
    char *env;

    gsize_opt = -1;
    for (nsx = 1; nsx <= nthread; nsx++)
    {
        if (nthread % nsx == 0)
        {
            for (nsy = 1; nsy <= nthread; nsy++)
            {
                if (nsx*nsy <= nthread && nthread % (nsx*nsy) == 0)
                {
                    nsz = nthread/(nsx*nsy);

                    /* Determine the number of grid points per thread */
                    gsize =
                        (div_round_up(n[XX], nsx) + ovl)*
                        (div_round_up(n[YY], nsy) + ovl)*
                        (div_round_up(n[ZZ], nsz) + ovl);

                    /* Minimize the number of grids points per thread
                     * and, secondarily, the number of cuts in minor dimensions.
                     */
                    if (gsize_opt == -1 ||
                        gsize < gsize_opt ||
                        (gsize == gsize_opt &&
                         (nsz < nsub[ZZ] || (nsz == nsub[ZZ] && nsy < nsub[YY]))))
                    {
                        nsub[XX]  = nsx;
                        nsub[YY]  = nsy;
                        nsub[ZZ]  = nsz;
                        gsize_opt = gsize;
                    }
                }
            }
        }
    }

    env = getenv("GMX_PME_THREAD_DIVISION");
    if (env != NULL)
    {
        sscanf(env, "%d %d %d", &nsub[XX], &nsub[YY], &nsub[ZZ]);
    }

    if (nsub[XX]*nsub[YY]*nsub[ZZ] != nthread)
    {
        gmx_fatal(FARGS, "PME grid thread division (%d x %d x %d) does not match the total number of threads (%d)", nsub[XX], nsub[YY], nsub[ZZ], nthread);
    }
}

static void pmegrids_init(pmegrids_t *grids,
                          int nx, int ny, int nz, int nz_base,
                          int pme_order,
                          gmx_bool bUseThreads,
                          int nthread,
                          int overlap_x,
                          int overlap_y)
{
    ivec n, n_base, g0, g1;
    int t, x, y, z, d, i, tfac;
    int max_comm_lines = -1;

    n[XX] = nx - (pme_order - 1);
    n[YY] = ny - (pme_order - 1);
    n[ZZ] = nz - (pme_order - 1);

    copy_ivec(n, n_base);
    n_base[ZZ] = nz_base;

    pmegrid_init(&grids->grid, 0, 0, 0, 0, 0, 0, n[XX], n[YY], n[ZZ], FALSE, pme_order,
                 NULL);

    grids->nthread = nthread;

    make_subgrid_division(n_base, pme_order-1, grids->nthread, grids->nc);

    if (bUseThreads)
    {
        ivec nst;
        int gridsize;

        for (d = 0; d < DIM; d++)
        {
            nst[d] = div_round_up(n[d], grids->nc[d]) + pme_order - 1;
        }
        set_grid_alignment(&nst[ZZ], pme_order);

        if (debug)
        {
            fprintf(debug, "pmegrid thread local division: %d x %d x %d\n",
                    grids->nc[XX], grids->nc[YY], grids->nc[ZZ]);
            fprintf(debug, "pmegrid %d %d %d max thread pmegrid %d %d %d\n",
                    nx, ny, nz,
                    nst[XX], nst[YY], nst[ZZ]);
        }

        snew(grids->grid_th, grids->nthread);
        t        = 0;
        gridsize = nst[XX]*nst[YY]*nst[ZZ];
        set_gridsize_alignment(&gridsize, pme_order);
        snew_aligned(grids->grid_all,
                     grids->nthread*gridsize+(grids->nthread+1)*GMX_CACHE_SEP,
                     SIMD4_ALIGNMENT);

        for (x = 0; x < grids->nc[XX]; x++)
        {
            for (y = 0; y < grids->nc[YY]; y++)
            {
                for (z = 0; z < grids->nc[ZZ]; z++)
                {
                    pmegrid_init(&grids->grid_th[t],
                                 x, y, z,
                                 (n[XX]*(x  ))/grids->nc[XX],
                                 (n[YY]*(y  ))/grids->nc[YY],
                                 (n[ZZ]*(z  ))/grids->nc[ZZ],
                                 (n[XX]*(x+1))/grids->nc[XX],
                                 (n[YY]*(y+1))/grids->nc[YY],
                                 (n[ZZ]*(z+1))/grids->nc[ZZ],
                                 TRUE,
                                 pme_order,
                                 grids->grid_all+GMX_CACHE_SEP+t*(gridsize+GMX_CACHE_SEP));
                    t++;
                }
            }
        }
    }
    else
    {
        grids->grid_th = NULL;
    }

    snew(grids->g2t, DIM);
    tfac = 1;
    for (d = DIM-1; d >= 0; d--)
    {
        snew(grids->g2t[d], n[d]);
        t = 0;
        for (i = 0; i < n[d]; i++)
        {
            /* The second check should match the parameters
             * of the pmegrid_init call above.
             */
            while (t + 1 < grids->nc[d] && i >= (n[d]*(t+1))/grids->nc[d])
            {
                t++;
            }
            grids->g2t[d][i] = t*tfac;
        }

        tfac *= grids->nc[d];

        switch (d)
        {
            case XX: max_comm_lines = overlap_x;     break;
            case YY: max_comm_lines = overlap_y;     break;
            case ZZ: max_comm_lines = pme_order - 1; break;
        }
        grids->nthread_comm[d] = 0;
        while ((n[d]*grids->nthread_comm[d])/grids->nc[d] < max_comm_lines &&
               grids->nthread_comm[d] < grids->nc[d])
        {
            grids->nthread_comm[d]++;
        }
        if (debug != NULL)
        {
            fprintf(debug, "pmegrid thread grid communication range in %c: %d\n",
                    'x'+d, grids->nthread_comm[d]);
        }
        /* It should be possible to make grids->nthread_comm[d]==grids->nc[d]
         * work, but this is not a problematic restriction.
         */
        if (grids->nc[d] > 1 && grids->nthread_comm[d] > grids->nc[d])
        {
            gmx_fatal(FARGS, "Too many threads for PME (%d) compared to the number of grid lines, reduce the number of threads doing PME", grids->nthread);
        }
    }
}


static void pmegrids_destroy(pmegrids_t *grids)
{
    int t;

    if (grids->grid.grid != NULL)
    {
        sfree(grids->grid.grid);

        if (grids->nthread > 0)
        {
            for (t = 0; t < grids->nthread; t++)
            {
                sfree(grids->grid_th[t].grid);
            }
            sfree(grids->grid_th);
        }
    }
}


static void realloc_work(pme_work_t *work, int nkx)
{
    int simd_width;

    if (nkx > work->nalloc)
    {
        work->nalloc = nkx;
        srenew(work->mhx, work->nalloc);
        srenew(work->mhy, work->nalloc);
        srenew(work->mhz, work->nalloc);
        srenew(work->m2, work->nalloc);
        /* Allocate an aligned pointer for SIMD operations, including extra
         * elements at the end for padding.
         */
#ifdef PME_SIMD_SOLVE
        simd_width = GMX_SIMD_REAL_WIDTH;
#else
        /* We can use any alignment, apart from 0, so we use 4 */
        simd_width = 4;
#endif
        sfree_aligned(work->denom);
        sfree_aligned(work->tmp1);
        sfree_aligned(work->tmp2);
        sfree_aligned(work->eterm);
        snew_aligned(work->denom, work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->tmp1,  work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->tmp2,  work->nalloc+simd_width, simd_width*sizeof(real));
        snew_aligned(work->eterm, work->nalloc+simd_width, simd_width*sizeof(real));
        srenew(work->m2inv, work->nalloc);
    }
}


static void free_work(pme_work_t *work)
{
    sfree(work->mhx);
    sfree(work->mhy);
    sfree(work->mhz);
    sfree(work->m2);
    sfree_aligned(work->denom);
    sfree_aligned(work->tmp1);
    sfree_aligned(work->tmp2);
    sfree_aligned(work->eterm);
    sfree(work->m2inv);
}


#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
gmx_inline static void calc_exponentials_q(int gmx_unused start, int end, real f, real *d_aligned, real *r_aligned, real *e_aligned)
{
    {
        const gmx_simd_real_t two = gmx_simd_set1_r(2.0);
        gmx_simd_real_t f_simd;
        gmx_simd_real_t lu;
        gmx_simd_real_t tmp_d1, d_inv, tmp_r, tmp_e;
        int kx;
        f_simd = gmx_simd_set1_r(f);
        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        for (kx = 0; kx < end; kx += GMX_SIMD_REAL_WIDTH)
        {
            tmp_d1   = gmx_simd_load_r(d_aligned+kx);
            d_inv    = gmx_simd_inv_r(tmp_d1);
            tmp_r    = gmx_simd_load_r(r_aligned+kx);
            tmp_r    = gmx_simd_exp_r(tmp_r);
            tmp_e    = gmx_simd_mul_r(f_simd, d_inv);
            tmp_e    = gmx_simd_mul_r(tmp_e, tmp_r);
            gmx_simd_store_r(e_aligned+kx, tmp_e);
        }
    }
}
#else
gmx_inline static void calc_exponentials_q(int start, int end, real f, real *d, real *r, real *e)
{
    int kx;
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0/d[kx];
    }
    for (kx = start; kx < end; kx++)
    {
        r[kx] = exp(r[kx]);
    }
    for (kx = start; kx < end; kx++)
    {
        e[kx] = f*r[kx]*d[kx];
    }
}
#endif

#if defined PME_SIMD_SOLVE
/* Calculate exponentials through SIMD */
gmx_inline static void calc_exponentials_lj(int gmx_unused start, int end, real *r_aligned, real *factor_aligned, real *d_aligned)
{
    gmx_simd_real_t tmp_r, tmp_d, tmp_fac, d_inv, tmp_mk;
    const gmx_simd_real_t sqr_PI = gmx_simd_sqrt_r(gmx_simd_set1_r(M_PI));
    int kx;
    for (kx = 0; kx < end; kx += GMX_SIMD_REAL_WIDTH)
    {
        /* We only need to calculate from start. But since start is 0 or 1
         * and we want to use aligned loads/stores, we always start from 0.
         */
        tmp_d = gmx_simd_load_r(d_aligned+kx);
        d_inv = gmx_simd_inv_r(tmp_d);
        gmx_simd_store_r(d_aligned+kx, d_inv);
        tmp_r = gmx_simd_load_r(r_aligned+kx);
        tmp_r = gmx_simd_exp_r(tmp_r);
        gmx_simd_store_r(r_aligned+kx, tmp_r);
        tmp_mk  = gmx_simd_load_r(factor_aligned+kx);
        tmp_fac = gmx_simd_mul_r(sqr_PI, gmx_simd_mul_r(tmp_mk, gmx_simd_erfc_r(tmp_mk)));
        gmx_simd_store_r(factor_aligned+kx, tmp_fac);
    }
}
#else
gmx_inline static void calc_exponentials_lj(int start, int end, real *r, real *tmp2, real *d)
{
    int kx;
    real mk;
    for (kx = start; kx < end; kx++)
    {
        d[kx] = 1.0/d[kx];
    }

    for (kx = start; kx < end; kx++)
    {
        r[kx] = exp(r[kx]);
    }

    for (kx = start; kx < end; kx++)
    {
        mk       = tmp2[kx];
        tmp2[kx] = sqrt(M_PI)*mk*gmx_erfc(mk);
    }
}
#endif

static int solve_pme_yzx(gmx_pme_t pme, t_complex *grid,
                         real ewaldcoeff, real vol,
                         gmx_bool bEnerVir,
                         int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    t_complex *p0;
    int     kx, ky, kz, maxkx, maxky, maxkz;
    int     nx, ny, nz, iyz0, iyz1, iyz, iy, iz, kxstart, kxend;
    real    mx, my, mz;
    real    factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real    ets2, struct2, vfactor, ets2vf;
    real    d1, d2, energy = 0;
    real    by, bz;
    real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real    rxx, ryx, ryy, rzx, rzy, rzz;
    pme_work_t *work;
    real    *mhx, *mhy, *mhz, *m2, *denom, *tmp1, *eterm, *m2inv;
    real    mhxk, mhyk, mhzk, m2k;
    real    corner_fac;
    ivec    complex_order;
    ivec    local_ndata, local_offset, local_size;
    real    elfac;

    elfac = ONE_4PI_EPS0/pme->epsilon_r;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_QA],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;
    maxkz = nz/2+1;

    work  = &pme->work[thread];
    mhx   = work->mhx;
    mhy   = work->mhy;
    mhz   = work->mhz;
    m2    = work->m2;
    denom = work->denom;
    tmp1  = work->tmp1;
    eterm = work->eterm;
    m2inv = work->m2inv;

    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz/local_ndata[ZZ];
        iz = iyz - iy*local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = M_PI*vol*pme->bsp_mod[YY][ky];

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme->bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        p0 = grid + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];

        /* We should skip the k-space point (0,0,0) */
        /* Note that since here x is the minor index, local_offset[XX]=0 */
        if (local_offset[XX] > 0 || ky > 0 || kz > 0)
        {
            kxstart = local_offset[XX];
        }
        else
        {
            kxstart = local_offset[XX] + 1;
            p0++;
        }
        kxend = local_offset[XX] + local_ndata[XX];

        if (bEnerVir)
        {
            /* More expensive inner loop, especially because of the storage
             * of the mh elements in array's.
             * Because x is the minor grid index, all mh elements
             * depend on kx for triclinic unit cells.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2inv[kx] = 1.0/m2[kx];
            }

            calc_exponentials_q(kxstart, kxend, elfac, denom, tmp1, eterm);

            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1      = p0->re;
                d2      = p0->im;

                p0->re  = d1*eterm[kx];
                p0->im  = d2*eterm[kx];

                struct2 = 2.0*(d1*d1+d2*d2);

                tmp1[kx] = eterm[kx]*struct2;
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2     = corner_fac*tmp1[kx];
                vfactor  = (factor*m2[kx] + 1.0)*2.0*m2inv[kx];
                energy  += ets2;

                ets2vf   = ets2*vfactor;
                virxx   += ets2vf*mhx[kx]*mhx[kx] - ets2;
                virxy   += ets2vf*mhx[kx]*mhy[kx];
                virxz   += ets2vf*mhx[kx]*mhz[kx];
                viryy   += ets2vf*mhy[kx]*mhy[kx] - ets2;
                viryz   += ets2vf*mhy[kx]*mhz[kx];
                virzz   += ets2vf*mhz[kx]*mhz[kx] - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             * In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */

            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
            }

            calc_exponentials_q(kxstart, kxend, elfac, denom, tmp1, eterm);

            for (kx = kxstart; kx < kxend; kx++, p0++)
            {
                d1      = p0->re;
                d2      = p0->im;

                p0->re  = d1*eterm[kx];
                p0->im  = d2*eterm[kx];
            }
        }
    }

    if (bEnerVir)
    {
        /* Update virial with local values.
         * The virial is symmetric by definition.
         * this virial seems ok for isotropic scaling, but I'm
         * experiencing problems on semiisotropic membranes.
         * IS THAT COMMENT STILL VALID??? (DvdS, 2001/02/07).
         */
        work->vir_q[XX][XX] = 0.25*virxx;
        work->vir_q[YY][YY] = 0.25*viryy;
        work->vir_q[ZZ][ZZ] = 0.25*virzz;
        work->vir_q[XX][YY] = work->vir_q[YY][XX] = 0.25*virxy;
        work->vir_q[XX][ZZ] = work->vir_q[ZZ][XX] = 0.25*virxz;
        work->vir_q[YY][ZZ] = work->vir_q[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        work->energy_q = 0.5*energy;
    }

    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}

static int solve_pme_lj_yzx(gmx_pme_t pme, t_complex **grid, gmx_bool bLB,
                            real ewaldcoeff, real vol,
                            gmx_bool bEnerVir, int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    int     ig, gcount;
    int     kx, ky, kz, maxkx, maxky, maxkz;
    int     nx, ny, nz, iy, iyz0, iyz1, iyz, iz, kxstart, kxend;
    real    mx, my, mz;
    real    factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real    ets2, ets2vf;
    real    eterm, vterm, d1, d2, energy = 0;
    real    by, bz;
    real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    real    rxx, ryx, ryy, rzx, rzy, rzz;
    real    *mhx, *mhy, *mhz, *m2, *denom, *tmp1, *tmp2;
    real    mhxk, mhyk, mhzk, m2k;
    real    mk;
    pme_work_t *work;
    real    corner_fac;
    ivec    complex_order;
    ivec    local_ndata, local_offset, local_size;
    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_C6A],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;
    maxkz = nz/2+1;

    work  = &pme->work[thread];
    mhx   = work->mhx;
    mhy   = work->mhy;
    mhz   = work->mhz;
    m2    = work->m2;
    denom = work->denom;
    tmp1  = work->tmp1;
    tmp2  = work->tmp2;

    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    for (iyz = iyz0; iyz < iyz1; iyz++)
    {
        iy = iyz/local_ndata[ZZ];
        iz = iyz - iy*local_ndata[ZZ];

        ky = iy + local_offset[YY];

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        by = 3.0*vol*pme->bsp_mod[YY][ky]
            / (M_PI*sqrt(M_PI)*ewaldcoeff*ewaldcoeff*ewaldcoeff);

        kz = iz + local_offset[ZZ];

        mz = kz;

        bz = pme->bsp_mod[ZZ][kz];

        /* 0.5 correction for corner points */
        corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        kxstart = local_offset[XX];
        kxend   = local_offset[XX] + local_ndata[XX];
        if (bEnerVir)
        {
            /* More expensive inner loop, especially because of the
             * storage of the mh elements in array's.  Because x is the
             * minor grid index, all mh elements depend on kx for
             * triclinic unit cells.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                mhx[kx]   = mhxk;
                mhy[kx]   = mhyk;
                mhz[kx]   = mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            calc_exponentials_lj(kxstart, kxend, tmp1, tmp2, denom);

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k   = factor*m2[kx];
                eterm = -((1.0 - 2.0*m2k)*tmp1[kx]
                          + 2.0*m2k*tmp2[kx]);
                vterm    = 3.0*(-tmp1[kx] + tmp2[kx]);
                tmp1[kx] = eterm*denom[kx];
                tmp2[kx] = vterm*denom[kx];
            }

            if (!bLB)
            {
                t_complex *p0;
                real       struct2;

                p0 = grid[0] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;

                    eterm   = tmp1[kx];
                    vterm   = tmp2[kx];
                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;

                    struct2 = 2.0*(d1*d1+d2*d2);

                    tmp1[kx] = eterm*struct2;
                    tmp2[kx] = vterm*struct2;
                }
            }
            else
            {
                real *struct2 = denom;
                real  str2;

                for (kx = kxstart; kx < kxend; kx++)
                {
                    struct2[kx] = 0.0;
                }
                /* Due to symmetry we only need to calculate 4 of the 7 terms */
                for (ig = 0; ig <= 3; ++ig)
                {
                    t_complex *p0, *p1;
                    real       scale;

                    p0    = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    p1    = grid[6-ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    scale = 2.0*lb_scale_factor_symm[ig];
                    for (kx = kxstart; kx < kxend; ++kx, ++p0, ++p1)
                    {
                        struct2[kx] += scale*(p0->re*p1->re + p0->im*p1->im);
                    }

                }
                for (ig = 0; ig <= 6; ++ig)
                {
                    t_complex *p0;

                    p0 = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                    for (kx = kxstart; kx < kxend; kx++, p0++)
                    {
                        d1     = p0->re;
                        d2     = p0->im;

                        eterm  = tmp1[kx];
                        p0->re = d1*eterm;
                        p0->im = d2*eterm;
                    }
                }
                for (kx = kxstart; kx < kxend; kx++)
                {
                    eterm    = tmp1[kx];
                    vterm    = tmp2[kx];
                    str2     = struct2[kx];
                    tmp1[kx] = eterm*str2;
                    tmp2[kx] = vterm*str2;
                }
            }

            for (kx = kxstart; kx < kxend; kx++)
            {
                ets2     = corner_fac*tmp1[kx];
                vterm    = 2.0*factor*tmp2[kx];
                energy  += ets2;
                ets2vf   = corner_fac*vterm;
                virxx   += ets2vf*mhx[kx]*mhx[kx] - ets2;
                virxy   += ets2vf*mhx[kx]*mhy[kx];
                virxz   += ets2vf*mhx[kx]*mhz[kx];
                viryy   += ets2vf*mhy[kx]*mhy[kx] - ets2;
                viryz   += ets2vf*mhy[kx]*mhz[kx];
                virzz   += ets2vf*mhz[kx]*mhz[kx] - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             *  In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */

            for (kx = kxstart; kx < maxkx; kx++)
            {
                mx = kx;

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            for (kx = maxkx; kx < kxend; kx++)
            {
                mx = (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                m2[kx]    = m2k;
                denom[kx] = bz*by*pme->bsp_mod[XX][kx];
                tmp1[kx]  = -factor*m2k;
                tmp2[kx]  = sqrt(factor*m2k);
            }

            calc_exponentials_lj(kxstart, kxend, tmp1, tmp2, denom);

            for (kx = kxstart; kx < kxend; kx++)
            {
                m2k    = factor*m2[kx];
                eterm  = -((1.0 - 2.0*m2k)*tmp1[kx]
                           + 2.0*m2k*tmp2[kx]);
                tmp1[kx] = eterm*denom[kx];
            }
            gcount = (bLB ? 7 : 1);
            for (ig = 0; ig < gcount; ++ig)
            {
                t_complex *p0;

                p0 = grid[ig] + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
                for (kx = kxstart; kx < kxend; kx++, p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;

                    eterm   = tmp1[kx];

                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;
                }
            }
        }
    }
    if (bEnerVir)
    {
        work->vir_lj[XX][XX] = 0.25*virxx;
        work->vir_lj[YY][YY] = 0.25*viryy;
        work->vir_lj[ZZ][ZZ] = 0.25*virzz;
        work->vir_lj[XX][YY] = work->vir_lj[YY][XX] = 0.25*virxy;
        work->vir_lj[XX][ZZ] = work->vir_lj[ZZ][XX] = 0.25*virxz;
        work->vir_lj[YY][ZZ] = work->vir_lj[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        work->energy_lj = 0.5*energy;
    }
    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}

static void get_pme_ener_vir_q(const gmx_pme_t pme, int nthread,
                               real *mesh_energy, matrix vir)
{
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    int thread;

    *mesh_energy = pme->work[0].energy_q;
    copy_mat(pme->work[0].vir_q, vir);

    for (thread = 1; thread < nthread; thread++)
    {
        *mesh_energy += pme->work[thread].energy_q;
        m_add(vir, pme->work[thread].vir_q, vir);
    }
}

static void get_pme_ener_vir_lj(const gmx_pme_t pme, int nthread,
                                real *mesh_energy, matrix vir)
{
    /* This function sums output over threads and should therefore
     * only be called after thread synchronization.
     */
    int thread;

    *mesh_energy = pme->work[0].energy_lj;
    copy_mat(pme->work[0].vir_lj, vir);

    for (thread = 1; thread < nthread; thread++)
    {
        *mesh_energy += pme->work[thread].energy_lj;
        m_add(vir, pme->work[thread].vir_lj, vir);
    }
}


#define DO_FSPLINE(order)                      \
    for (ithx = 0; (ithx < order); ithx++)              \
    {                                              \
        index_x = (i0+ithx)*pny*pnz;               \
        tx      = thx[ithx];                       \
        dx      = dthx[ithx];                      \
                                               \
        for (ithy = 0; (ithy < order); ithy++)          \
        {                                          \
            index_xy = index_x+(j0+ithy)*pnz;      \
            ty       = thy[ithy];                  \
            dy       = dthy[ithy];                 \
            fxy1     = fz1 = 0;                    \
                                               \
            for (ithz = 0; (ithz < order); ithz++)      \
            {                                      \
                gval  = grid[index_xy+(k0+ithz)];  \
                fxy1 += thz[ithz]*gval;            \
                fz1  += dthz[ithz]*gval;           \
            }                                      \
            fx += dx*ty*fxy1;                      \
            fy += tx*dy*fxy1;                      \
            fz += tx*ty*fz1;                       \
        }                                          \
    }


static void gather_f_bsplines(gmx_pme_t pme, real *grid,
                              gmx_bool bClearF, pme_atomcomm_t *atc,
                              splinedata_t *spline,
                              real scale)
{
    /* sum forces for local particles */
    int     nn, n, ithx, ithy, ithz, i0, j0, k0;
    int     index_x, index_xy;
    int     nx, ny, nz, pnx, pny, pnz;
    int *   idxptr;
    real    tx, ty, dx, dy, coefficient;
    real    fx, fy, fz, gval;
    real    fxy1, fz1;
    real    *thx, *thy, *thz, *dthx, *dthy, *dthz;
    int     norder;
    real    rxx, ryx, ryy, rzx, rzy, rzz;
    int     order;

    pme_spline_work_t *work;

#if defined PME_SIMD4_SPREAD_GATHER && !defined PME_SIMD4_UNALIGNED
    real           thz_buffer[GMX_SIMD4_WIDTH*3],  *thz_aligned;
    real           dthz_buffer[GMX_SIMD4_WIDTH*3], *dthz_aligned;

    thz_aligned  = gmx_simd4_align_r(thz_buffer);
    dthz_aligned = gmx_simd4_align_r(dthz_buffer);
#endif

    work = pme->spline_work;

    order = pme->pme_order;
    thx   = spline->theta[XX];
    thy   = spline->theta[YY];
    thz   = spline->theta[ZZ];
    dthx  = spline->dtheta[XX];
    dthy  = spline->dtheta[YY];
    dthz  = spline->dtheta[ZZ];
    nx    = pme->nkx;
    ny    = pme->nky;
    nz    = pme->nkz;
    pnx   = pme->pmegrid_nx;
    pny   = pme->pmegrid_ny;
    pnz   = pme->pmegrid_nz;

    rxx   = pme->recipbox[XX][XX];
    ryx   = pme->recipbox[YY][XX];
    ryy   = pme->recipbox[YY][YY];
    rzx   = pme->recipbox[ZZ][XX];
    rzy   = pme->recipbox[ZZ][YY];
    rzz   = pme->recipbox[ZZ][ZZ];

    for (nn = 0; nn < spline->n; nn++)
    {
        n           = spline->ind[nn];
        coefficient = scale*atc->coefficient[n];

        if (bClearF)
        {
            atc->f[n][XX] = 0;
            atc->f[n][YY] = 0;
            atc->f[n][ZZ] = 0;
        }
        if (coefficient != 0)
        {
            fx     = 0;
            fy     = 0;
            fz     = 0;
            idxptr = atc->idx[n];
            norder = nn*order;

            i0   = idxptr[XX];
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];

            /* Pointer arithmetic alert, next six statements */
            thx  = spline->theta[XX] + norder;
            thy  = spline->theta[YY] + norder;
            thz  = spline->theta[ZZ] + norder;
            dthx = spline->dtheta[XX] + norder;
            dthy = spline->dtheta[YY] + norder;
            dthz = spline->dtheta[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_GATHER_F_SIMD4_ORDER4
#else
#define PME_GATHER_F_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme_simd4.h"
#else
                    DO_FSPLINE(4);
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#define PME_GATHER_F_SIMD4_ALIGNED
#define PME_ORDER 5
#include "pme_simd4.h"
#else
                    DO_FSPLINE(5);
#endif
                    break;
                default:
                    DO_FSPLINE(order);
                    break;
            }

            atc->f[n][XX] += -coefficient*( fx*nx*rxx );
            atc->f[n][YY] += -coefficient*( fx*nx*ryx + fy*ny*ryy );
            atc->f[n][ZZ] += -coefficient*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
        }
    }
    /* Since the energy and not forces are interpolated
     * the net force might not be exactly zero.
     * This can be solved by also interpolating F, but
     * that comes at a cost.
     * A better hack is to remove the net force every
     * step, but that must be done at a higher level
     * since this routine doesn't see all atoms if running
     * in parallel. Don't know how important it is?  EL 990726
     */
}


static real gather_energy_bsplines(gmx_pme_t pme, real *grid,
                                   pme_atomcomm_t *atc)
{
    splinedata_t *spline;
    int     n, ithx, ithy, ithz, i0, j0, k0;
    int     index_x, index_xy;
    int *   idxptr;
    real    energy, pot, tx, ty, coefficient, gval;
    real    *thx, *thy, *thz;
    int     norder;
    int     order;

    spline = &atc->spline[0];

    order = pme->pme_order;

    energy = 0;
    for (n = 0; (n < atc->n); n++)
    {
        coefficient      = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = n*order;

            i0   = idxptr[XX];
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];

            /* Pointer arithmetic alert, next three statements */
            thx  = spline->theta[XX] + norder;
            thy  = spline->theta[YY] + norder;
            thz  = spline->theta[ZZ] + norder;

            pot = 0;
            for (ithx = 0; (ithx < order); ithx++)
            {
                index_x = (i0+ithx)*pme->pmegrid_ny*pme->pmegrid_nz;
                tx      = thx[ithx];

                for (ithy = 0; (ithy < order); ithy++)
                {
                    index_xy = index_x+(j0+ithy)*pme->pmegrid_nz;
                    ty       = thy[ithy];

                    for (ithz = 0; (ithz < order); ithz++)
                    {
                        gval  = grid[index_xy+(k0+ithz)];
                        pot  += tx*ty*thz[ithz]*gval;
                    }

                }
            }

            energy += pot*coefficient;
        }
    }

    return energy;
}

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE(order)                     \
    {                                              \
        int j, k, l;                                 \
        real dr, div;                               \
        real data[PME_ORDER_MAX];                  \
        real ddata[PME_ORDER_MAX];                 \
                                               \
        for (j = 0; (j < DIM); j++)                     \
        {                                          \
            dr  = xptr[j];                         \
                                               \
            /* dr is relative offset from lower cell limit */ \
            data[order-1] = 0;                     \
            data[1]       = dr;                          \
            data[0]       = 1 - dr;                      \
                                               \
            for (k = 3; (k < order); k++)               \
            {                                      \
                div       = 1.0/(k - 1.0);               \
                data[k-1] = div*dr*data[k-2];      \
                for (l = 1; (l < (k-1)); l++)           \
                {                                  \
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)* \
                                       data[k-l-1]);                \
                }                                  \
                data[0] = div*(1-dr)*data[0];      \
            }                                      \
            /* differentiate */                    \
            ddata[0] = -data[0];                   \
            for (k = 1; (k < order); k++)               \
            {                                      \
                ddata[k] = data[k-1] - data[k];    \
            }                                      \
                                               \
            div           = 1.0/(order - 1);                 \
            data[order-1] = div*dr*data[order-2];  \
            for (l = 1; (l < (order-1)); l++)           \
            {                                      \
                data[order-l-1] = div*((dr+l)*data[order-l-2]+    \
                                       (order-l-dr)*data[order-l-1]); \
            }                                      \
            data[0] = div*(1 - dr)*data[0];        \
                                               \
            for (k = 0; k < order; k++)                 \
            {                                      \
                theta[j][i*order+k]  = data[k];    \
                dtheta[j][i*order+k] = ddata[k];   \
            }                                      \
        }                                          \
    }

void make_bsplines(splinevec theta, splinevec dtheta, int order,
                   rvec fractx[], int nr, int ind[], real coefficient[],
                   gmx_bool bDoSplines)
{
    /* construct splines for local atoms */
    int  i, ii;
    real *xptr;

    for (i = 0; i < nr; i++)
    {
        /* With free energy we do not use the coefficient check.
         * In most cases this will be more efficient than calling make_bsplines
         * twice, since usually more than half the particles have non-zero coefficients.
         */
        ii = ind[i];
        if (bDoSplines || coefficient[ii] != 0.0)
        {
            xptr = fractx[ii];
            switch (order)
            {
                case 4:  CALC_SPLINE(4);     break;
                case 5:  CALC_SPLINE(5);     break;
                default: CALC_SPLINE(order); break;
            }
        }
    }
}


void make_dft_mod(real *mod, real *data, int ndata)
{
    int i, j;
    real sc, ss, arg;

    for (i = 0; i < ndata; i++)
    {
        sc = ss = 0;
        for (j = 0; j < ndata; j++)
        {
            arg = (2.0*M_PI*i*j)/ndata;
            sc += data[j]*cos(arg);
            ss += data[j]*sin(arg);
        }
        mod[i] = sc*sc+ss*ss;
    }
    for (i = 0; i < ndata; i++)
    {
        if (mod[i] < 1e-7)
        {
            mod[i] = (mod[i-1]+mod[i+1])*0.5;
        }
    }
}


static void make_bspline_moduli(splinevec bsp_mod,
                                int nx, int ny, int nz, int order)
{
    int nmax = max(nx, max(ny, nz));
    real *data, *ddata, *bsp_data;
    int i, k, l;
    real div;

    snew(data, order);
    snew(ddata, order);
    snew(bsp_data, nmax);

    data[order-1] = 0;
    data[1]       = 0;
    data[0]       = 1;

    for (k = 3; k < order; k++)
    {
        div       = 1.0/(k-1.0);
        data[k-1] = 0;
        for (l = 1; l < (k-1); l++)
        {
            data[k-l-1] = div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
        }
        data[0] = div*data[0];
    }
    /* differentiate */
    ddata[0] = -data[0];
    for (k = 1; k < order; k++)
    {
        ddata[k] = data[k-1]-data[k];
    }
    div           = 1.0/(order-1);
    data[order-1] = 0;
    for (l = 1; l < (order-1); l++)
    {
        data[order-l-1] = div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
    }
    data[0] = div*data[0];

    for (i = 0; i < nmax; i++)
    {
        bsp_data[i] = 0;
    }
    for (i = 1; i <= order; i++)
    {
        bsp_data[i] = data[i-1];
    }

    make_dft_mod(bsp_mod[XX], bsp_data, nx);
    make_dft_mod(bsp_mod[YY], bsp_data, ny);
    make_dft_mod(bsp_mod[ZZ], bsp_data, nz);

    sfree(data);
    sfree(ddata);
    sfree(bsp_data);
}


/* Return the P3M optimal influence function */
static double do_p3m_influence(double z, int order)
{
    double z2, z4;

    z2 = z*z;
    z4 = z2*z2;

    /* The formula and most constants can be found in:
     * Ballenegger et al., JCTC 8, 936 (2012)
     */
    switch (order)
    {
        case 2:
            return 1.0 - 2.0*z2/3.0;
            break;
        case 3:
            return 1.0 - z2 + 2.0*z4/15.0;
            break;
        case 4:
            return 1.0 - 4.0*z2/3.0 + 2.0*z4/5.0 + 4.0*z2*z4/315.0;
            break;
        case 5:
            return 1.0 - 5.0*z2/3.0 + 7.0*z4/9.0 - 17.0*z2*z4/189.0 + 2.0*z4*z4/2835.0;
            break;
        case 6:
            return 1.0 - 2.0*z2 + 19.0*z4/15.0 - 256.0*z2*z4/945.0 + 62.0*z4*z4/4725.0 + 4.0*z2*z4*z4/155925.0;
            break;
        case 7:
            return 1.0 - 7.0*z2/3.0 + 28.0*z4/15.0 - 16.0*z2*z4/27.0 + 26.0*z4*z4/405.0 - 2.0*z2*z4*z4/1485.0 + 4.0*z4*z4*z4/6081075.0;
        case 8:
            return 1.0 - 8.0*z2/3.0 + 116.0*z4/45.0 - 344.0*z2*z4/315.0 + 914.0*z4*z4/4725.0 - 248.0*z4*z4*z2/22275.0 + 21844.0*z4*z4*z4/212837625.0 - 8.0*z4*z4*z4*z2/638512875.0;
            break;
    }

    return 0.0;
}

/* Calculate the P3M B-spline moduli for one dimension */
static void make_p3m_bspline_moduli_dim(real *bsp_mod, int n, int order)
{
    double zarg, zai, sinzai, infl;
    int    maxk, i;

    if (order > 8)
    {
        gmx_fatal(FARGS, "The current P3M code only supports orders up to 8");
    }

    zarg = M_PI/n;

    maxk = (n + 1)/2;

    for (i = -maxk; i < 0; i++)
    {
        zai          = zarg*i;
        sinzai       = sin(zai);
        infl         = do_p3m_influence(sinzai, order);
        bsp_mod[n+i] = infl*infl*pow(sinzai/zai, -2.0*order);
    }
    bsp_mod[0] = 1.0;
    for (i = 1; i < maxk; i++)
    {
        zai        = zarg*i;
        sinzai     = sin(zai);
        infl       = do_p3m_influence(sinzai, order);
        bsp_mod[i] = infl*infl*pow(sinzai/zai, -2.0*order);
    }
}

/* Calculate the P3M B-spline moduli */
static void make_p3m_bspline_moduli(splinevec bsp_mod,
                                    int nx, int ny, int nz, int order)
{
    make_p3m_bspline_moduli_dim(bsp_mod[XX], nx, order);
    make_p3m_bspline_moduli_dim(bsp_mod[YY], ny, order);
    make_p3m_bspline_moduli_dim(bsp_mod[ZZ], nz, order);
}


static void setup_coordinate_communication(pme_atomcomm_t *atc)
{
    int nslab, n, i;
    int fw, bw;

    nslab = atc->nslab;

    n = 0;
    for (i = 1; i <= nslab/2; i++)
    {
        fw = (atc->nodeid + i) % nslab;
        bw = (atc->nodeid - i + nslab) % nslab;
        if (n < nslab - 1)
        {
            atc->node_dest[n] = fw;
            atc->node_src[n]  = bw;
            n++;
        }
        if (n < nslab - 1)
        {
            atc->node_dest[n] = bw;
            atc->node_src[n]  = fw;
            n++;
        }
    }
}

int gmx_pme_destroy(FILE *log, gmx_pme_t *pmedata)
{
    int thread, i;

    if (NULL != log)
    {
        fprintf(log, "Destroying PME data structures.\n");
    }

    sfree((*pmedata)->nnx);
    sfree((*pmedata)->nny);
    sfree((*pmedata)->nnz);

    for (i = 0; i < (*pmedata)->ngrids; ++i)
    {
        pmegrids_destroy(&(*pmedata)->pmegrid[i]);
        sfree((*pmedata)->fftgrid[i]);
        sfree((*pmedata)->cfftgrid[i]);
        gmx_parallel_3dfft_destroy((*pmedata)->pfft_setup[i]);
    }

    sfree((*pmedata)->lb_buf1);
    sfree((*pmedata)->lb_buf2);

    for (thread = 0; thread < (*pmedata)->nthread; thread++)
    {
        free_work(&(*pmedata)->work[thread]);
    }
    sfree((*pmedata)->work);

    sfree(*pmedata);
    *pmedata = NULL;

    return 0;
}

static int mult_up(int n, int f)
{
    return ((n + f - 1)/f)*f;
}


static double pme_load_imbalance(gmx_pme_t pme)
{
    int    nma, nmi;
    double n1, n2, n3;

    nma = pme->nnodes_major;
    nmi = pme->nnodes_minor;

    n1 = mult_up(pme->nkx, nma)*mult_up(pme->nky, nmi)*pme->nkz;
    n2 = mult_up(pme->nkx, nma)*mult_up(pme->nkz, nmi)*pme->nky;
    n3 = mult_up(pme->nky, nma)*mult_up(pme->nkz, nmi)*pme->nkx;

    /* pme_solve is roughly double the cost of an fft */

    return (n1 + n2 + 3*n3)/(double)(6*pme->nkx*pme->nky*pme->nkz);
}

static void init_atomcomm(gmx_pme_t pme, pme_atomcomm_t *atc,
                          int dimind, gmx_bool bSpread)
{
    int nk, k, s, thread;

    atc->dimind    = dimind;
    atc->nslab     = 1;
    atc->nodeid    = 0;
    atc->pd_nalloc = 0;
#ifdef GMX_MPI
    if (pme->nnodes > 1)
    {
        atc->mpi_comm = pme->mpi_comm_d[dimind];
        MPI_Comm_size(atc->mpi_comm, &atc->nslab);
        MPI_Comm_rank(atc->mpi_comm, &atc->nodeid);
    }
    if (debug)
    {
        fprintf(debug, "For PME atom communication in dimind %d: nslab %d rank %d\n", atc->dimind, atc->nslab, atc->nodeid);
    }
#endif

    atc->bSpread   = bSpread;
    atc->pme_order = pme->pme_order;

    if (atc->nslab > 1)
    {
        snew(atc->node_dest, atc->nslab);
        snew(atc->node_src, atc->nslab);
        setup_coordinate_communication(atc);

        snew(atc->count_thread, pme->nthread);
        for (thread = 0; thread < pme->nthread; thread++)
        {
            snew(atc->count_thread[thread], atc->nslab);
        }
        atc->count = atc->count_thread[0];
        snew(atc->rcount, atc->nslab);
        snew(atc->buf_index, atc->nslab);
    }

    atc->nthread = pme->nthread;
    if (atc->nthread > 1)
    {
        snew(atc->thread_plist, atc->nthread);
    }
    snew(atc->spline, atc->nthread);
    for (thread = 0; thread < atc->nthread; thread++)
    {
        if (atc->nthread > 1)
        {
            snew(atc->thread_plist[thread].n, atc->nthread+2*GMX_CACHE_SEP);
            atc->thread_plist[thread].n += GMX_CACHE_SEP;
        }
        snew(atc->spline[thread].thread_one, pme->nthread);
        atc->spline[thread].thread_one[thread] = 1;
    }
}

static void
init_overlap_comm(pme_overlap_t *  ol,
                  int              norder,
#ifdef GMX_MPI
                  MPI_Comm         comm,
#endif
                  int              nnodes,
                  int              nodeid,
                  int              ndata,
                  int              commplainsize)
{
    int lbnd, rbnd, maxlr, b, i;
    int exten;
    int nn, nk;
    pme_grid_comm_t *pgc;
    gmx_bool bCont;
    int fft_start, fft_end, send_index1, recv_index1;
#ifdef GMX_MPI
    MPI_Status stat;

    ol->mpi_comm = comm;
#endif

    ol->nnodes = nnodes;
    ol->nodeid = nodeid;

    /* Linear translation of the PME grid won't affect reciprocal space
     * calculations, so to optimize we only interpolate "upwards",
     * which also means we only have to consider overlap in one direction.
     * I.e., particles on this node might also be spread to grid indices
     * that belong to higher nodes (modulo nnodes)
     */

    snew(ol->s2g0, ol->nnodes+1);
    snew(ol->s2g1, ol->nnodes);
    if (debug)
    {
        fprintf(debug, "PME slab boundaries:");
    }
    for (i = 0; i < nnodes; i++)
    {
        /* s2g0 the local interpolation grid start.
         * s2g1 the local interpolation grid end.
         * Since in calc_pidx we divide particles, and not grid lines,
         * spatially uniform along dimension x or y, we need to round
         * s2g0 down and s2g1 up.
         */
        ol->s2g0[i] = ( i   *ndata + 0       )/nnodes;
        ol->s2g1[i] = ((i+1)*ndata + nnodes-1)/nnodes + norder - 1;

        if (debug)
        {
            fprintf(debug, "  %3d %3d", ol->s2g0[i], ol->s2g1[i]);
        }
    }
    ol->s2g0[nnodes] = ndata;
    if (debug)
    {
        fprintf(debug, "\n");
    }

    /* Determine with how many nodes we need to communicate the grid overlap */
    b = 0;
    do
    {
        b++;
        bCont = FALSE;
        for (i = 0; i < nnodes; i++)
        {
            if ((i+b <  nnodes && ol->s2g1[i] > ol->s2g0[i+b]) ||
                (i+b >= nnodes && ol->s2g1[i] > ol->s2g0[i+b-nnodes] + ndata))
            {
                bCont = TRUE;
            }
        }
    }
    while (bCont && b < nnodes);
    ol->noverlap_nodes = b - 1;

    snew(ol->send_id, ol->noverlap_nodes);
    snew(ol->recv_id, ol->noverlap_nodes);
    for (b = 0; b < ol->noverlap_nodes; b++)
    {
        ol->send_id[b] = (ol->nodeid + (b + 1)) % ol->nnodes;
        ol->recv_id[b] = (ol->nodeid - (b + 1) + ol->nnodes) % ol->nnodes;
    }
    snew(ol->comm_data, ol->noverlap_nodes);

    ol->send_size = 0;
    for (b = 0; b < ol->noverlap_nodes; b++)
    {
        pgc = &ol->comm_data[b];
        /* Send */
        fft_start        = ol->s2g0[ol->send_id[b]];
        fft_end          = ol->s2g0[ol->send_id[b]+1];
        if (ol->send_id[b] < nodeid)
        {
            fft_start += ndata;
            fft_end   += ndata;
        }
        send_index1       = ol->s2g1[nodeid];
        send_index1       = min(send_index1, fft_end);
        pgc->send_index0  = fft_start;
        pgc->send_nindex  = max(0, send_index1 - pgc->send_index0);
        ol->send_size    += pgc->send_nindex;

        /* We always start receiving to the first index of our slab */
        fft_start        = ol->s2g0[ol->nodeid];
        fft_end          = ol->s2g0[ol->nodeid+1];
        recv_index1      = ol->s2g1[ol->recv_id[b]];
        if (ol->recv_id[b] > nodeid)
        {
            recv_index1 -= ndata;
        }
        recv_index1      = min(recv_index1, fft_end);
        pgc->recv_index0 = fft_start;
        pgc->recv_nindex = max(0, recv_index1 - pgc->recv_index0);
    }

#ifdef GMX_MPI
    /* Communicate the buffer sizes to receive */
    for (b = 0; b < ol->noverlap_nodes; b++)
    {
        MPI_Sendrecv(&ol->send_size, 1, MPI_INT, ol->send_id[b], b,
                     &ol->comm_data[b].recv_size, 1, MPI_INT, ol->recv_id[b], b,
                     ol->mpi_comm, &stat);
    }
#endif

    /* For non-divisible grid we need pme_order iso pme_order-1 */
    snew(ol->sendbuf, norder*commplainsize);
    snew(ol->recvbuf, norder*commplainsize);
}

static void
make_gridindex5_to_localindex(int n, int local_start, int local_range,
                              int **global_to_local,
                              real **fraction_shift)
{
    int i;
    int * gtl;
    real * fsh;

    snew(gtl, 5*n);
    snew(fsh, 5*n);
    for (i = 0; (i < 5*n); i++)
    {
        /* Determine the global to local grid index */
        gtl[i] = (i - local_start + n) % n;
        /* For coordinates that fall within the local grid the fraction
         * is correct, we don't need to shift it.
         */
        fsh[i] = 0;
        if (local_range < n)
        {
            /* Due to rounding issues i could be 1 beyond the lower or
             * upper boundary of the local grid. Correct the index for this.
             * If we shift the index, we need to shift the fraction by
             * the same amount in the other direction to not affect
             * the weights.
             * Note that due to this shifting the weights at the end of
             * the spline might change, but that will only involve values
             * between zero and values close to the precision of a real,
             * which is anyhow the accuracy of the whole mesh calculation.
             */
            /* With local_range=0 we should not change i=local_start */
            if (i % n != local_start)
            {
                if (gtl[i] == n-1)
                {
                    gtl[i] = 0;
                    fsh[i] = -1;
                }
                else if (gtl[i] == local_range)
                {
                    gtl[i] = local_range - 1;
                    fsh[i] = 1;
                }
            }
        }
    }

    *global_to_local = gtl;
    *fraction_shift  = fsh;
}

static pme_spline_work_t *make_pme_spline_work(int gmx_unused order)
{
    pme_spline_work_t *work;

#ifdef PME_SIMD4_SPREAD_GATHER
    real             tmp[GMX_SIMD4_WIDTH*3], *tmp_aligned;
    gmx_simd4_real_t zero_S;
    gmx_simd4_real_t real_mask_S0, real_mask_S1;
    int              of, i;

    snew_aligned(work, 1, SIMD4_ALIGNMENT);

    tmp_aligned = gmx_simd4_align_r(tmp);

    zero_S = gmx_simd4_setzero_r();

    /* Generate bit masks to mask out the unused grid entries,
     * as we only operate on order of the 8 grid entries that are
     * load into 2 SIMD registers.
     */
    for (of = 0; of < 2*GMX_SIMD4_WIDTH-(order-1); of++)
    {
        for (i = 0; i < 2*GMX_SIMD4_WIDTH; i++)
        {
            tmp_aligned[i] = (i >= of && i < of+order ? -1.0 : 1.0);
        }
        real_mask_S0      = gmx_simd4_load_r(tmp_aligned);
        real_mask_S1      = gmx_simd4_load_r(tmp_aligned+GMX_SIMD4_WIDTH);
        work->mask_S0[of] = gmx_simd4_cmplt_r(real_mask_S0, zero_S);
        work->mask_S1[of] = gmx_simd4_cmplt_r(real_mask_S1, zero_S);
    }
#else
    work = NULL;
#endif

    return work;
}

void gmx_pme_check_restrictions(int pme_order,
                                int nkx, int nky, int nkz,
                                int nnodes_major,
                                int nnodes_minor,
                                gmx_bool bUseThreads,
                                gmx_bool bFatal,
                                gmx_bool *bValidSettings)
{
    if (pme_order > PME_ORDER_MAX)
    {
        if (!bFatal)
        {
            *bValidSettings = FALSE;
            return;
        }
        gmx_fatal(FARGS, "pme_order (%d) is larger than the maximum allowed value (%d). Modify and recompile the code if you really need such a high order.",
                  pme_order, PME_ORDER_MAX);
    }

    if (nkx <= pme_order*(nnodes_major > 1 ? 2 : 1) ||
        nky <= pme_order*(nnodes_minor > 1 ? 2 : 1) ||
        nkz <= pme_order)
    {
        if (!bFatal)
        {
            *bValidSettings = FALSE;
            return;
        }
        gmx_fatal(FARGS, "The PME grid sizes need to be larger than pme_order (%d) and for dimensions with domain decomposition larger than 2*pme_order",
                  pme_order);
    }

    /* Check for a limitation of the (current) sum_fftgrid_dd code.
     * We only allow multiple communication pulses in dim 1, not in dim 0.
     */
    if (bUseThreads && (nkx < nnodes_major*pme_order &&
                        nkx != nnodes_major*(pme_order - 1)))
    {
        if (!bFatal)
        {
            *bValidSettings = FALSE;
            return;
        }
        gmx_fatal(FARGS, "The number of PME grid lines per rank along x is %g. But when using OpenMP threads, the number of grid lines per rank along x should be >= pme_order (%d) or = pmeorder-1. To resolve this issue, use fewer ranks along x (and possibly more along y and/or z) by specifying -dd manually.",
                  nkx/(double)nnodes_major, pme_order);
    }

    if (bValidSettings != NULL)
    {
        *bValidSettings = TRUE;
    }

    return;
}

int gmx_pme_init(gmx_pme_t *         pmedata,
                 t_commrec *         cr,
                 int                 nnodes_major,
                 int                 nnodes_minor,
                 t_inputrec *        ir,
                 int                 homenr,
                 gmx_bool            bFreeEnergy_q,
                 gmx_bool            bFreeEnergy_lj,
                 gmx_bool            bReproducible,
                 int                 nthread)
{
    gmx_pme_t pme = NULL;

    int  use_threads, sum_use_threads, i;
    ivec ndata;

    if (debug)
    {
        fprintf(debug, "Creating PME data structures.\n");
    }
    snew(pme, 1);

    pme->sum_qgrid_tmp       = NULL;
    pme->sum_qgrid_dd_tmp    = NULL;
    pme->buf_nalloc          = 0;

    pme->nnodes              = 1;
    pme->bPPnode             = TRUE;

    pme->nnodes_major        = nnodes_major;
    pme->nnodes_minor        = nnodes_minor;

#ifdef GMX_MPI
    if (nnodes_major*nnodes_minor > 1)
    {
        pme->mpi_comm = cr->mpi_comm_mygroup;

        MPI_Comm_rank(pme->mpi_comm, &pme->nodeid);
        MPI_Comm_size(pme->mpi_comm, &pme->nnodes);
        if (pme->nnodes != nnodes_major*nnodes_minor)
        {
            gmx_incons("PME rank count mismatch");
        }
    }
    else
    {
        pme->mpi_comm = MPI_COMM_NULL;
    }
#endif

    if (pme->nnodes == 1)
    {
#ifdef GMX_MPI
        pme->mpi_comm_d[0] = MPI_COMM_NULL;
        pme->mpi_comm_d[1] = MPI_COMM_NULL;
#endif
        pme->ndecompdim   = 0;
        pme->nodeid_major = 0;
        pme->nodeid_minor = 0;
#ifdef GMX_MPI
        pme->mpi_comm_d[0] = pme->mpi_comm_d[1] = MPI_COMM_NULL;
#endif
    }
    else
    {
        if (nnodes_minor == 1)
        {
#ifdef GMX_MPI
            pme->mpi_comm_d[0] = pme->mpi_comm;
            pme->mpi_comm_d[1] = MPI_COMM_NULL;
#endif
            pme->ndecompdim   = 1;
            pme->nodeid_major = pme->nodeid;
            pme->nodeid_minor = 0;

        }
        else if (nnodes_major == 1)
        {
#ifdef GMX_MPI
            pme->mpi_comm_d[0] = MPI_COMM_NULL;
            pme->mpi_comm_d[1] = pme->mpi_comm;
#endif
            pme->ndecompdim   = 1;
            pme->nodeid_major = 0;
            pme->nodeid_minor = pme->nodeid;
        }
        else
        {
            if (pme->nnodes % nnodes_major != 0)
            {
                gmx_incons("For 2D PME decomposition, #PME ranks must be divisible by the number of ranks in the major dimension");
            }
            pme->ndecompdim = 2;

#ifdef GMX_MPI
            MPI_Comm_split(pme->mpi_comm, pme->nodeid % nnodes_minor,
                           pme->nodeid, &pme->mpi_comm_d[0]);  /* My communicator along major dimension */
            MPI_Comm_split(pme->mpi_comm, pme->nodeid/nnodes_minor,
                           pme->nodeid, &pme->mpi_comm_d[1]);  /* My communicator along minor dimension */

            MPI_Comm_rank(pme->mpi_comm_d[0], &pme->nodeid_major);
            MPI_Comm_size(pme->mpi_comm_d[0], &pme->nnodes_major);
            MPI_Comm_rank(pme->mpi_comm_d[1], &pme->nodeid_minor);
            MPI_Comm_size(pme->mpi_comm_d[1], &pme->nnodes_minor);
#endif
        }
        pme->bPPnode = (cr->duty & DUTY_PP);
    }

    pme->nthread = nthread;

    /* Check if any of the PME MPI ranks uses threads */
    use_threads = (pme->nthread > 1 ? 1 : 0);
#ifdef GMX_MPI
    if (pme->nnodes > 1)
    {
        MPI_Allreduce(&use_threads, &sum_use_threads, 1, MPI_INT,
                      MPI_SUM, pme->mpi_comm);
    }
    else
#endif
    {
        sum_use_threads = use_threads;
    }
    pme->bUseThreads = (sum_use_threads > 0);

    if (ir->ePBC == epbcSCREW)
    {
        gmx_fatal(FARGS, "pme does not (yet) work with pbc = screw");
    }

    pme->bFEP_q      = ((ir->efep != efepNO) && bFreeEnergy_q);
    pme->bFEP_lj     = ((ir->efep != efepNO) && bFreeEnergy_lj);
    pme->bFEP        = (pme->bFEP_q || pme->bFEP_lj);
    pme->nkx         = ir->nkx;
    pme->nky         = ir->nky;
    pme->nkz         = ir->nkz;
    pme->bP3M        = (ir->coulombtype == eelP3M_AD || getenv("GMX_PME_P3M") != NULL);
    pme->pme_order   = ir->pme_order;

    /* Always constant electrostatics coefficients */
    pme->epsilon_r   = ir->epsilon_r;

    /* Always constant LJ coefficients */
    pme->ljpme_combination_rule = ir->ljpme_combination_rule;

    /* If we violate restrictions, generate a fatal error here */
    gmx_pme_check_restrictions(pme->pme_order,
                               pme->nkx, pme->nky, pme->nkz,
                               pme->nnodes_major,
                               pme->nnodes_minor,
                               pme->bUseThreads,
                               TRUE,
                               NULL);

    if (pme->nnodes > 1)
    {
        double imbal;

#ifdef GMX_MPI
        MPI_Type_contiguous(DIM, mpi_type, &(pme->rvec_mpi));
        MPI_Type_commit(&(pme->rvec_mpi));
#endif

        /* Note that the coefficient spreading and force gathering, which usually
         * takes about the same amount of time as FFT+solve_pme,
         * is always fully load balanced
         * (unless the coefficient distribution is inhomogeneous).
         */

        imbal = pme_load_imbalance(pme);
        if (imbal >= 1.2 && pme->nodeid_major == 0 && pme->nodeid_minor == 0)
        {
            fprintf(stderr,
                    "\n"
                    "NOTE: The load imbalance in PME FFT and solve is %d%%.\n"
                    "      For optimal PME load balancing\n"
                    "      PME grid_x (%d) and grid_y (%d) should be divisible by #PME_ranks_x (%d)\n"
                    "      and PME grid_y (%d) and grid_z (%d) should be divisible by #PME_ranks_y (%d)\n"
                    "\n",
                    (int)((imbal-1)*100 + 0.5),
                    pme->nkx, pme->nky, pme->nnodes_major,
                    pme->nky, pme->nkz, pme->nnodes_minor);
        }
    }

    /* For non-divisible grid we need pme_order iso pme_order-1 */
    /* In sum_qgrid_dd x overlap is copied in place: take padding into account.
     * y is always copied through a buffer: we don't need padding in z,
     * but we do need the overlap in x because of the communication order.
     */
    init_overlap_comm(&pme->overlap[0], pme->pme_order,
#ifdef GMX_MPI
                      pme->mpi_comm_d[0],
#endif
                      pme->nnodes_major, pme->nodeid_major,
                      pme->nkx,
                      (div_round_up(pme->nky, pme->nnodes_minor)+pme->pme_order)*(pme->nkz+pme->pme_order-1));

    /* Along overlap dim 1 we can send in multiple pulses in sum_fftgrid_dd.
     * We do this with an offset buffer of equal size, so we need to allocate
     * extra for the offset. That's what the (+1)*pme->nkz is for.
     */
    init_overlap_comm(&pme->overlap[1], pme->pme_order,
#ifdef GMX_MPI
                      pme->mpi_comm_d[1],
#endif
                      pme->nnodes_minor, pme->nodeid_minor,
                      pme->nky,
                      (div_round_up(pme->nkx, pme->nnodes_major)+pme->pme_order+1)*pme->nkz);

    /* Double-check for a limitation of the (current) sum_fftgrid_dd code.
     * Note that gmx_pme_check_restrictions checked for this already.
     */
    if (pme->bUseThreads && pme->overlap[0].noverlap_nodes > 1)
    {
        gmx_incons("More than one communication pulse required for grid overlap communication along the major dimension while using threads");
    }

    snew(pme->bsp_mod[XX], pme->nkx);
    snew(pme->bsp_mod[YY], pme->nky);
    snew(pme->bsp_mod[ZZ], pme->nkz);

    /* The required size of the interpolation grid, including overlap.
     * The allocated size (pmegrid_n?) might be slightly larger.
     */
    pme->pmegrid_nx = pme->overlap[0].s2g1[pme->nodeid_major] -
        pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_ny = pme->overlap[1].s2g1[pme->nodeid_minor] -
        pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_nz_base = pme->nkz;
    pme->pmegrid_nz      = pme->pmegrid_nz_base + pme->pme_order - 1;
    set_grid_alignment(&pme->pmegrid_nz, pme->pme_order);

    pme->pmegrid_start_ix = pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_start_iy = pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_start_iz = 0;

    make_gridindex5_to_localindex(pme->nkx,
                                  pme->pmegrid_start_ix,
                                  pme->pmegrid_nx - (pme->pme_order-1),
                                  &pme->nnx, &pme->fshx);
    make_gridindex5_to_localindex(pme->nky,
                                  pme->pmegrid_start_iy,
                                  pme->pmegrid_ny - (pme->pme_order-1),
                                  &pme->nny, &pme->fshy);
    make_gridindex5_to_localindex(pme->nkz,
                                  pme->pmegrid_start_iz,
                                  pme->pmegrid_nz_base,
                                  &pme->nnz, &pme->fshz);

    pme->spline_work = make_pme_spline_work(pme->pme_order);

    ndata[0]    = pme->nkx;
    ndata[1]    = pme->nky;
    ndata[2]    = pme->nkz;
    /* It doesn't matter if we allocate too many grids here,
     * we only allocate and use the ones we need.
     */
    if (EVDW_PME(ir->vdwtype))
    {
        pme->ngrids = ((ir->ljpme_combination_rule == eljpmeLB) ? DO_Q_AND_LJ_LB : DO_Q_AND_LJ);
    }
    else
    {
        pme->ngrids = DO_Q;
    }
    snew(pme->fftgrid, pme->ngrids);
    snew(pme->cfftgrid, pme->ngrids);
    snew(pme->pfft_setup, pme->ngrids);

    for (i = 0; i < pme->ngrids; ++i)
    {
        if ((i <  DO_Q && EEL_PME(ir->coulombtype) && (i == 0 ||
                                                       bFreeEnergy_q)) ||
            (i >= DO_Q && EVDW_PME(ir->vdwtype) && (i == 2 ||
                                                    bFreeEnergy_lj ||
                                                    ir->ljpme_combination_rule == eljpmeLB)))
        {
            pmegrids_init(&pme->pmegrid[i],
                          pme->pmegrid_nx, pme->pmegrid_ny, pme->pmegrid_nz,
                          pme->pmegrid_nz_base,
                          pme->pme_order,
                          pme->bUseThreads,
                          pme->nthread,
                          pme->overlap[0].s2g1[pme->nodeid_major]-pme->overlap[0].s2g0[pme->nodeid_major+1],
                          pme->overlap[1].s2g1[pme->nodeid_minor]-pme->overlap[1].s2g0[pme->nodeid_minor+1]);
            /* This routine will allocate the grid data to fit the FFTs */
            gmx_parallel_3dfft_init(&pme->pfft_setup[i], ndata,
                                    &pme->fftgrid[i], &pme->cfftgrid[i],
                                    pme->mpi_comm_d,
                                    bReproducible, pme->nthread);

        }
    }

    if (!pme->bP3M)
    {
        /* Use plain SPME B-spline interpolation */
        make_bspline_moduli(pme->bsp_mod, pme->nkx, pme->nky, pme->nkz, pme->pme_order);
    }
    else
    {
        /* Use the P3M grid-optimized influence function */
        make_p3m_bspline_moduli(pme->bsp_mod, pme->nkx, pme->nky, pme->nkz, pme->pme_order);
    }

    /* Use atc[0] for spreading */
    init_atomcomm(pme, &pme->atc[0], nnodes_major > 1 ? 0 : 1, TRUE);
    if (pme->ndecompdim >= 2)
    {
        init_atomcomm(pme, &pme->atc[1], 1, FALSE);
    }

    if (pme->nnodes == 1)
    {
        pme->atc[0].n = homenr;
        pme_realloc_atomcomm_things(&pme->atc[0]);
    }

    pme->lb_buf1       = NULL;
    pme->lb_buf2       = NULL;
    pme->lb_buf_nalloc = 0;

    {
        int thread;

        /* Use fft5d, order after FFT is y major, z, x minor */

        snew(pme->work, pme->nthread);
        for (thread = 0; thread < pme->nthread; thread++)
        {
            realloc_work(&pme->work[thread], pme->nkx);
        }
    }

    *pmedata = pme;

    return 0;
}

static void reuse_pmegrids(const pmegrids_t *old, pmegrids_t *new)
{
    int d, t;

    for (d = 0; d < DIM; d++)
    {
        if (new->grid.n[d] > old->grid.n[d])
        {
            return;
        }
    }

    sfree_aligned(new->grid.grid);
    new->grid.grid = old->grid.grid;

    if (new->grid_th != NULL && new->nthread == old->nthread)
    {
        sfree_aligned(new->grid_all);
        for (t = 0; t < new->nthread; t++)
        {
            new->grid_th[t].grid = old->grid_th[t].grid;
        }
    }
}

int gmx_pme_reinit(gmx_pme_t *         pmedata,
                   t_commrec *         cr,
                   gmx_pme_t           pme_src,
                   const t_inputrec *  ir,
                   ivec                grid_size)
{
    t_inputrec irc;
    int homenr;
    int ret;

    irc     = *ir;
    irc.nkx = grid_size[XX];
    irc.nky = grid_size[YY];
    irc.nkz = grid_size[ZZ];

    if (pme_src->nnodes == 1)
    {
        homenr = pme_src->atc[0].n;
    }
    else
    {
        homenr = -1;
    }

    ret = gmx_pme_init(pmedata, cr, pme_src->nnodes_major, pme_src->nnodes_minor,
                       &irc, homenr, pme_src->bFEP_q, pme_src->bFEP_lj, FALSE, pme_src->nthread);

    if (ret == 0)
    {
        /* We can easily reuse the allocated pme grids in pme_src */
        reuse_pmegrids(&pme_src->pmegrid[PME_GRID_QA], &(*pmedata)->pmegrid[PME_GRID_QA]);
        /* We would like to reuse the fft grids, but that's harder */
    }

    return ret;
}


static void copy_local_grid(gmx_pme_t pme, pmegrids_t *pmegrids,
                            int grid_index, int thread, real *fftgrid)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  fft_my, fft_mz;
    int  nsx, nsy, nsz;
    ivec nf;
    int  offx, offy, offz, x, y, z, i0, i0t;
    int  d;
    pmegrid_t *pmegrid;
    real *grid_th;

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);
    fft_my = local_fft_size[YY];
    fft_mz = local_fft_size[ZZ];

    pmegrid = &pmegrids->grid_th[thread];

    nsx = pmegrid->s[XX];
    nsy = pmegrid->s[YY];
    nsz = pmegrid->s[ZZ];

    for (d = 0; d < DIM; d++)
    {
        nf[d] = min(pmegrid->n[d] - (pmegrid->order - 1),
                    local_fft_ndata[d] - pmegrid->offset[d]);
    }

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    /* Directly copy the non-overlapping parts of the local grids.
     * This also initializes the full grid.
     */
    grid_th = pmegrid->grid;
    for (x = 0; x < nf[XX]; x++)
    {
        for (y = 0; y < nf[YY]; y++)
        {
            i0  = ((offx + x)*fft_my + (offy + y))*fft_mz + offz;
            i0t = (x*nsy + y)*nsz;
            for (z = 0; z < nf[ZZ]; z++)
            {
                fftgrid[i0+z] = grid_th[i0t+z];
            }
        }
    }
}

static void
reduce_threadgrid_overlap(gmx_pme_t pme,
                          const pmegrids_t *pmegrids, int thread,
                          real *fftgrid, real *commbuf_x, real *commbuf_y,
                          int grid_index)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    int  fft_nx, fft_ny, fft_nz;
    int  fft_my, fft_mz;
    int  buf_my = -1;
    int  nsx, nsy, nsz;
    ivec localcopy_end;
    int  offx, offy, offz, x, y, z, i0, i0t;
    int  sx, sy, sz, fx, fy, fz, tx1, ty1, tz1, ox, oy, oz;
    gmx_bool bClearBufX, bClearBufY, bClearBufXY, bClearBuf;
    gmx_bool bCommX, bCommY;
    int  d;
    int  thread_f;
    const pmegrid_t *pmegrid, *pmegrid_g, *pmegrid_f;
    const real *grid_th;
    real *commbuf = NULL;

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);
    fft_nx = local_fft_ndata[XX];
    fft_ny = local_fft_ndata[YY];
    fft_nz = local_fft_ndata[ZZ];

    fft_my = local_fft_size[YY];
    fft_mz = local_fft_size[ZZ];

    /* This routine is called when all thread have finished spreading.
     * Here each thread sums grid contributions calculated by other threads
     * to the thread local grid volume.
     * To minimize the number of grid copying operations,
     * this routines sums immediately from the pmegrid to the fftgrid.
     */

    /* Determine which part of the full node grid we should operate on,
     * this is our thread local part of the full grid.
     */
    pmegrid = &pmegrids->grid_th[thread];

    for (d = 0; d < DIM; d++)
    {
        /* Determine up to where our thread needs to copy from the
         * thread-local charge spreading grid to the rank-local FFT grid.
         * This is up to our spreading grid end minus order-1 and
         * not beyond the local FFT grid.
         */
        localcopy_end[d] =
            min(pmegrid->offset[d]+pmegrid->n[d]-(pmegrid->order-1),
                local_fft_ndata[d]);
    }

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];


    bClearBufX  = TRUE;
    bClearBufY  = TRUE;
    bClearBufXY = TRUE;

    /* Now loop over all the thread data blocks that contribute
     * to the grid region we (our thread) are operating on.
     */
    /* Note that fft_nx/y is equal to the number of grid points
     * between the first point of our node grid and the one of the next node.
     */
    for (sx = 0; sx >= -pmegrids->nthread_comm[XX]; sx--)
    {
        fx     = pmegrid->ci[XX] + sx;
        ox     = 0;
        bCommX = FALSE;
        if (fx < 0)
        {
            fx    += pmegrids->nc[XX];
            ox    -= fft_nx;
            bCommX = (pme->nnodes_major > 1);
        }
        pmegrid_g = &pmegrids->grid_th[fx*pmegrids->nc[YY]*pmegrids->nc[ZZ]];
        ox       += pmegrid_g->offset[XX];
        /* Determine the end of our part of the source grid */
        if (!bCommX)
        {
            /* Use our thread local source grid and target grid part */
            tx1 = min(ox + pmegrid_g->n[XX], localcopy_end[XX]);
        }
        else
        {
            /* Use our thread local source grid and the spreading range */
            tx1 = min(ox + pmegrid_g->n[XX], pme->pme_order);
        }

        for (sy = 0; sy >= -pmegrids->nthread_comm[YY]; sy--)
        {
            fy     = pmegrid->ci[YY] + sy;
            oy     = 0;
            bCommY = FALSE;
            if (fy < 0)
            {
                fy    += pmegrids->nc[YY];
                oy    -= fft_ny;
                bCommY = (pme->nnodes_minor > 1);
            }
            pmegrid_g = &pmegrids->grid_th[fy*pmegrids->nc[ZZ]];
            oy       += pmegrid_g->offset[YY];
            /* Determine the end of our part of the source grid */
            if (!bCommY)
            {
                /* Use our thread local source grid and target grid part */
                ty1 = min(oy + pmegrid_g->n[YY], localcopy_end[YY]);
            }
            else
            {
                /* Use our thread local source grid and the spreading range */
                ty1 = min(oy + pmegrid_g->n[YY], pme->pme_order);
            }

            for (sz = 0; sz >= -pmegrids->nthread_comm[ZZ]; sz--)
            {
                fz = pmegrid->ci[ZZ] + sz;
                oz = 0;
                if (fz < 0)
                {
                    fz += pmegrids->nc[ZZ];
                    oz -= fft_nz;
                }
                pmegrid_g = &pmegrids->grid_th[fz];
                oz       += pmegrid_g->offset[ZZ];
                tz1       = min(oz + pmegrid_g->n[ZZ], localcopy_end[ZZ]);

                if (sx == 0 && sy == 0 && sz == 0)
                {
                    /* We have already added our local contribution
                     * before calling this routine, so skip it here.
                     */
                    continue;
                }

                thread_f = (fx*pmegrids->nc[YY] + fy)*pmegrids->nc[ZZ] + fz;

                pmegrid_f = &pmegrids->grid_th[thread_f];

                grid_th = pmegrid_f->grid;

                nsx = pmegrid_f->s[XX];
                nsy = pmegrid_f->s[YY];
                nsz = pmegrid_f->s[ZZ];

#ifdef DEBUG_PME_REDUCE
                printf("n%d t%d add %d  %2d %2d %2d  %2d %2d %2d  %2d-%2d %2d-%2d, %2d-%2d %2d-%2d, %2d-%2d %2d-%2d\n",
                       pme->nodeid, thread, thread_f,
                       pme->pmegrid_start_ix,
                       pme->pmegrid_start_iy,
                       pme->pmegrid_start_iz,
                       sx, sy, sz,
                       offx-ox, tx1-ox, offx, tx1,
                       offy-oy, ty1-oy, offy, ty1,
                       offz-oz, tz1-oz, offz, tz1);
#endif

                if (!(bCommX || bCommY))
                {
                    /* Copy from the thread local grid to the node grid */
                    for (x = offx; x < tx1; x++)
                    {
                        for (y = offy; y < ty1; y++)
                        {
                            i0  = (x*fft_my + y)*fft_mz;
                            i0t = ((x - ox)*nsy + (y - oy))*nsz - oz;
                            for (z = offz; z < tz1; z++)
                            {
                                fftgrid[i0+z] += grid_th[i0t+z];
                            }
                        }
                    }
                }
                else
                {
                    /* The order of this conditional decides
                     * where the corner volume gets stored with x+y decomp.
                     */
                    if (bCommY)
                    {
                        commbuf = commbuf_y;
                        /* The y-size of the communication buffer is set by
                         * the overlap of the grid part of our local slab
                         * with the part starting at the next slab.
                         */
                        buf_my  =
                            pme->overlap[1].s2g1[pme->nodeid_minor] -
                            pme->overlap[1].s2g0[pme->nodeid_minor+1];
                        if (bCommX)
                        {
                            /* We index commbuf modulo the local grid size */
                            commbuf += buf_my*fft_nx*fft_nz;

                            bClearBuf   = bClearBufXY;
                            bClearBufXY = FALSE;
                        }
                        else
                        {
                            bClearBuf  = bClearBufY;
                            bClearBufY = FALSE;
                        }
                    }
                    else
                    {
                        commbuf    = commbuf_x;
                        buf_my     = fft_ny;
                        bClearBuf  = bClearBufX;
                        bClearBufX = FALSE;
                    }

                    /* Copy to the communication buffer */
                    for (x = offx; x < tx1; x++)
                    {
                        for (y = offy; y < ty1; y++)
                        {
                            i0  = (x*buf_my + y)*fft_nz;
                            i0t = ((x - ox)*nsy + (y - oy))*nsz - oz;

                            if (bClearBuf)
                            {
                                /* First access of commbuf, initialize it */
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0+z]  = grid_th[i0t+z];
                                }
                            }
                            else
                            {
                                for (z = offz; z < tz1; z++)
                                {
                                    commbuf[i0+z] += grid_th[i0t+z];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


static void sum_fftgrid_dd(gmx_pme_t pme, real *fftgrid, int grid_index)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;
    pme_overlap_t *overlap;
    int  send_index0, send_nindex;
    int  recv_nindex;
#ifdef GMX_MPI
    MPI_Status stat;
#endif
    int  send_size_y, recv_size_y;
    int  ipulse, send_id, recv_id, datasize, gridsize, size_yx;
    real *sendptr, *recvptr;
    int  x, y, z, indg, indb;

    /* Note that this routine is only used for forward communication.
     * Since the force gathering, unlike the coefficient spreading,
     * can be trivially parallelized over the particles,
     * the backwards process is much simpler and can use the "old"
     * communication setup.
     */

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[grid_index],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    if (pme->nnodes_minor > 1)
    {
        /* Major dimension */
        overlap = &pme->overlap[1];

        if (pme->nnodes_major > 1)
        {
            size_yx = pme->overlap[0].comm_data[0].send_nindex;
        }
        else
        {
            size_yx = 0;
        }
        datasize = (local_fft_ndata[XX] + size_yx)*local_fft_ndata[ZZ];

        send_size_y = overlap->send_size;

        for (ipulse = 0; ipulse < overlap->noverlap_nodes; ipulse++)
        {
            send_id       = overlap->send_id[ipulse];
            recv_id       = overlap->recv_id[ipulse];
            send_index0   =
                overlap->comm_data[ipulse].send_index0 -
                overlap->comm_data[0].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            /* We don't use recv_index0, as we always receive starting at 0 */
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recv_size_y   = overlap->comm_data[ipulse].recv_size;

            sendptr = overlap->sendbuf + send_index0*local_fft_ndata[ZZ];
            recvptr = overlap->recvbuf;

            if (debug != NULL)
            {
                fprintf(debug, "PME fftgrid comm y %2d x %2d x %2d\n",
                        local_fft_ndata[XX], send_nindex, local_fft_ndata[ZZ]);
            }

#ifdef GMX_MPI
            MPI_Sendrecv(sendptr, send_size_y*datasize, GMX_MPI_REAL,
                         send_id, ipulse,
                         recvptr, recv_size_y*datasize, GMX_MPI_REAL,
                         recv_id, ipulse,
                         overlap->mpi_comm, &stat);
#endif

            for (x = 0; x < local_fft_ndata[XX]; x++)
            {
                for (y = 0; y < recv_nindex; y++)
                {
                    indg = (x*local_fft_size[YY] + y)*local_fft_size[ZZ];
                    indb = (x*recv_size_y        + y)*local_fft_ndata[ZZ];
                    for (z = 0; z < local_fft_ndata[ZZ]; z++)
                    {
                        fftgrid[indg+z] += recvptr[indb+z];
                    }
                }
            }

            if (pme->nnodes_major > 1)
            {
                /* Copy from the received buffer to the send buffer for dim 0 */
                sendptr = pme->overlap[0].sendbuf;
                for (x = 0; x < size_yx; x++)
                {
                    for (y = 0; y < recv_nindex; y++)
                    {
                        indg = (x*local_fft_ndata[YY] + y)*local_fft_ndata[ZZ];
                        indb = ((local_fft_ndata[XX] + x)*recv_size_y + y)*local_fft_ndata[ZZ];
                        for (z = 0; z < local_fft_ndata[ZZ]; z++)
                        {
                            sendptr[indg+z] += recvptr[indb+z];
                        }
                    }
                }
            }
        }
    }

    /* We only support a single pulse here.
     * This is not a severe limitation, as this code is only used
     * with OpenMP and with OpenMP the (PME) domains can be larger.
     */
    if (pme->nnodes_major > 1)
    {
        /* Major dimension */
        overlap = &pme->overlap[0];

        datasize = local_fft_ndata[YY]*local_fft_ndata[ZZ];
        gridsize = local_fft_size[YY] *local_fft_size[ZZ];

        ipulse = 0;

        send_id       = overlap->send_id[ipulse];
        recv_id       = overlap->recv_id[ipulse];
        send_nindex   = overlap->comm_data[ipulse].send_nindex;
        /* We don't use recv_index0, as we always receive starting at 0 */
        recv_nindex   = overlap->comm_data[ipulse].recv_nindex;

        sendptr = overlap->sendbuf;
        recvptr = overlap->recvbuf;

        if (debug != NULL)
        {
            fprintf(debug, "PME fftgrid comm x %2d x %2d x %2d\n",
                    send_nindex, local_fft_ndata[YY], local_fft_ndata[ZZ]);
        }

#ifdef GMX_MPI
        MPI_Sendrecv(sendptr, send_nindex*datasize, GMX_MPI_REAL,
                     send_id, ipulse,
                     recvptr, recv_nindex*datasize, GMX_MPI_REAL,
                     recv_id, ipulse,
                     overlap->mpi_comm, &stat);
#endif

        for (x = 0; x < recv_nindex; x++)
        {
            for (y = 0; y < local_fft_ndata[YY]; y++)
            {
                indg = (x*local_fft_size[YY]  + y)*local_fft_size[ZZ];
                indb = (x*local_fft_ndata[YY] + y)*local_fft_ndata[ZZ];
                for (z = 0; z < local_fft_ndata[ZZ]; z++)
                {
                    fftgrid[indg+z] += recvptr[indb+z];
                }
            }
        }
    }
}


static void spread_on_grid(gmx_pme_t pme,
                           pme_atomcomm_t *atc, pmegrids_t *grids,
                           gmx_bool bCalcSplines, gmx_bool bSpread,
                           real *fftgrid, gmx_bool bDoSplines, int grid_index)
{
    int nthread, thread;
#ifdef PME_TIME_THREADS
    gmx_cycles_t c1, c2, c3, ct1a, ct1b, ct1c;
    static double cs1     = 0, cs2 = 0, cs3 = 0;
    static double cs1a[6] = {0, 0, 0, 0, 0, 0};
    static int cnt        = 0;
#endif

    nthread = pme->nthread;
    assert(nthread > 0);

#ifdef PME_TIME_THREADS
    c1 = omp_cyc_start();
#endif
    if (bCalcSplines)
    {
#pragma omp parallel for num_threads(nthread) schedule(static)
        for (thread = 0; thread < nthread; thread++)
        {
            int start, end;

            start = atc->n* thread   /nthread;
            end   = atc->n*(thread+1)/nthread;

            /* Compute fftgrid index for all atoms,
             * with help of some extra variables.
             */
            calc_interpolation_idx(pme, atc, start, grid_index, end, thread);
        }
    }
#ifdef PME_TIME_THREADS
    c1   = omp_cyc_end(c1);
    cs1 += (double)c1;
#endif

#ifdef PME_TIME_THREADS
    c2 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        splinedata_t *spline;
        pmegrid_t *grid = NULL;

        /* make local bsplines  */
        if (grids == NULL || !pme->bUseThreads)
        {
            spline = &atc->spline[0];

            spline->n = atc->n;

            if (bSpread)
            {
                grid = &grids->grid;
            }
        }
        else
        {
            spline = &atc->spline[thread];

            if (grids->nthread == 1)
            {
                /* One thread, we operate on all coefficients */
                spline->n = atc->n;
            }
            else
            {
                /* Get the indices our thread should operate on */
                make_thread_local_ind(atc, thread, spline);
            }

            grid = &grids->grid_th[thread];
        }

        if (bCalcSplines)
        {
            make_bsplines(spline->theta, spline->dtheta, pme->pme_order,
                          atc->fractx, spline->n, spline->ind, atc->coefficient, bDoSplines);
        }

        if (bSpread)
        {
            /* put local atoms on grid. */
#ifdef PME_TIME_SPREAD
            ct1a = omp_cyc_start();
#endif
            spread_coefficients_bsplines_thread(grid, atc, spline, pme->spline_work);

            if (pme->bUseThreads)
            {
                copy_local_grid(pme, grids, grid_index, thread, fftgrid);
            }
#ifdef PME_TIME_SPREAD
            ct1a          = omp_cyc_end(ct1a);
            cs1a[thread] += (double)ct1a;
#endif
        }
    }
#ifdef PME_TIME_THREADS
    c2   = omp_cyc_end(c2);
    cs2 += (double)c2;
#endif

    if (bSpread && pme->bUseThreads)
    {
#ifdef PME_TIME_THREADS
        c3 = omp_cyc_start();
#endif
#pragma omp parallel for num_threads(grids->nthread) schedule(static)
        for (thread = 0; thread < grids->nthread; thread++)
        {
            reduce_threadgrid_overlap(pme, grids, thread,
                                      fftgrid,
                                      pme->overlap[0].sendbuf,
                                      pme->overlap[1].sendbuf,
                                      grid_index);
        }
#ifdef PME_TIME_THREADS
        c3   = omp_cyc_end(c3);
        cs3 += (double)c3;
#endif

        if (pme->nnodes > 1)
        {
            /* Communicate the overlapping part of the fftgrid.
             * For this communication call we need to check pme->bUseThreads
             * to have all ranks communicate here, regardless of pme->nthread.
             */
            sum_fftgrid_dd(pme, fftgrid, grid_index);
        }
    }

#ifdef PME_TIME_THREADS
    cnt++;
    if (cnt % 20 == 0)
    {
        printf("idx %.2f spread %.2f red %.2f",
               cs1*1e-9, cs2*1e-9, cs3*1e-9);
#ifdef PME_TIME_SPREAD
        for (thread = 0; thread < nthread; thread++)
        {
            printf(" %.2f", cs1a[thread]*1e-9);
        }
#endif
        printf("\n");
    }
#endif
}


static void dump_grid(FILE *fp,
                      int sx, int sy, int sz, int nx, int ny, int nz,
                      int my, int mz, const real *g)
{
    int x, y, z;

    for (x = 0; x < nx; x++)
    {
        for (y = 0; y < ny; y++)
        {
            for (z = 0; z < nz; z++)
            {
                fprintf(fp, "%2d %2d %2d %6.3f\n",
                        sx+x, sy+y, sz+z, g[(x*my + y)*mz + z]);
            }
        }
    }
}

static void dump_local_fftgrid(gmx_pme_t pme, const real *fftgrid)
{
    ivec local_fft_ndata, local_fft_offset, local_fft_size;

    gmx_parallel_3dfft_real_limits(pme->pfft_setup[PME_GRID_QA],
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    dump_grid(stderr,
              pme->pmegrid_start_ix,
              pme->pmegrid_start_iy,
              pme->pmegrid_start_iz,
              pme->pmegrid_nx-pme->pme_order+1,
              pme->pmegrid_ny-pme->pme_order+1,
              pme->pmegrid_nz-pme->pme_order+1,
              local_fft_size[YY],
              local_fft_size[ZZ],
              fftgrid);
}


void gmx_pme_calc_energy(gmx_pme_t pme, int n, rvec *x, real *q, real *V)
{
    pme_atomcomm_t *atc;
    pmegrids_t *grid;

    if (pme->nnodes > 1)
    {
        gmx_incons("gmx_pme_calc_energy called in parallel");
    }
    if (pme->bFEP_q > 1)
    {
        gmx_incons("gmx_pme_calc_energy with free energy");
    }

    atc            = &pme->atc_energy;
    atc->nthread   = 1;
    if (atc->spline == NULL)
    {
        snew(atc->spline, atc->nthread);
    }
    atc->nslab     = 1;
    atc->bSpread   = TRUE;
    atc->pme_order = pme->pme_order;
    atc->n         = n;
    pme_realloc_atomcomm_things(atc);
    atc->x           = x;
    atc->coefficient = q;

    /* We only use the A-charges grid */
    grid = &pme->pmegrid[PME_GRID_QA];

    /* Only calculate the spline coefficients, don't actually spread */
    spread_on_grid(pme, atc, NULL, TRUE, FALSE, pme->fftgrid[PME_GRID_QA], FALSE, PME_GRID_QA);

    *V = gather_energy_bsplines(pme, grid->grid.grid, atc);
}


static void reset_pmeonly_counters(gmx_wallcycle_t wcycle,
                                   gmx_walltime_accounting_t walltime_accounting,
                                   t_nrnb *nrnb, t_inputrec *ir,
                                   gmx_int64_t step)
{
    /* Reset all the counters related to performance over the run */
    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    init_nrnb(nrnb);
    if (ir->nsteps >= 0)
    {
        /* ir->nsteps is not used here, but we update it for consistency */
        ir->nsteps -= step - ir->init_step;
    }
    ir->init_step = step;
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_start(walltime_accounting);
}


static void gmx_pmeonly_switch(int *npmedata, gmx_pme_t **pmedata,
                               ivec grid_size,
                               t_commrec *cr, t_inputrec *ir,
                               gmx_pme_t *pme_ret)
{
    int ind;
    gmx_pme_t pme = NULL;

    ind = 0;
    while (ind < *npmedata)
    {
        pme = (*pmedata)[ind];
        if (pme->nkx == grid_size[XX] &&
            pme->nky == grid_size[YY] &&
            pme->nkz == grid_size[ZZ])
        {
            *pme_ret = pme;

            return;
        }

        ind++;
    }

    (*npmedata)++;
    srenew(*pmedata, *npmedata);

    /* Generate a new PME data structure, copying part of the old pointers */
    gmx_pme_reinit(&((*pmedata)[ind]), cr, pme, ir, grid_size);

    *pme_ret = (*pmedata)[ind];
}

int gmx_pmeonly(gmx_pme_t pme,
                t_commrec *cr,    t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                real ewaldcoeff_q, real ewaldcoeff_lj,
                t_inputrec *ir)
{
    int npmedata;
    gmx_pme_t *pmedata;
    gmx_pme_pp_t pme_pp;
    int  ret;
    int  natoms;
    matrix box;
    rvec *x_pp      = NULL, *f_pp = NULL;
    real *chargeA   = NULL, *chargeB = NULL;
    real *c6A       = NULL, *c6B = NULL;
    real *sigmaA    = NULL, *sigmaB = NULL;
    real lambda_q   = 0;
    real lambda_lj  = 0;
    int  maxshift_x = 0, maxshift_y = 0;
    real energy_q, energy_lj, dvdlambda_q, dvdlambda_lj;
    matrix vir_q, vir_lj;
    float cycles;
    int  count;
    gmx_bool bEnerVir;
    int pme_flags;
    gmx_int64_t step, step_rel;
    ivec grid_switch;

    /* This data will only use with PME tuning, i.e. switching PME grids */
    npmedata = 1;
    snew(pmedata, npmedata);
    pmedata[0] = pme;

    pme_pp = gmx_pme_pp_init(cr);

    init_nrnb(mynrnb);

    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* The reason for having a loop here is PME grid tuning/switching */
        do
        {
            /* Domain decomposition */
            ret = gmx_pme_recv_coeffs_coords(pme_pp,
                                             &natoms,
                                             &chargeA, &chargeB,
                                             &c6A, &c6B,
                                             &sigmaA, &sigmaB,
                                             box, &x_pp, &f_pp,
                                             &maxshift_x, &maxshift_y,
                                             &pme->bFEP_q, &pme->bFEP_lj,
                                             &lambda_q, &lambda_lj,
                                             &bEnerVir,
                                             &pme_flags,
                                             &step,
                                             grid_switch, &ewaldcoeff_q, &ewaldcoeff_lj);

            if (ret == pmerecvqxSWITCHGRID)
            {
                /* Switch the PME grid to grid_switch */
                gmx_pmeonly_switch(&npmedata, &pmedata, grid_switch, cr, ir, &pme);
            }

            if (ret == pmerecvqxRESETCOUNTERS)
            {
                /* Reset the cycle and flop counters */
                reset_pmeonly_counters(wcycle, walltime_accounting, mynrnb, ir, step);
            }
        }
        while (ret == pmerecvqxSWITCHGRID || ret == pmerecvqxRESETCOUNTERS);

        if (ret == pmerecvqxFINISH)
        {
            /* We should stop: break out of the loop */
            break;
        }

        step_rel = step - ir->init_step;

        if (count == 0)
        {
            wallcycle_start(wcycle, ewcRUN);
            walltime_accounting_start(walltime_accounting);
        }

        wallcycle_start(wcycle, ewcPMEMESH);

        dvdlambda_q  = 0;
        dvdlambda_lj = 0;
        clear_mat(vir_q);
        clear_mat(vir_lj);

        gmx_pme_do(pme, 0, natoms, x_pp, f_pp,
                   chargeA, chargeB, c6A, c6B, sigmaA, sigmaB, box,
                   cr, maxshift_x, maxshift_y, mynrnb, wcycle,
                   vir_q, ewaldcoeff_q, vir_lj, ewaldcoeff_lj,
                   &energy_q, &energy_lj, lambda_q, lambda_lj, &dvdlambda_q, &dvdlambda_lj,
                   pme_flags | GMX_PME_DO_ALL_F | (bEnerVir ? GMX_PME_CALC_ENER_VIR : 0));

        cycles = wallcycle_stop(wcycle, ewcPMEMESH);

        gmx_pme_send_force_vir_ener(pme_pp,
                                    f_pp, vir_q, energy_q, vir_lj, energy_lj,
                                    dvdlambda_q, dvdlambda_lj, cycles);

        count++;
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);

    walltime_accounting_end(walltime_accounting);

    return 0;
}

static void
calc_initial_lb_coeffs(gmx_pme_t pme, real *local_c6, real *local_sigma)
{
    int  i;

    for (i = 0; i < pme->atc[0].n; ++i)
    {
        real sigma4;

        sigma4                     = local_sigma[i];
        sigma4                     = sigma4*sigma4;
        sigma4                     = sigma4*sigma4;
        pme->atc[0].coefficient[i] = local_c6[i] / sigma4;
    }
}

static void
calc_next_lb_coeffs(gmx_pme_t pme, real *local_sigma)
{
    int  i;

    for (i = 0; i < pme->atc[0].n; ++i)
    {
        pme->atc[0].coefficient[i] *= local_sigma[i];
    }
}

static void
do_redist_pos_coeffs(gmx_pme_t pme, t_commrec *cr, int start, int homenr,
                     gmx_bool bFirst, rvec x[], real *data)
{
    int      d;
    pme_atomcomm_t *atc;
    atc = &pme->atc[0];

    for (d = pme->ndecompdim - 1; d >= 0; d--)
    {
        int             n_d;
        rvec           *x_d;
        real           *param_d;

        if (d == pme->ndecompdim - 1)
        {
            n_d     = homenr;
            x_d     = x + start;
            param_d = data;
        }
        else
        {
            n_d     = pme->atc[d + 1].n;
            x_d     = atc->x;
            param_d = atc->coefficient;
        }
        atc      = &pme->atc[d];
        atc->npd = n_d;
        if (atc->npd > atc->pd_nalloc)
        {
            atc->pd_nalloc = over_alloc_dd(atc->npd);
            srenew(atc->pd, atc->pd_nalloc);
        }
        pme_calc_pidx_wrapper(n_d, pme->recipbox, x_d, atc);
        where();
        /* Redistribute x (only once) and qA/c6A or qB/c6B */
        if (DOMAINDECOMP(cr))
        {
            dd_pmeredist_pos_coeffs(pme, n_d, bFirst, x_d, param_d, atc);
        }
    }
}

int gmx_pme_do(gmx_pme_t pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real *chargeA,   real *chargeB,
               real *c6A,       real *c6B,
               real *sigmaA,    real *sigmaB,
               matrix box, t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix vir_q,      real ewaldcoeff_q,
               matrix vir_lj,   real ewaldcoeff_lj,
               real *energy_q,  real *energy_lj,
               real lambda_q, real lambda_lj,
               real *dvdlambda_q, real *dvdlambda_lj,
               int flags)
{
    int     d, i, j, k, ntot, npme, grid_index, max_grid_index;
    int     nx, ny, nz;
    int     n_d, local_ny;
    pme_atomcomm_t *atc = NULL;
    pmegrids_t *pmegrid = NULL;
    real    *grid       = NULL;
    real    *ptr;
    rvec    *x_d, *f_d;
    real    *coefficient = NULL;
    real    energy_AB[4];
    matrix  vir_AB[4];
    real    scale, lambda;
    gmx_bool bClearF;
    gmx_parallel_3dfft_t pfft_setup;
    real *  fftgrid;
    t_complex * cfftgrid;
    int     thread;
    gmx_bool bFirst, bDoSplines;
    int fep_state;
    int fep_states_lj           = pme->bFEP_lj ? 2 : 1;
    const gmx_bool bCalcEnerVir = flags & GMX_PME_CALC_ENER_VIR;
    const gmx_bool bCalcF       = flags & GMX_PME_CALC_F;

    assert(pme->nnodes > 0);
    assert(pme->nnodes == 1 || pme->ndecompdim > 0);

    if (pme->nnodes > 1)
    {
        atc      = &pme->atc[0];
        atc->npd = homenr;
        if (atc->npd > atc->pd_nalloc)
        {
            atc->pd_nalloc = over_alloc_dd(atc->npd);
            srenew(atc->pd, atc->pd_nalloc);
        }
        for (d = pme->ndecompdim-1; d >= 0; d--)
        {
            atc           = &pme->atc[d];
            atc->maxshift = (atc->dimind == 0 ? maxshift_x : maxshift_y);
        }
    }
    else
    {
        atc = &pme->atc[0];
        /* This could be necessary for TPI */
        pme->atc[0].n = homenr;
        if (DOMAINDECOMP(cr))
        {
            pme_realloc_atomcomm_things(atc);
        }
        atc->x = x;
        atc->f = f;
    }

    m_inv_ur0(box, pme->recipbox);
    bFirst = TRUE;

    /* For simplicity, we construct the splines for all particles if
     * more than one PME calculations is needed. Some optimization
     * could be done by keeping track of which atoms have splines
     * constructed, and construct new splines on each pass for atoms
     * that don't yet have them.
     */

    bDoSplines = pme->bFEP || ((flags & GMX_PME_DO_COULOMB) && (flags & GMX_PME_DO_LJ));

    /* We need a maximum of four separate PME calculations:
     * grid_index=0: Coulomb PME with charges from state A
     * grid_index=1: Coulomb PME with charges from state B
     * grid_index=2: LJ PME with C6 from state A
     * grid_index=3: LJ PME with C6 from state B
     * For Lorentz-Berthelot combination rules, a separate loop is used to
     * calculate all the terms
     */

    /* If we are doing LJ-PME with LB, we only do Q here */
    max_grid_index = (pme->ljpme_combination_rule == eljpmeLB) ? DO_Q : DO_Q_AND_LJ;

    for (grid_index = 0; grid_index < max_grid_index; ++grid_index)
    {
        /* Check if we should do calculations at this grid_index
         * If grid_index is odd we should be doing FEP
         * If grid_index < 2 we should be doing electrostatic PME
         * If grid_index >= 2 we should be doing LJ-PME
         */
        if ((grid_index <  DO_Q && (!(flags & GMX_PME_DO_COULOMB) ||
                                    (grid_index == 1 && !pme->bFEP_q))) ||
            (grid_index >= DO_Q && (!(flags & GMX_PME_DO_LJ) ||
                                    (grid_index == 3 && !pme->bFEP_lj))))
        {
            continue;
        }
        /* Unpack structure */
        pmegrid    = &pme->pmegrid[grid_index];
        fftgrid    = pme->fftgrid[grid_index];
        cfftgrid   = pme->cfftgrid[grid_index];
        pfft_setup = pme->pfft_setup[grid_index];
        switch (grid_index)
        {
            case 0: coefficient = chargeA + start; break;
            case 1: coefficient = chargeB + start; break;
            case 2: coefficient = c6A + start; break;
            case 3: coefficient = c6B + start; break;
        }

        grid = pmegrid->grid.grid;

        if (debug)
        {
            fprintf(debug, "PME: number of ranks = %d, rank = %d\n",
                    cr->nnodes, cr->nodeid);
            fprintf(debug, "Grid = %p\n", (void*)grid);
            if (grid == NULL)
            {
                gmx_fatal(FARGS, "No grid!");
            }
        }
        where();

        if (pme->nnodes == 1)
        {
            atc->coefficient = coefficient;
        }
        else
        {
            wallcycle_start(wcycle, ewcPME_REDISTXF);
            do_redist_pos_coeffs(pme, cr, start, homenr, bFirst, x, coefficient);
            where();

            wallcycle_stop(wcycle, ewcPME_REDISTXF);
        }

        if (debug)
        {
            fprintf(debug, "Rank= %6d, pme local particles=%6d\n",
                    cr->nodeid, atc->n);
        }

        if (flags & GMX_PME_SPREAD)
        {
            wallcycle_start(wcycle, ewcPME_SPREADGATHER);

            /* Spread the coefficients on a grid */
            spread_on_grid(pme, &pme->atc[0], pmegrid, bFirst, TRUE, fftgrid, bDoSplines, grid_index);

            if (bFirst)
            {
                inc_nrnb(nrnb, eNR_WEIGHTS, DIM*atc->n);
            }
            inc_nrnb(nrnb, eNR_SPREADBSP,
                     pme->pme_order*pme->pme_order*pme->pme_order*atc->n);

            if (!pme->bUseThreads)
            {
                wrap_periodic_pmegrid(pme, grid);

                /* sum contributions to local grid from other nodes */
#ifdef GMX_MPI
                if (pme->nnodes > 1)
                {
                    gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_FORWARD);
                    where();
                }
#endif

                copy_pmegrid_to_fftgrid(pme, grid, fftgrid, grid_index);
            }

            wallcycle_stop(wcycle, ewcPME_SPREADGATHER);

            /*
               dump_local_fftgrid(pme,fftgrid);
               exit(0);
             */
        }

        /* Here we start a large thread parallel region */
#pragma omp parallel num_threads(pme->nthread) private(thread)
        {
            thread = gmx_omp_get_thread_num();
            if (flags & GMX_PME_SOLVE)
            {
                int loop_count;

                /* do 3d-fft */
                if (thread == 0)
                {
                    wallcycle_start(wcycle, ewcPME_FFT);
                }
                gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_REAL_TO_COMPLEX,
                                           thread, wcycle);
                if (thread == 0)
                {
                    wallcycle_stop(wcycle, ewcPME_FFT);
                }
                where();

                /* solve in k-space for our local cells */
                if (thread == 0)
                {
                    wallcycle_start(wcycle, (grid_index < DO_Q ? ewcPME_SOLVE : ewcLJPME));
                }
                if (grid_index < DO_Q)
                {
                    loop_count =
                        solve_pme_yzx(pme, cfftgrid, ewaldcoeff_q,
                                      box[XX][XX]*box[YY][YY]*box[ZZ][ZZ],
                                      bCalcEnerVir,
                                      pme->nthread, thread);
                }
                else
                {
                    loop_count =
                        solve_pme_lj_yzx(pme, &cfftgrid, FALSE, ewaldcoeff_lj,
                                         box[XX][XX]*box[YY][YY]*box[ZZ][ZZ],
                                         bCalcEnerVir,
                                         pme->nthread, thread);
                }

                if (thread == 0)
                {
                    wallcycle_stop(wcycle, (grid_index < DO_Q ? ewcPME_SOLVE : ewcLJPME));
                    where();
                    inc_nrnb(nrnb, eNR_SOLVEPME, loop_count);
                }
            }

            if (bCalcF)
            {
                /* do 3d-invfft */
                if (thread == 0)
                {
                    where();
                    wallcycle_start(wcycle, ewcPME_FFT);
                }
                gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_COMPLEX_TO_REAL,
                                           thread, wcycle);
                if (thread == 0)
                {
                    wallcycle_stop(wcycle, ewcPME_FFT);

                    where();

                    if (pme->nodeid == 0)
                    {
                        ntot  = pme->nkx*pme->nky*pme->nkz;
                        npme  = ntot*log((real)ntot)/log(2.0);
                        inc_nrnb(nrnb, eNR_FFT, 2*npme);
                    }

                    /* Note: this wallcycle region is closed below
                       outside an OpenMP region, so take care if
                       refactoring code here. */
                    wallcycle_start(wcycle, ewcPME_SPREADGATHER);
                }

                copy_fftgrid_to_pmegrid(pme, fftgrid, grid, grid_index, pme->nthread, thread);
            }
        }
        /* End of thread parallel section.
         * With MPI we have to synchronize here before gmx_sum_qgrid_dd.
         */

        if (bCalcF)
        {
            /* distribute local grid to all nodes */
#ifdef GMX_MPI
            if (pme->nnodes > 1)
            {
                gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_BACKWARD);
            }
#endif
            where();

            unwrap_periodic_pmegrid(pme, grid);

            /* interpolate forces for our local atoms */

            where();

            /* If we are running without parallelization,
             * atc->f is the actual force array, not a buffer,
             * therefore we should not clear it.
             */
            lambda  = grid_index < DO_Q ? lambda_q : lambda_lj;
            bClearF = (bFirst && PAR(cr));
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
            for (thread = 0; thread < pme->nthread; thread++)
            {
                gather_f_bsplines(pme, grid, bClearF, atc,
                                  &atc->spline[thread],
                                  pme->bFEP ? (grid_index % 2 == 0 ? 1.0-lambda : lambda) : 1.0);
            }

            where();

            inc_nrnb(nrnb, eNR_GATHERFBSP,
                     pme->pme_order*pme->pme_order*pme->pme_order*pme->atc[0].n);
            /* Note: this wallcycle region is opened above inside an OpenMP
               region, so take care if refactoring code here. */
            wallcycle_stop(wcycle, ewcPME_SPREADGATHER);
        }

        if (bCalcEnerVir)
        {
            /* This should only be called on the master thread
             * and after the threads have synchronized.
             */
            if (grid_index < 2)
            {
                get_pme_ener_vir_q(pme, pme->nthread, &energy_AB[grid_index], vir_AB[grid_index]);
            }
            else
            {
                get_pme_ener_vir_lj(pme, pme->nthread, &energy_AB[grid_index], vir_AB[grid_index]);
            }
        }
        bFirst = FALSE;
    } /* of grid_index-loop */

    /* For Lorentz-Berthelot combination rules in LJ-PME, we need to calculate
     * seven terms. */

    if ((flags & GMX_PME_DO_LJ) && pme->ljpme_combination_rule == eljpmeLB)
    {
        /* Loop over A- and B-state if we are doing FEP */
        for (fep_state = 0; fep_state < fep_states_lj; ++fep_state)
        {
            real *local_c6 = NULL, *local_sigma = NULL, *RedistC6 = NULL, *RedistSigma = NULL;
            if (pme->nnodes == 1)
            {
                if (pme->lb_buf1 == NULL)
                {
                    pme->lb_buf_nalloc = pme->atc[0].n;
                    snew(pme->lb_buf1, pme->lb_buf_nalloc);
                }
                pme->atc[0].coefficient = pme->lb_buf1;
                switch (fep_state)
                {
                    case 0:
                        local_c6      = c6A;
                        local_sigma   = sigmaA;
                        break;
                    case 1:
                        local_c6      = c6B;
                        local_sigma   = sigmaB;
                        break;
                    default:
                        gmx_incons("Trying to access wrong FEP-state in LJ-PME routine");
                }
            }
            else
            {
                atc = &pme->atc[0];
                switch (fep_state)
                {
                    case 0:
                        RedistC6      = c6A;
                        RedistSigma   = sigmaA;
                        break;
                    case 1:
                        RedistC6      = c6B;
                        RedistSigma   = sigmaB;
                        break;
                    default:
                        gmx_incons("Trying to access wrong FEP-state in LJ-PME routine");
                }
                wallcycle_start(wcycle, ewcPME_REDISTXF);

                do_redist_pos_coeffs(pme, cr, start, homenr, bFirst, x, RedistC6);
                if (pme->lb_buf_nalloc < atc->n)
                {
                    pme->lb_buf_nalloc = atc->nalloc;
                    srenew(pme->lb_buf1, pme->lb_buf_nalloc);
                    srenew(pme->lb_buf2, pme->lb_buf_nalloc);
                }
                local_c6 = pme->lb_buf1;
                for (i = 0; i < atc->n; ++i)
                {
                    local_c6[i] = atc->coefficient[i];
                }
                where();

                do_redist_pos_coeffs(pme, cr, start, homenr, FALSE, x, RedistSigma);
                local_sigma = pme->lb_buf2;
                for (i = 0; i < atc->n; ++i)
                {
                    local_sigma[i] = atc->coefficient[i];
                }
                where();

                wallcycle_stop(wcycle, ewcPME_REDISTXF);
            }
            calc_initial_lb_coeffs(pme, local_c6, local_sigma);

            /*Seven terms in LJ-PME with LB, grid_index < 2 reserved for electrostatics*/
            for (grid_index = 2; grid_index < 9; ++grid_index)
            {
                /* Unpack structure */
                pmegrid    = &pme->pmegrid[grid_index];
                fftgrid    = pme->fftgrid[grid_index];
                cfftgrid   = pme->cfftgrid[grid_index];
                pfft_setup = pme->pfft_setup[grid_index];
                calc_next_lb_coeffs(pme, local_sigma);
                grid = pmegrid->grid.grid;
                where();

                if (flags & GMX_PME_SPREAD)
                {
                    wallcycle_start(wcycle, ewcPME_SPREADGATHER);
                    /* Spread the c6 on a grid */
                    spread_on_grid(pme, &pme->atc[0], pmegrid, bFirst, TRUE, fftgrid, bDoSplines, grid_index);

                    if (bFirst)
                    {
                        inc_nrnb(nrnb, eNR_WEIGHTS, DIM*atc->n);
                    }

                    inc_nrnb(nrnb, eNR_SPREADBSP,
                             pme->pme_order*pme->pme_order*pme->pme_order*atc->n);
                    if (pme->nthread == 1)
                    {
                        wrap_periodic_pmegrid(pme, grid);
                        /* sum contributions to local grid from other nodes */
#ifdef GMX_MPI
                        if (pme->nnodes > 1)
                        {
                            gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_FORWARD);
                            where();
                        }
#endif
                        copy_pmegrid_to_fftgrid(pme, grid, fftgrid, grid_index);
                    }
                    wallcycle_stop(wcycle, ewcPME_SPREADGATHER);
                }
                /*Here we start a large thread parallel region*/
#pragma omp parallel num_threads(pme->nthread) private(thread)
                {
                    thread = gmx_omp_get_thread_num();
                    if (flags & GMX_PME_SOLVE)
                    {
                        /* do 3d-fft */
                        if (thread == 0)
                        {
                            wallcycle_start(wcycle, ewcPME_FFT);
                        }

                        gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_REAL_TO_COMPLEX,
                                                   thread, wcycle);
                        if (thread == 0)
                        {
                            wallcycle_stop(wcycle, ewcPME_FFT);
                        }
                        where();
                    }
                }
                bFirst = FALSE;
            }
            if (flags & GMX_PME_SOLVE)
            {
                /* solve in k-space for our local cells */
#pragma omp parallel num_threads(pme->nthread) private(thread)
                {
                    int loop_count;
                    thread = gmx_omp_get_thread_num();
                    if (thread == 0)
                    {
                        wallcycle_start(wcycle, ewcLJPME);
                    }

                    loop_count =
                        solve_pme_lj_yzx(pme, &pme->cfftgrid[2], TRUE, ewaldcoeff_lj,
                                         box[XX][XX]*box[YY][YY]*box[ZZ][ZZ],
                                         bCalcEnerVir,
                                         pme->nthread, thread);
                    if (thread == 0)
                    {
                        wallcycle_stop(wcycle, ewcLJPME);
                        where();
                        inc_nrnb(nrnb, eNR_SOLVEPME, loop_count);
                    }
                }
            }

            if (bCalcEnerVir)
            {
                /* This should only be called on the master thread and
                 * after the threads have synchronized.
                 */
                get_pme_ener_vir_lj(pme, pme->nthread, &energy_AB[2+fep_state], vir_AB[2+fep_state]);
            }

            if (bCalcF)
            {
                bFirst = !(flags & GMX_PME_DO_COULOMB);
                calc_initial_lb_coeffs(pme, local_c6, local_sigma);
                for (grid_index = 8; grid_index >= 2; --grid_index)
                {
                    /* Unpack structure */
                    pmegrid    = &pme->pmegrid[grid_index];
                    fftgrid    = pme->fftgrid[grid_index];
                    cfftgrid   = pme->cfftgrid[grid_index];
                    pfft_setup = pme->pfft_setup[grid_index];
                    grid       = pmegrid->grid.grid;
                    calc_next_lb_coeffs(pme, local_sigma);
                    where();
#pragma omp parallel num_threads(pme->nthread) private(thread)
                    {
                        thread = gmx_omp_get_thread_num();
                        /* do 3d-invfft */
                        if (thread == 0)
                        {
                            where();
                            wallcycle_start(wcycle, ewcPME_FFT);
                        }

                        gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_COMPLEX_TO_REAL,
                                                   thread, wcycle);
                        if (thread == 0)
                        {
                            wallcycle_stop(wcycle, ewcPME_FFT);

                            where();

                            if (pme->nodeid == 0)
                            {
                                ntot  = pme->nkx*pme->nky*pme->nkz;
                                npme  = ntot*log((real)ntot)/log(2.0);
                                inc_nrnb(nrnb, eNR_FFT, 2*npme);
                            }
                            wallcycle_start(wcycle, ewcPME_SPREADGATHER);
                        }

                        copy_fftgrid_to_pmegrid(pme, fftgrid, grid, grid_index, pme->nthread, thread);

                    } /*#pragma omp parallel*/

                    /* distribute local grid to all nodes */
#ifdef GMX_MPI
                    if (pme->nnodes > 1)
                    {
                        gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_BACKWARD);
                    }
#endif
                    where();

                    unwrap_periodic_pmegrid(pme, grid);

                    /* interpolate forces for our local atoms */
                    where();
                    bClearF = (bFirst && PAR(cr));
                    scale   = pme->bFEP ? (fep_state < 1 ? 1.0-lambda_lj : lambda_lj) : 1.0;
                    scale  *= lb_scale_factor[grid_index-2];
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                    for (thread = 0; thread < pme->nthread; thread++)
                    {
                        gather_f_bsplines(pme, grid, bClearF, &pme->atc[0],
                                          &pme->atc[0].spline[thread],
                                          scale);
                    }
                    where();

                    inc_nrnb(nrnb, eNR_GATHERFBSP,
                             pme->pme_order*pme->pme_order*pme->pme_order*pme->atc[0].n);
                    wallcycle_stop(wcycle, ewcPME_SPREADGATHER);

                    bFirst = FALSE;
                } /* for (grid_index = 8; grid_index >= 2; --grid_index) */
            }     /* if (bCalcF) */
        }         /* for (fep_state = 0; fep_state < fep_states_lj; ++fep_state) */
    }             /* if ((flags & GMX_PME_DO_LJ) && pme->ljpme_combination_rule == eljpmeLB) */

    if (bCalcF && pme->nnodes > 1)
    {
        wallcycle_start(wcycle, ewcPME_REDISTXF);
        for (d = 0; d < pme->ndecompdim; d++)
        {
            atc = &pme->atc[d];
            if (d == pme->ndecompdim - 1)
            {
                n_d = homenr;
                f_d = f + start;
            }
            else
            {
                n_d = pme->atc[d+1].n;
                f_d = pme->atc[d+1].f;
            }
            if (DOMAINDECOMP(cr))
            {
                dd_pmeredist_f(pme, atc, n_d, f_d,
                               d == pme->ndecompdim-1 && pme->bPPnode);
            }
        }

        wallcycle_stop(wcycle, ewcPME_REDISTXF);
    }
    where();

    if (bCalcEnerVir)
    {
        if (flags & GMX_PME_DO_COULOMB)
        {
            if (!pme->bFEP_q)
            {
                *energy_q = energy_AB[0];
                m_add(vir_q, vir_AB[0], vir_q);
            }
            else
            {
                *energy_q       = (1.0-lambda_q)*energy_AB[0] + lambda_q*energy_AB[1];
                *dvdlambda_q   += energy_AB[1] - energy_AB[0];
                for (i = 0; i < DIM; i++)
                {
                    for (j = 0; j < DIM; j++)
                    {
                        vir_q[i][j] += (1.0-lambda_q)*vir_AB[0][i][j] +
                            lambda_q*vir_AB[1][i][j];
                    }
                }
            }
            if (debug)
            {
                fprintf(debug, "Electrostatic PME mesh energy: %g\n", *energy_q);
            }
        }
        else
        {
            *energy_q = 0;
        }

        if (flags & GMX_PME_DO_LJ)
        {
            if (!pme->bFEP_lj)
            {
                *energy_lj = energy_AB[2];
                m_add(vir_lj, vir_AB[2], vir_lj);
            }
            else
            {
                *energy_lj     = (1.0-lambda_lj)*energy_AB[2] + lambda_lj*energy_AB[3];
                *dvdlambda_lj += energy_AB[3] - energy_AB[2];
                for (i = 0; i < DIM; i++)
                {
                    for (j = 0; j < DIM; j++)
                    {
                        vir_lj[i][j] += (1.0-lambda_lj)*vir_AB[2][i][j] + lambda_lj*vir_AB[3][i][j];
                    }
                }
            }
            if (debug)
            {
                fprintf(debug, "Lennard-Jones PME mesh energy: %g\n", *energy_lj);
            }
        }
        else
        {
            *energy_lj = 0;
        }
    }
    return 0;
}
