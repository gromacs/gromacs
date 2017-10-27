/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

/* TODO This file is a temporary holding area for stuff local to the
 * PME code, before it acquires some more normal ewald/file.c and
 * ewald/file.h structure.  In future clean up, get rid of this file,
 * to build more normal. */

#ifndef GMX_EWALD_PME_INTERNAL_H
#define GMX_EWALD_PME_INTERNAL_H

#include "config.h"

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"

#include "pme-gpu-types.h"

//! A repeat of typedef from parallel_3dfft.h
typedef struct gmx_parallel_3dfft *gmx_parallel_3dfft_t;

struct t_commrec;
struct t_inputrec;
struct PmeGpu;

//@{
//! Grid indices for A state for charge and Lennard-Jones C6
#define PME_GRID_QA    0
#define PME_GRID_C6A   2
//@}

//@{
/*! \brief Flags that indicate the number of PME grids in use */
#define DO_Q           2 /* Electrostatic grids have index q<2 */
#define DO_Q_AND_LJ    4 /* non-LB LJ grids have index 2 <= q < 4 */
#define DO_Q_AND_LJ_LB 9 /* With LB rules we need a total of 2+7 grids */
//@}

/*! \brief Pascal triangle coefficients scaled with (1/2)^6 for LJ-PME with LB-rules */
static const real lb_scale_factor[] = {
    1.0/64, 6.0/64, 15.0/64, 20.0/64,
    15.0/64, 6.0/64, 1.0/64
};

/*! \brief Pascal triangle coefficients used in solve_pme_lj_yzx, only need to do 4 calculations due to symmetry */
static const real lb_scale_factor_symm[] = { 2.0/64, 12.0/64, 30.0/64, 20.0/64 };

/*! \brief We only define a maximum to be able to use local arrays without allocation.
 * An order larger than 12 should never be needed, even for test cases.
 * If needed it can be changed here.
 */
#define PME_ORDER_MAX 12

/*! \brief As gmx_pme_init, but takes most settings, except the grid/Ewald coefficients, from pme_src.
 * This is only called when the PME cut-off/grid size changes.
 */
void gmx_pme_reinit(struct gmx_pme_t **pmedata,
                    t_commrec *        cr,
                    struct gmx_pme_t * pme_src,
                    const t_inputrec * ir,
                    const ivec         grid_size,
                    real               ewaldcoeff_q,
                    real               ewaldcoeff_lj);


/* Temporary suppression until these structs become opaque and don't live in
 * a header that is included by other headers. Also, until then I have no
 * idea what some of the names mean. */

//! @cond Doxygen_Suppress

/*! \brief Data structure for grid communication */
typedef struct {
    int send_index0;
    int send_nindex;
    int recv_index0;
    int recv_nindex;
    int recv_size;   /* Receive buffer width, used with OpenMP */
} pme_grid_comm_t;

/*! \brief Data structure for grid overlap communication */
typedef struct {
#if GMX_MPI
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

/*! \brief Data structure for organizing particle allocation to threads */
typedef struct {
    int *n;      /* Cumulative counts of the number of particles per thread */
    int  nalloc; /* Allocation size of i */
    int *i;      /* Particle indices ordered on thread index (n) */
} thread_plist_t;

/*! \brief Helper typedef for spline vectors */
typedef real *splinevec[DIM];

/*! \brief Data structure for beta-spline interpolation */
typedef struct {
    int       n;
    int      *ind;
    splinevec theta;
    real     *ptr_theta_z;
    splinevec dtheta;
    real     *ptr_dtheta_z;
} splinedata_t;

/*! \brief Data structure for coordinating transfer between PP and PME ranks*/
struct pme_atomcomm_t{
    int      dimind;        /* The index of the dimension, 0=x, 1=y */
    int      nslab;
    int      nodeid;
#if GMX_MPI
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
};

/*! \brief Data structure for a single PME grid */
struct pmegrid_t{
    ivec  ci;     /* The spatial location of this grid         */
    ivec  n;      /* The used size of *grid, including order-1 */
    ivec  offset; /* The grid offset from the full node grid   */
    int   order;  /* PME spreading order                       */
    ivec  s;      /* The allocated size of *grid, s >= n       */
    real *grid;   /* The grid local thread, size n             */
};

/*! \brief Data structures for PME grids */
struct pmegrids_t{
    pmegrid_t  grid;         /* The full node grid (non thread-local)            */
    int        nthread;      /* The number of threads operating on this grid     */
    ivec       nc;           /* The local spatial decomposition over the threads */
    pmegrid_t *grid_th;      /* Array of grids for each thread                   */
    real      *grid_all;     /* Allocated array for the grids in *grid_th        */
    int       *g2t[DIM];     /* The grid to thread index                         */
    ivec       nthread_comm; /* The number of threads to communicate with        */
};

/*! \brief Data structure for spline-interpolation working buffers */
struct pme_spline_work;

/*! \brief Data structure for working buffers */
struct pme_solve_work_t;

/*! \brief Master PME data structure */
struct gmx_pme_t {
    int           ndecompdim; /* The number of decomposition dimensions */
    int           nodeid;     /* Our nodeid in mpi->mpi_comm */
    int           nodeid_major;
    int           nodeid_minor;
    int           nnodes;    /* The number of nodes doing PME */
    int           nnodes_major;
    int           nnodes_minor;

    MPI_Comm      mpi_comm;
    MPI_Comm      mpi_comm_d[2]; /* Indexed on dimension, 0=x, 1=y */
#if GMX_MPI
    MPI_Datatype  rvec_mpi;      /* the pme vector's MPI type */
#endif

    gmx_bool   bUseThreads;   /* Does any of the PME ranks have nthread>1 ?  */
    int        nthread;       /* The number of threads doing PME on our rank */

    gmx_bool   bPPnode;       /* Node also does particle-particle forces */
    bool       doCoulomb;     /* Apply PME to electrostatics */
    bool       doLJ;          /* Apply PME to Lennard-Jones r^-6 interactions */
    gmx_bool   bFEP;          /* Compute Free energy contribution */
    gmx_bool   bFEP_q;
    gmx_bool   bFEP_lj;
    int        nkx, nky, nkz; /* Grid dimensions */
    gmx_bool   bP3M;          /* Do P3M: optimize the influence function */
    int        pme_order;
    real       ewaldcoeff_q;  /* Ewald splitting coefficient for Coulomb */
    real       ewaldcoeff_lj; /* Ewald splitting coefficient for r^-6 */
    real       epsilon_r;


    enum PmeRunMode runMode; /* Which codepath is the PME runner taking - CPU, GPU, mixed;
                              * TODO: this is the information that should be owned by the task scheduler,
                              * and ideally not be duplicated here.
                              */

    PmeGpu      *gpu;        /* A pointer to the GPU data.
                              * TODO: this should be unique or a shared pointer.
                              * Currently in practice there is a single gmx_pme_t instance while a code
                              * is partially set up for many of them. The PME tuning calls gmx_pme_reinit()
                              * which fully reinitializes the one and only PME structure anew while maybe
                              * keeping the old grid buffers if they were already large enough.
                              * This small choice should be made clear in the later refactoring -
                              * do we store many PME objects for different grid sizes,
                              * or a single PME object that handles different grid sizes gracefully.
                              */


    class EwaldBoxZScaler *boxScaler;   /**< The scaling data Ewald uses with walls (set at pme_init constant for the entire run) */


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
    pme_spline_work          *spline_work;

    real                    **fftgrid; /* Grids for FFT. With 1D FFT decomposition this can be a pointer */
    /* inside the interpolation grid, but separate for 2D PME decomp. */
    int                       fftgrid_nx, fftgrid_ny, fftgrid_nz;

    t_complex               **cfftgrid; /* Grids for complex FFT data */

    int                       cfftgrid_nx, cfftgrid_ny, cfftgrid_nz;

    gmx_parallel_3dfft_t     *pfft_setup;

    int                      *nnx, *nny, *nnz;
    real                     *fshx, *fshy, *fshz;

    pme_atomcomm_t            atc[2]; /* Indexed on decomposition index */
    matrix                    recipbox;
    real                      boxVolume;
    splinevec                 bsp_mod;
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
    struct pme_solve_work_t *solve_work;

    /* Work data for sum_qgrid */
    real *   sum_qgrid_tmp;
    real *   sum_qgrid_dd_tmp;
};

//! @endcond

/*! \brief
 * Finds out if PME is currently running on GPU.
 * TODO: should this be removed eventually?
 *
 * \param[in] pme  The PME structure.
 * \returns        True if PME runs on GPU currently, false otherwise.
 */
inline bool pme_gpu_active(const gmx_pme_t *pme)
{
    return (pme != nullptr) && (pme->runMode != PmeRunMode::CPU);
}

/*! \brief Tell our PME-only node to switch to a new grid size */
void gmx_pme_send_switchgrid(t_commrec *cr,
                             ivec       grid_size,
                             real       ewaldcoeff_q,
                             real       ewaldcoeff_lj);

#endif
