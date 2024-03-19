/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include <memory>
#include <vector>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/defaultinitializationallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/unique_cptr.h"

//! A repeat of typedef from parallel_3dfft.h
typedef struct gmx_parallel_3dfft* gmx_parallel_3dfft_t;

struct t_commrec;
struct t_inputrec;
struct PmeGpu;
class EwaldBoxZScaler;
enum class PmeRunMode;
enum class LongRangeVdW : int;
class PmeSolve;

//! The number of grids for LJ-PME with LB combination rules
static constexpr int sc_numGridsLJLB = 7;

/*! \brief Pascal triangle coefficients scaled with (1/2)^6 for LJ-PME with LB-rules */
static const real lb_scale_factor[] = { 1.0 / 64,  6.0 / 64, 15.0 / 64, 20.0 / 64,
                                        15.0 / 64, 6.0 / 64, 1.0 / 64 };

/*! \brief Pascal triangle coefficients used in solve_pme_lj_yzx, only need to do 4 calculations due to symmetry */
static const real lb_scale_factor_symm[] = { 2.0 / 64, 12.0 / 64, 30.0 / 64, 20.0 / 64 };

/*! \brief We only define a maximum to be able to use local arrays without allocation.
 * An order larger than 12 should never be needed, even for test cases.
 * If needed it can be changed here.
 */
#define PME_ORDER_MAX 12


/* Temporary suppression until these structs become opaque and don't live in
 * a header that is included by other headers. Also, until then I have no
 * idea what some of the names mean. */

//! @cond Doxygen_Suppress

template<typename T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;

template<typename T>
using FastVector = std::vector<T, gmx::DefaultInitializationAllocator<T>>;

// Storage for all PME coefficient grids (not the FFT grids)
struct PmeGridsStorage
{
    /* Storage for Coulomb coefficient grids.
     * The major vector is for A and B charges, B is only present with perturbation.
     * The middle vector has either size 1 grid (with single OpenMP thread) or 1 + #threads entries.
     * The minor vector is for a grid, aligned for in case aligned SIMD load/stores are used.
     */
    std::vector<std::vector<AlignedVector<real>>> coulomb;
    /* Storage for LJ coefficient grids.
     * The major vector has size 1 with geometric combination rules and 6 with LB combination rules.
     * The middle vector has either size 1 grid (with single OpenMP thread) or 1 + #threads entries.
     * The minor vector is for a grid, aligned for in case aligned SIMD load/stores are used.
     */
    std::vector<std::vector<AlignedVector<real>>> lj;
};

/*! \brief Data structure for grid communication */
struct pme_grid_comm_t
{
    int send_id; //!< Source rank id
    int send_index0;
    int send_nindex;
    int recv_id; //!< Destination rank id
    int recv_index0;
    int recv_nindex;
    int recv_size = 0; //!< Receive buffer width, used with OpenMP
};

/*! \brief Data structure for grid overlap communication in a single dimension */
struct pme_overlap_t
{
    MPI_Comm                     mpi_comm;  //!< MPI communcator
    int                          nnodes;    //!< Number of ranks
    int                          nodeid;    //!< Unique rank identifcator
    std::vector<int>             s2g0;      //!< The local interpolation grid start
    std::vector<int>             s2g1;      //!< The local interpolation grid end
    int                          send_size; //!< Send buffer width, used with OpenMP
    std::vector<pme_grid_comm_t> comm_data; //!< All the individual communication data for each rank
    std::vector<real>            sendbuf;   //!< Shared buffer for sending
    std::vector<real>            recvbuf;   //!< Shared buffer for receiving
};

/*! \brief Data structure for organizing particle allocation to threads */
struct AtomToThreadMap
{
    //! Cumulative counts of the number of particles per thread
    int* n = nullptr;
    //! Storage buffer for n
    std::vector<int> nBuffer;
    //! Particle indices ordered on thread index (n)
    FastVector<int> i;
};

/*! \internal
 * \brief Coefficients for theta or dtheta
 */
class SplineCoefficients
{
public:
    //! Reallocate for use with up to nalloc coefficients
    void realloc(int nalloc);

    //! Pointers to the coefficient buffer for x, y, z
    std::array<real*, DIM> coefficients = { { nullptr } };

private:
    //! Storage for x coefficients
    std::vector<real> bufferX_;
    //! Storage for y coefficients
    std::vector<real> bufferY_;
    //! Storage for z coefficients, aligned for SIMD load
    AlignedVector<real> bufferZ_;
};

/*! \brief Data structure for beta-spline interpolation */
struct splinedata_t
{
    int                n = 0;
    FastVector<int>    ind;
    SplineCoefficients theta;
    SplineCoefficients dtheta;
    int                nalloc = 0;
};

/*! \brief PME slab MPI communication setup */
struct SlabCommSetup
{
    //! The nodes to send x and q to with DD
    int node_dest;
    //! The nodes to receive x and q from with DD
    int node_src;
    //! The number of atoms to receive
    int rcount;
};

/*! \internal
 * \brief Data structure for coordinating transfers between PME ranks along one dimension
 *
 * Also used for passing coordinates, coefficients and forces to and from PME routines.
 */
class PmeAtomComm
{
public:
    //! Constructor, \p PmeMpiCommunicator is the communicator for this dimension
    PmeAtomComm(MPI_Comm PmeMpiCommunicator, int numThreads, int pmeOrder, int dimIndex, bool doSpread);

    //! Set the atom count and when necessary resizes atom buffers
    void setNumAtoms(int numAtoms);

    //! Returns the atom count
    int numAtoms() const { return numAtoms_; }

    //! Returns the number of atoms to send to each rank
    gmx::ArrayRef<int> sendCount()
    {
        GMX_ASSERT(!count_thread.empty(), "Need at least one thread_count");
        return count_thread[0];
    }

    //! The index of the dimension, 0=x, 1=y
    int dimind = 0;
    //! The number of slabs and ranks this dimension is decomposed over
    int nslab = 1;
    //! Our slab index, our MPI rank is equal to this number
    int slabIndex = 0;
    //! Communicator for this dimension
    MPI_Comm mpi_comm;

    //! Communication setup for each slab, ordered as alternating forward/backward and increasing slab shift
    std::vector<SlabCommSetup> slabCommSetup;
    //! The maximum communication distance counted in MPI ranks
    int maxshift = 0;
    //! Working buffer indices, indexed with the slab index
    std::vector<int> bufferIndices;

    //! The target slab index for each particle
    FastVector<int> pd;
    //! Target particle counts for each slab, for each thread
    std::vector<std::vector<int>> count_thread;

private:
    //! The number of atoms
    int numAtoms_ = 0;

public:
    //! The coordinates
    gmx::ArrayRef<const gmx::RVec> x;
    //! The coefficient, charges or LJ C6
    gmx::ArrayRef<const real> coefficient;
    //! The forces
    gmx::ArrayRef<gmx::RVec> f;
    //! Coordinate buffer, used only with nslab > 1
    FastVector<gmx::RVec> xBuffer;
    //! Coefficient buffer, used only with nslab > 1
    FastVector<real> coefficientBuffer;
    //! Force buffer, used only with nslab > 1
    FastVector<gmx::RVec> fBuffer;
    //! Tells whether these coordinates are used for spreading
    bool bSpread;
    //! The PME order
    int pme_order;
    //! The grid index per atom
    FastVector<gmx::IVec> idx;
    //! Fractional atom coordinates relative to the lower cell boundary
    FastVector<gmx::RVec> fractx;

    //! The number of threads to use in PME
    int nthread;
    //! Thread index for each atom
    FastVector<int>              thread_idx;
    std::vector<AtomToThreadMap> threadMap;
    std::vector<splinedata_t>    spline;
};

/*! \brief Data structure for a single PME grid */
struct pmegrid_t
{
    ivec ci;     /* The spatial location of this grid         */
    ivec n;      /* The used size of *grid, including order-1 */
    ivec offset; /* The grid offset from the full node grid   */
    int  order;  /* PME spreading order                       */
    ivec s;      /* The allocated size of *grid, s >= n       */
    // The local grid, used size n. Data owned by gmx_pme_t::pmeGridsStorage.
    gmx::ArrayRef<real> grid;
};

/*! \brief Data structures for PME grids */
struct pmegrids_t
{
    pmegrid_t              grid;    /* The full node grid (non thread-local)            */
    int                    nthread; /* The number of threads operating on this grid     */
    ivec                   nc;      /* The local spatial decomposition over the threads */
    std::vector<pmegrid_t> grid_th; /* Array of grids for each thread                   */

    std::array<std::vector<int>, DIM> g2t; /* The grid to thread index                         */
    ivec nthread_comm;                     /* The number of threads to communicate with        */
};

/*! \brief Wrapper for gmx_parallel_3dfft_destroy to use as destructor with return type void */
void parallel_3dfft_destroy(gmx_parallel_3dfft* pfft_setup);

//! The data for PME spread/gather plus FFT grids for one set of coefficients
struct PmeAndFftGrids
{
    // Grids on which we do spreading/interpolation, includes overlap
    pmegrids_t pmeGrids;

    // Pointer to the FFT grid, allocated along with pfft_setup
    real* fftgrid = nullptr;

    // Grid for complex FFT data
    t_complex* cfftgrid = nullptr;

    // The setup for the parallel 3D FFT
    gmx::unique_cptr<gmx_parallel_3dfft, parallel_3dfft_destroy> pfft_setup;
};

/*! \brief Data structure for spline-interpolation working buffers */
struct pme_spline_work;

/*! \brief Data structure for working buffers */
struct pme_solve_work_t;

/*! \brief Main PME data structure */
struct gmx_pme_t
{                   //NOLINT(clang-analyzer-optin.performance.Padding)
    int ndecompdim; /* The number of decomposition dimensions */
    int nodeid;     /* Our nodeid in mpi->mpi_comm */
    int nodeid_major;
    int nodeid_minor;
    int nnodes; /* The number of nodes doing PME */
    int nnodes_major;
    int nnodes_minor;

    MPI_Comm mpi_comm;
    MPI_Comm mpi_comm_d[2]; /* Indexed on dimension, 0=x, 1=y */
#if GMX_MPI
    MPI_Datatype rvec_mpi; /* the pme vector's MPI type */
#endif

    bool bUseThreads; /* Does any of the PME ranks have nthread>1 ?  */
    int  nthread;     /* The number of threads doing PME on our rank */

    bool bPPnode;   /* Node also does particle-particle forces */
    bool doCoulomb; /* Apply PME to electrostatics */
    bool doLJ;      /* Apply PME to Lennard-Jones r^-6 interactions */
    bool bFEP;      /* Compute Free energy contribution */
    bool bFEP_q;
    bool bFEP_lj;
    int  nkx, nky, nkz; /* Grid dimensions */
    bool bP3M;          /* Do P3M: optimize the influence function */
    int  pme_order;
    real ewaldcoeff_q;  /* Ewald splitting coefficient for Coulomb */
    real ewaldcoeff_lj; /* Ewald splitting coefficient for r^-6 */
    real epsilon_r;
    int  pmeGpuGridHalo = 0; /* Size of the grid halo region with PME GPU decomposition */
    real haloExtentForAtomDisplacement = .0; /* extent of halo region in nm to account for atom */


    enum PmeRunMode runMode; /* Which codepath is the PME runner taking - CPU, GPU, mixed;
                              * TODO: this is the information that should be owned by the task
                              * scheduler, and ideally not be duplicated here.
                              */

    PmeGpu* gpu; /* A pointer to the GPU data.
                  * TODO: this should be unique or a shared pointer.
                  * Currently in practice there is a single gmx_pme_t instance while a code
                  * is partially set up for many of them. The PME tuning calls gmx_pme_reinit()
                  * which fully reinitializes the one and only PME structure anew while maybe
                  * keeping the old grid buffers if they were already large enough.
                  * This small choice should be made clear in the later refactoring -
                  * do we store many PME objects for different grid sizes,
                  * or a single PME object that handles different grid sizes gracefully.
                  */


    std::unique_ptr<EwaldBoxZScaler> boxScaler; /**< The scaling data Ewald uses with walls (set at pme_init constant for the entire run) */


    LongRangeVdW ljpme_combination_rule; /* Type of combination rule in LJ-PME */

    /* The PME coefficient spreading grid sizes/strides, includes pme_order-1 */
    int pmegrid_nx, pmegrid_ny, pmegrid_nz;
    /* pmegrid_nz might be larger than strictly necessary to ensure
     * memory alignment, pmegrid_nz_base gives the real base size.
     */
    int pmegrid_nz_base;
    /* The local PME grid starting indices */
    int pmegrid_start_ix, pmegrid_start_iy, pmegrid_start_iz;

    /* Work data for spreading and gathering */
    std::unique_ptr<pme_spline_work> spline_work;

    // Storage for all grids in gridCoulomb and GridsLJ, can be shared with other gmx_pme_t objects
    std::shared_ptr<PmeGridsStorage> pmeGridsStorage;

    // PME and FFT grids for Coulomb, one without FEP, two with FEP
    std::vector<PmeAndFftGrids> gridsCoulomb;
    // PME and FFT grids for LJ
    std::vector<PmeAndFftGrids> gridsLJ;

    // Utility struct with references to grids to be used for combined Coulomb+LJ loop
    struct GridsRef
    {
        PmeAndFftGrids& grids;      // The PME and FFT grids
        bool            isCoulomb;  // true: Coulomb grids, false: LJ grids
        int             gridsIndex; // gridsIndex=0: A-coefficients, gridsIndex=1: B-coefficients
    };

    // List of references to Coulomb and/or LJ grids for convenient looping
    std::vector<GridsRef> gridsRefs;

    // List of pointers pointing to the same data as cfftgrids in grids
    std::vector<t_complex*> cfftgrids;

    std::vector<int> nnx;
    std::vector<int> nny;
    std::vector<int> nnz;

    std::vector<real> fshx;
    std::vector<real> fshy;
    std::vector<real> fshz;

    std::vector<PmeAtomComm> atc; /* Indexed on decomposition index */
    matrix                   recipbox;
    real                     boxVolume;
    // The B-spline moduli coefficients
    std::array<std::vector<real>, DIM> bsp_mod;
    /* Buffers to store data for local atoms for L-B combination rule
     * calculations in LJ-PME. lb_buf1 stores either the coefficients
     * for spreading/gathering (in serial), or the C6 coefficient for
     * local atoms (in parallel).  lb_buf2 is only used in parallel,
     * and stores the sigma values for local atoms. */
    FastVector<real> lb_buf1;
    FastVector<real> lb_buf2;

    std::array<pme_overlap_t, 2> overlap; /* Indexed on dimension, 0=x, 1=y */

    /* Atom step for energy only calculation in gmx_pme_calc_energy() */
    std::unique_ptr<PmeAtomComm> atc_energy;

    /* Communication buffers */
    std::vector<gmx::RVec> bufv; /* Communication buffer */
    std::vector<real>      bufr; /* Communication buffer */

    /* thread local work data for solve_pme */
    std::unique_ptr<PmeSolve> pmeSolve;
};

//! @endcond

#endif
