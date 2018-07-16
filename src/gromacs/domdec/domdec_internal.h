/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * \brief Declares implementation functions and types for the domain
 * decomposition module.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_INTERNAL_H
#define GMX_DOMDEC_DOMDEC_INTERNAL_H

#include "config.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/block.h"

struct t_commrec;

/*! \cond INTERNAL */

struct BalanceRegion;

struct gmx_domdec_ind_t
{
    /* The numbers of charge groups to send and receive for each cell
     * that requires communication, the last entry contains the total
     * number of atoms that needs to be communicated.
     */
    int              nsend[DD_MAXIZONE+2];
    int              nrecv[DD_MAXIZONE+2];
    /* The charge groups to send */
    std::vector<int> index;
    /* The atom range for non-in-place communication */
    int              cell2at0[DD_MAXIZONE];
    int              cell2at1[DD_MAXIZONE];
};

struct gmx_domdec_comm_dim_t
{
    /* Returns the number of grid pulses (the number of domains in the halo along this dimension) */
    int numPulses() const
    {
        return ind.size();
    }

    int                           np_dlb;         /* For dlb, for use with edlbAUTO          */
    std::vector<gmx_domdec_ind_t> ind;            /* The indices to communicate, size np     */
    bool                          receiveInPlace; /* Can we receive data in place?            */
};

/*! \brief Load balancing data along a dim used on the master rank of that dim */
struct RowMaster
{
    struct Bounds
    {
        real cellFracLowerMax; /**< State var.: max lower bound., incl. neighbors */
        real cellFracUpperMin; /**< State var.: min upper bound., incl. neighbors */
        real boundMin;         /**< Temp. var.: lower limit for cell boundary */
        real boundMax;         /**< Temp. var.: upper limit for cell boundary */
    };

    std::vector<bool>   isCellMin;    /**< Temp. var.: is this cell size at the limit */
    std::vector<real>   cellFrac;     /**< State var.: cell boundaries, box relative */
    std::vector<real>   oldCellFrac;  /**< Temp. var.: old cell size */
    std::vector<Bounds> bounds;       /**< Cell bounds */
    bool                dlbIsLimited; /**< State var.: is DLB limited in this row */
    std::vector<real>   buf_ncd;      /**< Temp. var.  */
};

/*! \brief Struct for managing cell sizes with DLB along a dimension */
struct DDCellsizesWithDlb
{
    /* Cell sizes for dynamic load balancing */
    std::unique_ptr<RowMaster> rowMaster;    /**< Cell row root struct, only available on the first rank in a row */
    std::vector<real>          fracRow;      /**< The cell sizes, in fractions, along a row, not available on the first rank in a row */
    real                       fracLower;    /**< The lower corner, in fractions, in triclinic space */
    real                       fracUpper;    /**< The upper corner, in fractions, in triclinic space */
    real                       fracLowerMax; /**< The maximum lower corner among all our neighbors */
    real                       fracUpperMin; /**< The minimum upper corner among all our neighbors */
};

/*! \brief Struct for compute load commuication
 *
 * Here floats are accurate enough, since these variables
 * only influence the load balancing, not the actual MD results.
 */
typedef struct
{
    int    nload;     /**< The number of load recordings */
    float *load;      /**< Scan of the sum of load over dimensions */
    float  sum;       /**< The sum of the load over the ranks up to our current dimension */
    float  max;       /**< The maximum over the ranks contributing to \p sum */
    float  sum_m;     /**< Like \p sum, but takes the maximum when the load balancing is limited */
    float  cvol_min;  /**< Minimum cell volume, relative to the box */
    float  mdf;       /**< The PP time during which PME can overlap */
    float  pme;       /**< The PME-only rank load */
    int    flags;     /**< Bit flags that tell if DLB was limited, per dimension */
} domdec_load_t;

/*! \brief Data needed to sort an atom to the desired location in the local state */
typedef struct
{
    int  nsc;     /**< Neighborsearch grid cell index */
    int  ind_gl;  /**< Global atom/charge group index */
    int  ind;     /**< Local atom/charge group index */
} gmx_cgsort_t;

/*! \brief Temporary buffers for sorting atoms */
typedef struct
{
    std::vector<gmx_cgsort_t> sorted;     /**< Sorted array of indices */
    std::vector<gmx_cgsort_t> stationary; /**< Array of stationary atom/charge group indices */
    std::vector<gmx_cgsort_t> moved;      /**< Array of moved atom/charge group indices */
    std::vector<int>          intBuffer;  /**< Integer buffer for sorting */
} gmx_domdec_sort_t;

/*! \brief Manages atom ranges and order for the local state atom vectors */
class DDAtomRanges
{
    public:
        /*! \brief The local state atom order
         *
         * This enum determines the order of the atoms in the local state.
         * ddnatHOME and ddnatZONE should be first and second,
         * the others can be ordered as wanted.
         */
        enum class Type : int
        {
            Home,        /**< The home atoms */
            Zones,       /**< All zones in the eighth shell */
            Vsites,      /**< Atoms communicated for virtual sites */
            Constraints, /**< Atoms communicated for constraints */
            Number       /**< Not a count, only present for convenience */
        };

        /*! \brief Returns the start atom index for range \p rangeType */
        int start(Type rangeType) const
        {
            if (rangeType == Type::Home)
            {
                return 0;
            }
            else
            {
                return end_[static_cast<int>(rangeType) - 1];
            }
        }

        /*! \brief Returns the end atom index for range \p rangeType */
        int end(Type rangeType) const
        {
            return end_[static_cast<int>(rangeType)];
        }

        /*! \brief Returns the number of home atoms */
        int numHomeAtoms() const
        {
            return end_[static_cast<int>(Type::Home)];
        }

        /*! \brief Returns the total number of atoms */
        int numAtomsTotal() const
        {
            return end_[static_cast<int>(Type::Number) - 1];
        }

        /*! \brief Sets the end index of range \p rangeType to \p end
         *
         * This should be called either with Type::Home or with a type
         * that is larger than that passed in the previous call to setEnd.
         * A release assertion for these conditions are present.
         */
        void setEnd(Type rangeType,
                    int  end)
        {
            GMX_RELEASE_ASSERT(rangeType == Type::Home || rangeType > lastTypeSet_, "Can only set either home or a larger type than the last one");

            for (int i = static_cast<int>(rangeType); i < static_cast<int>(Type::Number); i++)
            {
                end_[i] = end;
            }

            lastTypeSet_ = rangeType;
        }

    private:
        /*! \brief The list of end atom indices */
        std::array<int, static_cast<int>(Type::Number)> end_;
        Type                                            lastTypeSet_ = Type::Number;
};

/*! \brief Enum of dynamic load balancing states
 *
 * Allowed DLB states and transitions
 * - intialization at startup:
 *                             -> edlbsOffUser ("-dlb no")
 *                             -> edlbsOnUser  ("-dlb yes")
 *                             -> edlbsOffCanTurnOn ("-dlb auto")
 *
 * - in automatic mode (i.e. initial state edlbsOffCanTurnOn):
 *   edlbsOffCanTurnOn         -> edlbsOnCanTurnOff
 *   edlbsOffCanTurnOn         -> edlbsOffForever
 *   edlbsOffCanTurnOn         -> edlbsOffTemporarilyLocked
 *   edlbsOffTemporarilyLocked -> edlbsOffCanTurnOn
 *   edlbsOnCanTurnOff         -> edlbsOffCanTurnOn
 */
enum {
    edlbsOffUser,              /**< DLB is permanently off per user request */
    edlbsOffForever,           /**< DLB is off due to a runtime condition (not supported or causes performance loss) and will never be turned on */
    edlbsOffCanTurnOn,         /**< DLB is off and will turn on on imbalance */
    edlbsOffTemporarilyLocked, /**< DLB is off and temporarily can't turn on */
    edlbsOnCanTurnOff,         /**< DLB is on and can turn off when slow */
    edlbsOnUser,               /**< DLB is permanently on per user request */
    edlbsNR                    /**< The number of DLB states */
};

/*! \brief The PME domain decomposition for one dimension */
typedef struct
{
    int      dim;       /**< The dimension */
    gmx_bool dim_match; /**< Tells if DD and PME dims match */
    int      nslab;     /**< The number of PME ranks/domains in this dimension */
    real    *slb_dim_f; /**< Cell sizes for determining the PME comm. with SLB */
    int     *pp_min;    /**< The minimum pp node location, size nslab */
    int     *pp_max;    /**< The maximum pp node location, size nslab */
    int      maxshift;  /**< The maximum shift for coordinate redistribution in PME */
} gmx_ddpme_t;

struct gmx_ddzone_t
{
    real min0;    /* The minimum bottom of this zone                        */
    real max1;    /* The maximum top of this zone                           */
    real min1;    /* The minimum top of this zone                           */
    real mch0;    /* The maximum bottom communicaton height for this zone   */
    real mch1;    /* The maximum top communicaton height for this zone      */
    real p1_0;    /* The bottom value of the first cell in this zone        */
    real p1_1;    /* The top value of the first cell in this zone           */
};

/*! \brief The number of reals in gmx_ddzone_t */
constexpr int c_ddzoneNumReals = 7;

/*! \brief Forward declaration */
template<typename T> class DDBufferAccess;

/*! \brief Temporary storage container that minimizes (re)allocation and clearing
 *
 * This is only the storage, actual access happens through DDBufferAccess.
 * All methods check if the buffer is (not) in use.
 */
template<typename T>
class DDBuffer
{
    private:
        /*! \brief Returns a buffer of size \p numElements, the elements are undefined */
        gmx::ArrayRef<T> resize(size_t numElements)
        {
            GMX_ASSERT(isInUse_, "Should only operate on acquired buffers");

            if (numElements > buffer_.size())
            {
                buffer_.resize(numElements);
            }

            return gmx::arrayRefFromArray(buffer_.data(), numElements);
        }

        /*! \brief Acquire the buffer for use with size set to \p numElements, the elements are undefined */
        gmx::ArrayRef<T> acquire(size_t numElements)
        {
            GMX_RELEASE_ASSERT(!isInUse_, "Should only request free buffers");
            isInUse_ = true;

            return resize(numElements);
        }

        /*! \brief Releases the buffer, buffer_ should not be used after this */
        void release()
        {
            GMX_RELEASE_ASSERT(isInUse_, "Should only release buffers in use");
            isInUse_ = false;
        }

        std::vector<T> buffer_;          /**< The actual memory buffer */
        bool           isInUse_ = false; /**< Flag that tells whether the buffer is in use */

        friend class DDBufferAccess<T>;
};

/*! \brief Class that manages access to a temporary memory buffer */
template<typename T>
class DDBufferAccess
{
    public:
        /*! \brief Constructor, returns a buffer of size \p numElements, element values are undefined
         *
         * \note The actual memory buffer \p ddBuffer can not be used to
         *       create other DDBufferAccess objects until the one created
         *       here is destroyed.
         */
        DDBufferAccess(DDBuffer<T> &ddBuffer,
                       size_t       numElements) :
            ddBuffer_(ddBuffer)
        {
            buffer = ddBuffer_.acquire(numElements);
        }

        /*! \brief Destructor */
        ~DDBufferAccess()
        {
            ddBuffer_.release();
        }

        /*! \brief Resizes the buffer to \p numElements, new elements are undefined
         *
         * \note The buffer arrayref is updated after this call.
         */
        void resize(size_t numElements)
        {
            buffer = ddBuffer_.resize(numElements);
        }

    private:
        DDBuffer<T>      &ddBuffer_; /**< Reference to the storage class */
    public:
        gmx::ArrayRef<T>  buffer;    /**< The access to the memory buffer */
};

/*! brief Temporary buffer for setting up communiation over one pulse and all zones in the halo */
struct dd_comm_setup_work_t
{
    std::vector<int>       localAtomGroupBuffer; /**< The local atom group indices to send */
    std::vector<int>       atomGroupBuffer;      /**< Buffer for collecting the global atom group indices to send */
    std::vector<gmx::RVec> positionBuffer;       /**< Buffer for collecting the atom group positions to send */
    int                    nat;                  /**< The number of atoms contained in the atom groups to send */
    int                    nsend_zone;           /**< The number of atom groups to send for the last zone */
};

/*! \brief Struct for domain decomposition communication
 *
 * This struct contains most information about domain decomposition
 * communication setup, some communication buffers, some statistics
 * and also the setup for the communication between particle-particle
 * and PME only ranks.
 *
 * All arrays are indexed with 0 to dd->ndim (not Cartesian indexing),
 * unless stated otherwise.
 */
struct gmx_domdec_comm_t
{
    /* PME and Cartesian communicator stuff */
    int         npmedecompdim;     /**< The number of decomposition dimensions for PME, 0: no PME */
    int         npmenodes;         /**< The number of ranks doing PME (PP/PME or only PME) */
    int         npmenodes_x;       /**< The number of PME ranks/domains along x */
    int         npmenodes_y;       /**< The number of PME ranks/domains along y */
    gmx_bool    bCartesianPP_PME;  /**< Use Cartesian communication between PP and PME ranks */
    ivec        ntot;              /**< Cartesian grid for combinted PP+PME ranks */
    int         cartpmedim;        /**< The number of dimensions for the PME setup that are Cartesian */
    int        *pmenodes;          /**< The PME ranks, size npmenodes */
    int        *ddindex2simnodeid; /**< The Cartesian index to sim rank conversion, used with bCartesianPP_PME */
    gmx_ddpme_t ddpme[2];          /**< The 1D or 2D PME domain decomposition setup */

    /* The DD particle-particle nodes only */
    gmx_bool bCartesianPP;        /**< Use a Cartesian communicator for PP */
    int     *ddindex2ddnodeid;    /**< The Cartesian index to DD rank conversion, used with bCartesianPP */

    /* The DLB state, used for reloading old states, during e.g. EM */
    t_block cgs_gl;               /**< The global charge groups, this defined the DD state (except for the DLB state) */

    /* Charge group / atom sorting */
    std::unique_ptr<gmx_domdec_sort_t> sort; /**< Data structure for cg/atom sorting */

    /* Are there charge groups? */
    gmx_bool bCGs;                /**< True when there are charge groups */

    gmx_bool bInterCGBondeds;     /**< Are there inter-cg bonded interactions? */
    gmx_bool bInterCGMultiBody;   /**< Are there inter-cg multi-body interactions? */

    /* Data for the optional bonded interaction atom communication range */
    gmx_bool  bBondComm;          /**< Only communicate atoms beyond the non-bonded cut-off when they are involved in bonded interactions with non-local atoms */
    t_blocka *cglink;             /**< Links between cg's through bonded interactions */
    char     *bLocalCG;           /**< Local cg availability, TODO: remove when group scheme is removed */

    /* The DLB state, possible values are defined above */
    int      dlbState;
    /* With dlbState=edlbsOffCanTurnOn, should we check if to DLB on at the next DD? */
    gmx_bool bCheckWhetherToTurnDlbOn;
    /* The first DD count since we are running without DLB */
    int      ddPartioningCountFirstDlbOff = 0;

    /* Cell sizes for static load balancing, first index cartesian */
    real **slb_frac;

    /* The width of the communicated boundaries */
    real     cutoff_mbody;        /**< Cut-off for multi-body interactions, also 2-body bonded when \p cutoff_mody > \p cutoff */
    real     cutoff;              /**< Cut-off for non-bonded/2-body interactions */
    rvec     cellsize_min;        /**< The minimum guaranteed cell-size, Cartesian indexing */
    rvec     cellsize_min_dlb;    /**< The minimum guaranteed cell-size with dlb=auto */
    real     cellsize_limit;      /**< The lower limit for the DD cell size with DLB */
    gmx_bool bVacDLBNoLimit;      /**< Effectively no NB cut-off limit with DLB for systems without PBC? */

    /** With PME load balancing we set limits on DLB */
    gmx_bool bPMELoadBalDLBLimits;
    /** DLB needs to take into account that we want to allow this maximum
     *  cut-off (for PME load balancing), this could limit cell boundaries.
     */
    real PMELoadBal_max_cutoff;

    ivec tric_dir;                /**< tric_dir from \p gmx_domdec_box_t is only stored here because dd_get_ns_ranges needs it */
    rvec box0;                    /**< box lower corner, required with dim's without pbc and -gcom */
    rvec box_size;                /**< box size, required with dim's without pbc and -gcom */

    rvec cell_x0;                 /**< The DD cell lower corner, in triclinic space */
    rvec cell_x1;                 /**< The DD cell upper corner, in triclinic space */

    rvec old_cell_x0;             /**< The old \p cell_x0, to check cg displacements */
    rvec old_cell_x1;             /**< The old \p cell_x1, to check cg displacements */

    /** The communication setup and charge group boundaries for the zones */
    gmx_domdec_zones_t zones;

    /* The zone limits for DD dimensions 1 and 2 (not 0), determined from
     * cell boundaries of neighboring cells for staggered grids when using
     * dynamic load balancing.
     */
    gmx_ddzone_t zone_d1[2];          /**< Zone limits for dim 1 with staggered grids */
    gmx_ddzone_t zone_d2[2][2];       /**< Zone limits for dim 2 with staggered grids */

    /** The coordinate/force communication setup and indices */
    gmx_domdec_comm_dim_t cd[DIM];
    /** The maximum number of cells to communicate with in one dimension */
    int                   maxpulse;

    /** Which cg distribution is stored on the master node,
     *  stored as DD partitioning call count.
     */
    int64_t master_cg_ddp_count;

    /** The number of cg's received from the direct neighbors */
    int  zone_ncg1[DD_MAXZONE];

    /** The atom ranges in the local state */
    DDAtomRanges atomRanges;

    /** Array for signalling if atoms have moved to another domain */
    std::vector<int> movedBuffer;

    /** Communication int buffer for general use */
    DDBuffer<int> intBuffer;

    /** Communication rvec buffer for general use */
    DDBuffer<gmx::RVec> rvecBuffer;

    /* Temporary storage for thread parallel communication setup */
    std::vector<dd_comm_setup_work_t> dth; /**< Thread-local work data */

    /* Communication buffer only used with multiple grid pulses */
    DDBuffer<gmx::RVec> rvecBuffer2; /**< Another rvec comm. buffer */

    /* Communication buffers for local redistribution */
    std::array<std::vector<int>, DIM*2>       cggl_flag;  /**< Charge group flag comm. buffers */
    std::array<std::vector<gmx::RVec>, DIM*2> cgcm_state; /**< Charge group center comm. buffers */

    /* Cell sizes for dynamic load balancing */
    std::vector<DDCellsizesWithDlb> cellsizesWithDlb;

    /* Stuff for load communication */
    gmx_bool        bRecordLoad;         /**< Should we record the load */
    domdec_load_t  *load;                /**< The recorded load data */
    int             nrank_gpu_shared;    /**< The number of MPI ranks sharing the GPU our rank is using */
#if GMX_MPI
    MPI_Comm       *mpi_comm_load;       /**< The MPI load communicator */
    MPI_Comm        mpi_comm_gpu_shared; /**< The MPI load communicator for ranks sharing a GPU */
#endif

    /* Information for managing the dynamic load balancing */
    int            dlb_scale_lim;      /**< Maximum DLB scaling per load balancing step in percent */

    BalanceRegion *balanceRegion;      /**< Struct for timing the force load balancing region */

    /* Cycle counters over nstlist steps */
    float  cycl[ddCyclNr];             /**< Total cycles counted */
    int    cycl_n[ddCyclNr];           /**< The number of cycle recordings */
    float  cycl_max[ddCyclNr];         /**< The maximum cycle count */
    /** Flop counter (0=no,1=yes,2=with (eFlop-1)*5% noise */
    int    eFlop;
    double flop;                       /**< Total flops counted */
    int    flop_n;                     /**< The number of flop recordings */
    /** How many times did we have load measurements */
    int    n_load_have;
    /** How many times have we collected the load measurements */
    int    n_load_collect;

    /* Cycle count history for DLB checks */
    float       cyclesPerStepBeforeDLB;     /**< The averaged cycles per step over the last nstlist step before turning on DLB */
    float       cyclesPerStepDlbExpAverage; /**< The running average of the cycles per step during DLB */
    bool        haveTurnedOffDlb;           /**< Have we turned off DLB (after turning DLB on)? */
    int64_t     dlbSlowerPartitioningCount; /**< The DD step at which we last measured that DLB off was faster than DLB on, 0 if there was no such step */

    /* Statistics for atoms */
    double sum_nat[static_cast<int>(DDAtomRanges::Type::Number)]; /**< The atoms per range, summed over the steps */

    /* Statistics for calls and times */
    int    ndecomp;                    /**< The number of partioning calls */
    int    nload;                      /**< The number of load recordings */
    double load_step;                  /**< Total MD step time */
    double load_sum;                   /**< Total PP force time */
    double load_max;                   /**< Max \p load_sum over the ranks */
    ivec   load_lim;                   /**< Was load balancing limited, per DD dim */
    double load_mdf;                   /**< Total time on PP done during PME overlap time */
    double load_pme;                   /**< Total time on our PME-only rank */

    /** The last partition step */
    int64_t partition_step;

    /* Debugging */
    int  nstDDDump;                    /**< Step interval for dumping the local+non-local atoms to pdb */
    int  nstDDDumpGrid;                /**< Step interval for duming the DD grid to pdb */
    int  DD_debug;                     /**< DD debug print level: 0, 1, 2 */
};

/*! \brief DD zone permutation
 *
 * Zone permutation from the Cartesian x-major/z-minor order to an order
 * that leads to consecutive charge groups for neighbor searching.
 * TODO: remove when the group scheme is removed
 */
static const int zone_perm[3][4] = { {0, 0, 0, 0}, {1, 0, 0, 0}, {3, 0, 1, 2} };

/*! \brief DD zone reordering to Cartesian order
 *
 * Index to reorder the zone such that the end up in Cartesian order
 * with dimension index 0 major and dimension index 2 minor.
 */
static const int zone_reorder_cartesian[DD_MAXZONE] = { 0, 1, 3, 2, 5, 4, 6, 7 };

/* dd_zo and dd_zp3 is set up such that i zones with non-zero
 * components see only j zones with that component 0.
 */

/*! \brief Returns the DD cut-off distance for multi-body interactions */
real dd_cutoff_multibody(const gmx_domdec_t *dd);

/*! \brief Returns the DD cut-off distance for two-body interactions */
real dd_cutoff_twobody(const gmx_domdec_t *dd);

/*! \brief Returns the domain index given the number of domains and the domain coordinates
 *
 * This order is required to minimize the coordinate communication in PME
 * which uses decomposition in the x direction.
 */
static inline int dd_index(const ivec numDomains,
                           const ivec domainCoordinates)
{
    return ((domainCoordinates[XX]*numDomains[YY] + domainCoordinates[YY])*numDomains[ZZ]) + domainCoordinates[ZZ];
};

/*! Returns the size of the buffer to hold fractional cell boundaries for DD dimension index dimIndex */
static inline int ddCellFractionBufferSize(const gmx_domdec_t *dd,
                                           int                 dimIndex)
{
    return dd->nc[dd->dim[dimIndex ]] + 1 + dimIndex*2 + 1 + dimIndex;
}

/*! \brief Maximum number of ranks for using send/recv for state scattering and gathering
 *
 * Use separate MPI send and receive commands
 * when #ranks <= c_maxNumRanksUseSendRecvForScatterAndGather
 * This saves memory (and some copying for small #ranks).
 * For high parallelization scatter and gather calls are used.
 */
static constexpr int c_maxNumRanksUseSendRecvForScatterAndGather = 4;

/*! \brief Make DD cells larger by this factor than the limit to avoid rounding issues */
static constexpr double DD_CELL_MARGIN       = 1.0001;

/*! \brief Factor for checking DD cell size limitation during DLB, should be in between 1 and DD_CELL_MARGIN */
static constexpr double DD_CELL_MARGIN2      = 1.00005;

/*! \brief With pressure scaling, keep cell sizes 2% above the limit to allow for some scaling */
static constexpr double DD_PRES_SCALE_MARGIN = 1.02;

/*! \endcond */

#endif
