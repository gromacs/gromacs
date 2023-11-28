/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief Declares implementation functions and types for the domain
 * decomposition module.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_INTERNAL_H
#define GMX_DOMDEC_DOMDEC_INTERNAL_H

#include "config.h"

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/listoflists.h"

struct t_commrec;

/*! \cond INTERNAL */

#define DD_NLOAD_MAX 9

namespace gmx
{
enum class DdRankOrder : int;
}
// namespace


//! Indices to communicate in a dimension
struct gmx_domdec_ind_t
{
    //! @{
    /*! \brief The numbers of charge groups to send and receive for each
     * cell that requires communication, the last entry contains the total
     * number of atoms that needs to be communicated.
     */
    int nsend[DD_MAXIZONE + 2] = {};
    int nrecv[DD_MAXIZONE + 2] = {};
    //! @}
    //! The charge groups to send
    std::vector<int> index;
    //! @{
    /* The atom range for non-in-place communication */
    int cell2at0[DD_MAXIZONE] = {};
    int cell2at1[DD_MAXIZONE] = {};
    //! @}
};

//! Things relating to index communication
struct gmx_domdec_comm_dim_t
{
    /* Returns the number of grid pulses (the number of domains in the halo along this dimension) */
    int numPulses() const { return ind.size(); }

    /**< For dlb, for use with edlbAUTO          */
    int np_dlb = 0;
    /**< The indices to communicate, size np     */
    std::vector<gmx_domdec_ind_t> ind;
    /**< Can we receive data in place?            */
    bool receiveInPlace = false;
};

/*! \brief Load balancing data along a dim used on the main rank of that dim */
struct RowCoordinator
{
    struct Bounds
    {
        /**< State var.: max lower bound., incl. neighbors */
        real cellFracLowerMax = 0;
        /**< State var.: min upper bound., incl. neighbors */
        real cellFracUpperMin = 0;
        /**< Temp. var.: lower limit for cell boundary */
        real boundMin = 0;
        /**< Temp. var.: upper limit for cell boundary */
        real boundMax = 0;
    };

    /**< Temp. var.: is this cell size at the limit */
    std::vector<bool> isCellMin;
    /**< State var.: cell boundaries, box relative */
    std::vector<real> cellFrac;
    /**< Temp. var.: old cell size */
    std::vector<real> oldCellFrac;
    /**< Cell bounds */
    std::vector<Bounds> bounds;
    /**< State var.: is DLB limited in this row */
    bool dlbIsLimited = false;
    /**< Temp. var.  */
    std::vector<real> buf_ncd;
};

/*! \brief Struct for managing cell sizes with DLB along a dimension */
struct DDCellsizesWithDlb
{
    /**< Cell row root struct, only available on the first rank in a row */
    std::unique_ptr<RowCoordinator> rowCoordinator;
    /**< The cell sizes, in fractions, along a row, not available on the first rank in a row */
    std::vector<real> fracRow;
    /**< The lower corner, in fractions, in triclinic space */
    real fracLower = 0;
    /**< The upper corner, in fractions, in triclinic space */
    real fracUpper = 0;
    /**< The maximum lower corner among all our neighbors */
    real fracLowerMax = 0;
    /**< The minimum upper corner among all our neighbors */
    real fracUpperMin = 0;
};

/*! \brief Struct for compute load commuication
 *
 * Here floats are accurate enough, since these variables
 * only influence the load balancing, not the actual MD results.
 */
typedef struct domdec_load
{
    /**< The number of load recordings */
    int nload = 0;
    /**< Scan of the sum of load over dimensions */
    std::vector<float> load;
    /**< The sum of the load over the ranks up to our current dimension */
    float sum = 0;
    /**< The maximum over the ranks contributing to \p sum */
    float max = 0;
    /**< Like \p sum, but takes the maximum when the load balancing is limited */
    float sum_m = 0;
    /**< Minimum cell volume, relative to the box */
    float cvol_min = 0;
    /**< The PP time during which PME can overlap */
    float mdf = 0;
    /**< The PME-only rank load */
    float pme = 0;
    /**< Bit flags that tell if DLB was limited, per dimension */
    int flags = 0;
} domdec_load_t;

/*! \brief Data needed to sort an atom to the desired location in the local state */
typedef struct gmx_cgsort
{
    /**< Local atom/charge group index */
    int ind = 0;
} gmx_cgsort_t;

/*! \brief Temporary buffers for sorting atoms */
typedef struct gmx_domdec_sort
{
    /**< Sorted array of indices */
    std::vector<gmx_cgsort_t> sorted;
    /**< Array of stationary atom/charge group indices */
    std::vector<gmx_cgsort_t> stationary;
    /**< Array of moved atom/charge group indices */
    std::vector<gmx_cgsort_t> moved;
    /**< Integer buffer for sorting */
    std::vector<int> intBuffer;
    /**< Int64 buffer for sorting */
    std::vector<int64_t> int64Buffer;
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
    int end(Type rangeType) const { return end_[static_cast<int>(rangeType)]; }

    /*! \brief Returns the number of home atoms */
    int numHomeAtoms() const { return end_[static_cast<int>(Type::Home)]; }

    /*! \brief Returns the total number of atoms */
    int numAtomsTotal() const { return end_[static_cast<int>(Type::Number) - 1]; }

    /*! \brief Sets the end index of range \p rangeType to \p end
     *
     * This should be called either with Type::Home or with a type
     * that is larger than that passed in the previous call to setEnd.
     * A release assertion for these conditions are present.
     */
    void setEnd(Type rangeType, int end)
    {
        GMX_RELEASE_ASSERT(rangeType == Type::Home || rangeType > lastTypeSet_,
                           "Can only set either home or a larger type than the last one");

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
 * - initialization at startup:
 *                             -> offUser ("-dlb no")
 *                             -> onUser  ("-dlb yes")
 *                             -> offCanTurnOn ("-dlb auto")
 *
 * - in automatic mode (i.e. initial state offCanTurnOn):
 *   offCanTurnOn         -> onCanTurnOff
 *   offCanTurnOn         -> offForever
 *   offCanTurnOn         -> offTemporarilyLocked
 *   offTemporarilyLocked -> offCanTurnOn
 *   onCanTurnOff         -> offCanTurnOn
 */
enum class DlbState
{
    offUser,    /**< DLB is permanently off per user request */
    offForever, /**< DLB is off due to a runtime condition (not supported or causes performance loss) and will never be turned on */
    offCanTurnOn,         /**< DLB is off and will turn on on imbalance */
    offTemporarilyLocked, /**< DLB is off and temporarily can't turn on */
    onCanTurnOff,         /**< DLB is on and can turn off when slow */
    onUser,               /**< DLB is permanently on per user request */
    Count                 /**< The number of DLB states */
};

/*! \brief The PME domain decomposition for one dimension */
typedef struct gmx_ddpme
{
    /**< The dimension */
    int dim = 0;
    /**< Tells if DD and PME dims match */
    gmx_bool dim_match = false;
    /**< The number of PME ranks/domains in this dimension */
    int nslab = 0;
    /**< Cell sizes for determining the PME comm. with SLB */
    std::vector<real> slb_dim_f;
    /**< The minimum pp node location, size nslab */
    std::vector<int> pp_min;
    /**< The maximum pp node location, size nslab */
    std::vector<int> pp_max;
    /**< The maximum shift for coordinate redistribution in PME */
    int maxshift = 0;
} gmx_ddpme_t;

struct gmx_ddzone_t
{
    /**< The minimum bottom of this zone                        */
    real min0 = 0;
    /**< The maximum top of this zone                           */
    real max1 = 0;
    /**< The minimum top of this zone                           */
    real min1 = 0;
    /**< The maximum bottom communicaton height for this zone   */
    real mch0 = 0;
    /**< The maximum top communicaton height for this zone      */
    real mch1 = 0;
    /**< The bottom value of the first cell in this zone        */
    real p1_0 = 0;
    /**< The top value of the first cell in this zone           */
    real p1_1 = 0;
    /**< Bool disguised as a real, 1 when the above data has been set. 0 otherwise */
    real dataSet = 0;
};

/*! \brief The number of reals in gmx_ddzone_t */
constexpr int c_ddzoneNumReals = 8;

template<typename T>
class DDBufferAccess;

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
    DDBufferAccess(DDBuffer<T>& ddBuffer, size_t numElements) : ddBuffer_(ddBuffer)
    {
        buffer = ddBuffer_.acquire(numElements);
    }

    ~DDBufferAccess() { ddBuffer_.release(); }

    /*! \brief Resizes the buffer to \p numElements, new elements are undefined
     *
     * \note The buffer arrayref is updated after this call.
     */
    void resize(size_t numElements) { buffer = ddBuffer_.resize(numElements); }

private:
    DDBuffer<T>& ddBuffer_; /**< Reference to the storage class */
public:
    gmx::ArrayRef<T> buffer; /**< The access to the memory buffer */
};

/*! \brief Temporary buffer for setting up communiation over one pulse and all zones in the halo */
struct dd_comm_setup_work_t
{
    /**< The local atom group indices to send */
    std::vector<int> localAtomGroupBuffer;
    /**< Buffer for collecting the global atom group indices to send */
    std::vector<int> atomGroupBuffer;
    /**< Buffer for collecting the atom group positions to send */
    std::vector<gmx::RVec> positionBuffer;
    /**< The number of atoms contained in the atom groups to send */
    int nat = 0;
    /**< The number of atom groups to send for the last zone */
    int nsend_zone = 0;
};

/*! \brief Information about the simulated system */
struct DDSystemInfo
{
    //! True when update groups are used
    bool useUpdateGroups = false;
    //! Update atom grouping for each molecule type
    gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType;
    //! The maximum radius over all update groups
    real maxUpdateGroupRadius;

    //! Are molecules always whole, i.e. not broken by PBC?
    bool moleculesAreAlwaysWhole = false;
    //! Are there inter-domain bonded interactions?
    bool haveInterDomainBondeds = false;
    //! Are there inter-domain multi-body interactions?
    bool haveInterDomainMultiBodyBondeds = false;

    //! Cut-off for multi-body interactions
    real minCutoffForMultiBody = 0;
    //! Cut-off for non-bonded/2-body interactions
    real cutoff = 0;
    //! The lower limit for the DD cell size
    real cellsizeLimit = 0;

    //! Can atoms connected by constraints be assigned to different domains?
    bool mayHaveSplitConstraints = false;
    //! Can atoms connected by settles be assigned to different domains?
    bool mayHaveSplitSettles = false;
    //! Estimated communication range needed for constraints
    real constraintCommunicationRange = 0;

    //! Whether to only communicate atoms beyond the non-bonded cut-off when they are involved in bonded interactions with non-local atoms
    bool filterBondedCommunication = false;
    //! Whether to increase the multi-body cut-off beyond the minimum required
    bool increaseMultiBodyCutoff = false;

    //! Whether we have continous box deformation
    bool haveBoxDeformation;
    //! The box deformation rate in units of 1/ps
    matrix boxDeformationRate;
};

/*! \brief Settings that affect the behavior of the domain decomposition
 *
 * These settings depend on options chosen by the user, set by enviroment
 * variables, as well as hardware support. The initial DLB state also
 * depends on the integrator.
 *
 * Note: Settings that depend on the simulated system are in DDSystemInfo.
 */
struct DDSettings
{
    //! Use MPI_Sendrecv communication instead of non-blocking calls
    bool useSendRecv2 = false;

    /* Information for managing the dynamic load balancing */
    //! Maximum DLB scaling per load balancing step in percent
    int dlb_scale_lim = 0;
    //! Flop counter (0=no,1=yes,2=with (eFlop-1)*5% noise
    int eFlop = 0;

    //! Whether to order the DD dimensions from z to x
    bool useDDOrderZYX = false;

    //! Whether to use MPI Cartesian reordering of communicators, when supported (almost never)
    bool useCartesianReorder = true;

    //! Whether we should record the load
    bool recordLoad = false;

    /* Debugging */
    //! Step interval for dumping the local+non-local atoms to pdb
    int nstDDDump = 0;
    //! Step interval for duming the DD grid to pdb
    int nstDDDumpGrid = 0;
    //! DD debug print level: 0, 1, 2
    int DD_debug = 0;

    //! The DLB state at the start of the run
    DlbState initialDlbState = DlbState::offCanTurnOn;
};

/*! \brief Information on how the DD ranks are set up */
// The following suppression suppresses an error: "declaration uses
// identifier '__i0', which is a reserved identifier" which does not
// make sense from the code, but is not yet a known clang-tidy bug.
struct DDRankSetup //NOLINT(bugprone-reserved-identifier,google-readability-braces-around-statements,readability-braces-around-statements)
{
    /**< The rank ordering */
    gmx::DdRankOrder rankOrder;

    /**< The number of particle-particle (non PME-only) ranks */
    int numPPRanks = 0;
    /**< The DD PP grid */
    ivec numPPCells = { 0, 0, 0 };

    /* PME and Cartesian communicator stuff */
    bool usePmeOnlyRanks = false;
    /**< The number of decomposition dimensions for PME, 0: no PME */
    int npmedecompdim = 0;
    /**< The number of ranks doing PME (PP/PME or only PME) */
    int numRanksDoingPme = 0;
    /**< The number of PME ranks/domains along x */
    int npmenodes_x = 0;
    /**< The number of PME ranks/domains along y */
    int npmenodes_y = 0;
    /**< The 1D or 2D PME domain decomposition setup */
    gmx_ddpme_t ddpme[2];
};

/*! \brief Information on Cartesian MPI setup of the DD ranks */
struct CartesianRankSetup
{
    /**< Use Cartesian communication between PP and PME ranks */
    bool bCartesianPP_PME = false;
    /**< Cartesian grid for combinted PP+PME ranks */
    ivec ntot = {};
    /**< The number of dimensions for the PME setup that are Cartesian */
    int cartpmedim = 0;
    /**< The Cartesian index to sim rank conversion, used with bCartesianPP_PME */
    std::vector<int> ddindex2simnodeid;

    /* The DD particle-particle nodes only */
    /**< Use a Cartesian communicator for PP */
    bool bCartesianPP = false;
    /**< The Cartesian index to DD rank conversion, used with bCartesianPP */
    std::vector<int> ddindex2ddnodeid;
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
struct gmx_domdec_comm_t // NOLINT (clang-analyzer-optin.performance.Padding)
{
    /**< Constant parameters that control DD behavior */
    DDSettings ddSettings;

    /**< Information on how the DD ranks are set up */
    DDRankSetup ddRankSetup;
    /**< Information on the Cartesian part of the DD rank setup */
    CartesianRankSetup cartesianRankSetup;

    /* Charge group / atom sorting */
    /**< Data structure for cg/atom sorting */
    std::unique_ptr<gmx_domdec_sort_t> sort;

    //! Centers of mass of local update groups
    std::unique_ptr<gmx::UpdateGroupsCog> updateGroupsCog;

    /* Data for the optional filtering of communication of atoms for bonded interactions */
    /**< Links between atoms through bonded interactions */
    std::unique_ptr<gmx::ListOfLists<int>> bondedLinks;

    /* The DLB state, possible values are defined above */
    DlbState dlbState;
    /* With dlbState=DlbState::offCanTurnOn, should we check if to DLB on at the next DD? */
    gmx_bool bCheckWhetherToTurnDlbOn = false;
    /* The first DD count since we are running without DLB */
    int ddPartioningCountFirstDlbOff = 0;

    /* Cell sizes for static load balancing, first index cartesian */
    std::array<std::vector<real>, DIM> slb_frac;

    /**< Information about the simulated system */
    DDSystemInfo systemInfo;

    /* The width of the communicated boundaries */
    /**< Cut-off for multi-body interactions, also 2-body bonded when \p cutoff_mody > \p cutoff */
    real cutoff_mbody = 0;
    /**< The minimum guaranteed cell-size, Cartesian indexing */
    gmx::RVec cellsize_min = { 0, 0, 0 };
    /**< The minimum guaranteed cell-size with dlb=auto */
    gmx::RVec cellsize_min_dlb = { 0, 0, 0 };
    /**< The lower limit for the DD cell size with DLB */
    real cellsize_limit = 0;
    /**< Effectively no NB cut-off limit with DLB for systems without PBC? */
    bool bVacDLBNoLimit = false;

    /** With PME load balancing we set limits on DLB */
    bool bPMELoadBalDLBLimits = false;
    /** DLB needs to take into account that we want to allow this maximum
     *  cut-off (for PME load balancing), this could limit cell boundaries.
     */
    real PMELoadBal_max_cutoff = 0;

    /**< box lower corner, required with dim's without pbc and -gcom */
    gmx::RVec box0 = { 0, 0, 0 };
    /**< box size, required with dim's without pbc and -gcom */
    gmx::RVec box_size = { 0, 0, 0 };

    /**< The DD cell lower corner, in triclinic space */
    gmx::RVec cell_x0 = { 0, 0, 0 };
    /**< The DD cell upper corner, in triclinic space */
    gmx::RVec cell_x1 = { 0, 0, 0 };

    /**< The old \p cell_x0, to check cg displacements */
    gmx::RVec old_cell_x0 = { 0, 0, 0 };
    /**< The old \p cell_x1, to check cg displacements */
    gmx::RVec old_cell_x1 = { 0, 0, 0 };

    /** The communication setup and charge group boundaries for the zones */
    gmx_domdec_zones_t zones;

    /* The zone limits for DD dimensions 1 and 2 (not 0), determined from
     * cell boundaries of neighboring cells for staggered grids when using
     * dynamic load balancing.
     */
    /**< Zone limits for dim 1 with staggered grids */
    std::array<gmx_ddzone_t, 2> zone_d1;
    /**< Zone limits for dim 2 with staggered grids */
    gmx_ddzone_t zone_d2[2][2];

    /** The coordinate/force communication setup and indices */
    std::array<gmx_domdec_comm_dim_t, DIM> cd;
    /** Restricts the maximum number of cells to communicate with in one dimension
     *
     * Dynamic load balancing is not permitted to change sizes if it
     * would violate this restriction. */
    int maxpulse = 0;

    /** The step interval for algorithms that require global communication
     *  such as DLB and the computation the extent of unbound dimensions
     * (i.e. dimensions without PBC and without walls).
     */
    int nstDDGlobalComm = 0;

    /** Which cg distribution is stored on the main node,
     *  stored as DD partitioning call count.
     */
    int64_t main_cg_ddp_count = 0;

    /** The number of cg's received from the direct neighbors */
    std::array<int, DD_MAXZONE> zone_ncg1 = { 0 };

    /** The atom ranges in the local state */
    DDAtomRanges atomRanges;

    /** Array for signalling if atoms have moved to another domain */
    std::vector<int> movedBuffer;

    /** Communication int buffer for general use */
    DDBuffer<int> intBuffer;

    /** Communication rvec buffer for general use */
    DDBuffer<gmx::RVec> rvecBuffer;

    /* Temporary storage for thread parallel communication setup */
    /**< Thread-local work data */
    std::vector<dd_comm_setup_work_t> dth;

    /* Communication buffer only used with multiple grid pulses */
    /**< Another rvec comm. buffer */
    DDBuffer<gmx::RVec> rvecBuffer2;

    /* Communication buffers for local redistribution */
    /**< Charge group flag comm. buffers */
    std::array<std::vector<int>, DIM * 2> cggl_flag;
    /**< Charge group center comm. buffers */
    std::array<std::vector<gmx::RVec>, DIM * 2> cgcm_state;

    /* Cell sizes for dynamic load balancing */
    std::vector<DDCellsizesWithDlb> cellsizesWithDlb;

    /* Stuff for load communication */
    /**< The recorded load data */
    std::vector<domdec_load_t> load;
    /**< The number of MPI ranks sharing the GPU our rank is using */
    int nrank_gpu_shared = 0;
#if GMX_MPI
    /**< The MPI load communicator */
    std::vector<MPI_Comm> mpi_comm_load;
    /**< The MPI load communicator for ranks sharing a GPU */
    MPI_Comm mpi_comm_gpu_shared;
#endif

    /**< Struct for timing the force load balancing region */
    BalanceRegion balanceRegion;

#if GMX_MPI
    /**< MPI data type corresponding to rvec */
    MPI_Datatype mpiRVec;
#endif

    /* Cycle counters over nstlist steps */
    /**< Total cycles counted */
    std::array<float, ddCyclNr> cycl = { 0 };
    /**< The number of cycle recordings */
    std::array<int, ddCyclNr> cycl_n = { 0 };
    /**< The maximum cycle count */
    std::array<float, ddCyclNr> cycl_max = { 0 };
    /**< Total flops counted */
    double flop = 0.0;
    /**< The number of flop recordings */
    int flop_n = 0;
    /** How many times did we have load measurements */
    int n_load_have = 0;
    /** How many times have we collected the load measurements */
    int n_load_collect = 0;

    /* Cycle count history for DLB checks */
    /**< The averaged cycles per step over the last nstlist step before turning on DLB */
    float cyclesPerStepBeforeDLB = 0;
    /**< The running average of the cycles per step during DLB */
    float cyclesPerStepDlbExpAverage = 0;
    /**< Have we turned off DLB (after turning DLB on)? */
    bool haveTurnedOffDlb = false;
    /**< The DD step at which we last measured that DLB off was faster than DLB on, 0 if there was no such step */
    int64_t dlbSlowerPartitioningCount = 0;

    /* Statistics for atoms */
    /**< The atoms per range, summed over the steps */
    double sum_nat[static_cast<int>(DDAtomRanges::Type::Number)] = {};

    /* Statistics for calls and times */
    /**< The number of partioning calls */
    int ndecomp = 0;
    /**< The number of load recordings */
    int nload = 0;
    /**< Total MD step time */
    double load_step = 0.0;
    /**< Total PP force time */
    double load_sum = 0.0;
    /**< Max \p load_sum over the ranks */
    double load_max = 0.0;
    /**< Was load balancing limited, per DD dim */
    gmx::IVec load_lim = { 0, 0, 0 };
    /**< Total time on PP done during PME overlap time */
    double load_mdf = 0.0;
    /**< Total time on our PME-only rank */
    double load_pme = 0.0;

    /** The last partition step */
    int64_t partition_step = INT_MIN;
};

/*! \brief DD zone permutation
 *
 * Zone permutation from the Cartesian x-major/z-minor order to an order
 * that leads to consecutive charge groups for neighbor searching.
 * TODO: It should be possible to remove this now that the group scheme is removed
 */
static const int zone_perm[3][4] = { { 0, 0, 0, 0 }, { 1, 0, 0, 0 }, { 3, 0, 1, 2 } };

/* dd_zo and dd_zp3 is set up such that i zones with non-zero
 * components see only j zones with that component 0.
 */

/*! \brief Returns the DD cut-off distance for multi-body interactions */
real dd_cutoff_multibody(const gmx_domdec_t* dd);

/*! \brief Returns the DD cut-off distance for two-body interactions */
real dd_cutoff_twobody(const gmx_domdec_t* dd);

/*! \brief Returns the domain index given the number of domains and the domain coordinates
 *
 * This order is required to minimize the coordinate communication in PME
 * which uses decomposition in the x direction.
 */
static inline int dd_index(const ivec numDomains, const ivec domainCoordinates)
{
    return ((domainCoordinates[XX] * numDomains[YY] + domainCoordinates[YY]) * numDomains[ZZ])
           + domainCoordinates[ZZ];
};

/*! Returns the size of the buffer to hold fractional cell boundaries for DD dimension index dimIndex */
static inline int ddCellFractionBufferSize(const gmx_domdec_t* dd, int dimIndex)
{
    return dd->numCells[dd->dim[dimIndex]] + 1 + dimIndex * 2 + 1 + dimIndex;
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
static constexpr double DD_CELL_MARGIN = 1.0001;

/*! \brief Factor for checking DD cell size limitation during DLB, should be in between 1 and DD_CELL_MARGIN */
static constexpr double DD_CELL_MARGIN2 = 1.00005;

/*! \brief With pressure scaling, keep cell sizes 2% above the limit to allow for some scaling */
static constexpr double DD_PRES_SCALE_MARGIN = 1.02;

/*! \endcond */

#endif
