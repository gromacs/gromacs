/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2005- The GROMACS Authors
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
/*! \defgroup module_domdec Spatial domain decomposition (for parallelization over MPI)
 * \ingroup group_mdrun
 *
 * \brief Manages the decomposition of the simulation volume over MPI
 * ranks to try to distribute work evenly with minimal communication
 * overheads.
 *
 * \todo Get domdec stuff out of mdtypes/commrec.h
 *
 * \author Berk Hess <hess@kth.se>
 *
 */

/*! \libinternal \file
 *
 * \brief This file declares functions for mdrun to call to manage the
 * details of its domain decomposition.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_H
#define GMX_DOMDEC_DOMDEC_H

#include <cstddef>

#include <vector>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_ddbox_t;
struct gmx_domdec_zones_t;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_wallcycle;
enum class PbcType : int;
class t_state;
class DeviceContext;
class GpuEventSynchronizer;

namespace gmx
{
struct AtomInfoWithinMoleculeBlock;
class DeviceStreamManager;
class DomdecZones;
class ForceWithShiftForces;
class MDLogger;
class RangePartitioning;
class VirtualSitesHandler;
template<typename>
class ArrayRef;
template<typename, size_t>
class FixedCapacityVector;
} // namespace gmx

/*! \brief Returns the global topology atom number belonging to local atom index i.
 *
 * This function is intended for writing ASCII output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */
int ddglatnr(const gmx_domdec_t* dd, int i);

/*! \brief Store the global cg indices of the home cgs in state,
 *
 * This means it can be reset, even after a new DD partitioning.
 */
void dd_store_state(const gmx_domdec_t& dd, t_state* state);

/*! \brief Returns a const reference to the gmx::DomdecZones object */
const gmx::DomdecZones& getDomdecZones(const gmx_domdec_t& dd);

/*! \brief Returns the range for atoms in zones*/
int dd_numAtomsZones(const gmx_domdec_t& dd);

/*! \brief Returns the number of home atoms */
int dd_numHomeAtoms(const gmx_domdec_t& dd);

/*! \brief Returns the atom range in the local state for atoms that need to be present in mdatoms */
int dd_natoms_mdatoms(const gmx_domdec_t& dd);

/*! \brief Returns the atom range in the local state for atoms involved in virtual sites */
int dd_natoms_vsite(const gmx_domdec_t& dd);

/*! \brief Sets the atom range for atom in the local state for atoms received in constraints communication */
void dd_get_constraint_range(const gmx_domdec_t& dd, int* at_start, int* at_end);

/*! \libinternal \brief Struct for passing around the number of PME domains */
struct NumPmeDomains
{
    int x; //!< The number of PME domains along dimension x
    int y; //!< The number of PME domains along dimension y
};

/*! \brief Returns the number of PME domains, can be called with dd=NULL */
NumPmeDomains getNumPmeDomains(const gmx_domdec_t* dd);

/*! \brief Returns the set of DD ranks that communicate with pme node cr->nodeid */
std::vector<int> get_pme_ddranks(const t_commrec* cr, int pmenodeid);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim x */
int dd_pme_maxshift_x(const gmx_domdec_t& dd);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim y */
int dd_pme_maxshift_y(const gmx_domdec_t& dd);

/*! \brief Return whether update groups are used */
bool ddUsesUpdateGroups(const gmx_domdec_t& dd);

/*! \brief Returns whether molecules are always whole, i.e. not broken by PBC */
bool dd_moleculesAreAlwaysWhole(const gmx_domdec_t& dd);

/*! \brief Returns if we need to do pbc for calculating bonded interactions */
bool dd_bonded_molpbc(const gmx_domdec_t& dd, PbcType pbcType);

/*! \brief Change the DD non-bonded communication cut-off.
 *
 * This could fail when trying to increase the cut-off,
 * then FALSE will be returned and the cut-off is not modified.
 *
 * \param[in] cr               Communication recrod
 * \param[in] box              Box matrix, used for computing the dimensions of the system
 * \param[in] x                Position vector, used for computing the dimensions of the system
 * \param[in] cutoffRequested  The requested atom to atom cut-off distance, usually the pair-list
 *                             cutoff distance
 * \param[in] checkGpuDdLimitation Whether to check the GPU DD support limitation
 */
bool change_dd_cutoff(t_commrec*                     cr,
                      const matrix                   box,
                      gmx::ArrayRef<const gmx::RVec> x,
                      real                           cutoffRequested,
                      bool                           checkGpuDdLimitation);

/*! \brief Set up communication for averaging GPU wait times over domains
 *
 * When domains (PP MPI ranks) share a GPU, the individual GPU wait times
 * are meaningless, as it depends on the order in which tasks on the same
 * GPU finish. Therefore there wait times need to be averaged over the ranks
 * sharing the same GPU. This function sets up the communication for that.
 */
void dd_setup_dlb_resource_sharing(const t_commrec* cr, int gpu_id);

/*! \brief Cycle counter indices used internally in the domain decomposition */
enum
{
    ddCyclStep,
    ddCyclPPduringPME,
    ddCyclF,
    ddCyclWaitGPU,
    ddCyclPME,
    ddCyclNr
};

/*! \brief Add the wallcycle count to the DD counter */
void dd_cycles_add(const gmx_domdec_t* dd, float cycles, int ddCycl);

/*! \brief Communicate the coordinates to the neighboring cells and do pbc. */
void dd_move_x(struct gmx_domdec_t* dd, const matrix box, gmx::ArrayRef<gmx::RVec> x, gmx_wallcycle* wcycle);

/*! \brief Sum the forces over the neighboring cells.
 *
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */
void dd_move_f(struct gmx_domdec_t* dd, gmx::ForceWithShiftForces* forceWithShiftForces, gmx_wallcycle* wcycle);

/*! \brief Reset all the statistics and counters for total run counting */
void reset_dd_statistics_counters(struct gmx_domdec_t* dd);

/* In domdec_con.c */

/*! \brief Communicates the virtual site forces, reduces the shift forces when \p fshift != NULL */
void dd_move_f_vsites(const gmx_domdec_t& dd, gmx::ArrayRef<gmx::RVec> f, gmx::ArrayRef<gmx::RVec> fshift);

/*! \brief Clears the forces for virtual sites */
void dd_clear_f_vsites(const gmx_domdec_t& dd, gmx::ArrayRef<gmx::RVec> f);

/*! \brief Move x0 and also x1 if x1!=NULL. bX1IsCoord tells if to do PBC on x1 */
void dd_move_x_constraints(struct gmx_domdec_t*     dd,
                           const matrix             box,
                           gmx::ArrayRef<gmx::RVec> x0,
                           gmx::ArrayRef<gmx::RVec> x1,
                           bool                     bX1IsCoord);

/*! \brief Communicates the coordinates involved in virtual sites */
void dd_move_x_vsites(const gmx_domdec_t& dd, const matrix box, gmx::ArrayRef<gmx::RVec> x);
/*! \brief Communicates the positions and velocities involved in virtual sites */
void dd_move_x_and_v_vsites(const gmx_domdec_t&      dd,
                            const matrix             box,
                            gmx::ArrayRef<gmx::RVec> x,
                            gmx::ArrayRef<gmx::RVec> v);

/*! \brief Returns the local atom count array for all constraints
 *
 * The local atom count for a constraint, possible values 2/1/0, is needed
 * to avoid not/double-counting contributions linked to the Lagrange
 * multiplier, such as the virial and free-energy derivatives.
 *
 * \note When \p dd = nullptr, an empty reference is returned.
 */
gmx::ArrayRef<const int> dd_constraints_nlocalatoms(const gmx_domdec_t* dd);

/*! \brief Construct local state */
void dd_init_local_state(const gmx_domdec_t& dd, const t_state* state_global, t_state* local_state);

/*! \brief Construct the GPU halo exchange object(s).
 *
 * \param[in] cr                  The commrec object.
 * \param[in] deviceStreamManager Manager of the GPU context and streams.
 * \param[in] wcycle              The wallclock counter.
 */
void constructGpuHaloExchange(const t_commrec&                cr,
                              const gmx::DeviceStreamManager& deviceStreamManager,
                              gmx_wallcycle*                  wcycle);

/*! \brief
 * (Re-) Initialization for GPU halo exchange
 * \param [in] cr                   The commrec object
 * \param [in] d_coordinatesBuffer  pointer to coordinates buffer in GPU memory
 * \param [in] d_forcesBuffer       pointer to forces buffer in GPU memory
 */
void reinitGpuHaloExchange(const t_commrec&        cr,
                           DeviceBuffer<gmx::RVec> d_coordinatesBuffer,
                           DeviceBuffer<gmx::RVec> d_forcesBuffer);


/*! \brief GPU halo exchange of coordinates buffer.
 * \param [in] cr                The commrec object
 * \param [in] box               Coordinate box (from which shifts will be constructed)
 * \param [in] dependencyEvent   Dependency event for this operation
 * \returns                      Event recorded when this operation has been launched
 */
GpuEventSynchronizer* communicateGpuHaloCoordinates(const t_commrec&      cr,
                                                    const matrix          box,
                                                    GpuEventSynchronizer* dependencyEvent);

/*! \brief  Wait for copy of nonlocal part of coordinate array from GPU to CPU
 * following coordinate halo exchange
 * \param [in] cr   The commrec object
 * \param [in] accumulateForces  True if forces should accumulate, otherwise they are set
 * \param [in] dependencyEvents  Dependency events for this operation
 */
void communicateGpuHaloForces(const t_commrec&                                    cr,
                              bool                                                accumulateForces,
                              gmx::FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents);

/*! \brief Wraps the \c positions so that atoms from the same
 * update group share the same periodic image wrt \c box.
 *
 * When DD and update groups are in use, the simulation main rank
 * should call this to ensure that e.g. when restarting a simulation
 * that did not use update groups that the coordinates satisfy the new
 * requirements.
 *
 * This function can probably be removed when even single-rank
 * simulations use domain decomposition, because then the choice of
 * whether update groups are used is probably going to be the same
 * regardless of the rank count.
 *
 * \param[in]    dd         The DD manager
 * \param[in]    mtop       The system topology
 * \param[in]    box        The global system box
 * \param[in]    positions  The global system positions
 */
void putUpdateGroupAtomsInSamePeriodicImage(const gmx_domdec_t&      dd,
                                            const gmx_mtop_t&        mtop,
                                            const matrix             box,
                                            gmx::ArrayRef<gmx::RVec> positions);

#endif
