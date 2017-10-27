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
/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_H
#define GMX_EWALD_PME_H

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct interaction_const_t;
struct t_commrec;
struct t_inputrec;
struct t_nrnb;
struct PmeGpu;
struct gmx_wallclock_gpu_pme_t;
struct gmx_device_info_t;
struct gmx_pme_t;

enum class GpuTaskCompletion;

namespace gmx
{
class ForceWithVirial;
class MDLogger;
}

enum {
    GMX_SUM_GRID_FORWARD, GMX_SUM_GRID_BACKWARD
};

/*! \brief Possible PME codepaths on a rank.
 * \todo: make this enum class with gmx_pme_t C++ refactoring
 */
enum PmeRunMode
{
    None,    //!< No PME task is done
    CPU,     //!< Whole PME computation is done on CPU
    GPU,     //!< Whole PME computation is done on GPU
    Mixed,   //!< Mixed mode: only spread and gather run on GPU; FFT and solving are done on CPU.
};

//! PME gathering output forces treatment
enum class PmeForceOutputHandling
{
    Set,             /**< Gather simply writes into provided force buffer */
    ReduceWithInput, /**< Gather adds its output to the buffer.
                        On GPU, that means additional H2D copy before the kernel launch. */
};

/*! \brief Return the smallest allowed PME grid size for \p pmeOrder */
int minimalPmeGridSize(int pmeOrder);

/*! \brief Check restrictions on pme_order and the PME grid nkx,nky,nkz.
 *
 * With errorsAreFatal=true, an exception or fatal error is generated
 * on violation of restrictions.
 * With errorsAreFatal=false, false is returned on violation of restrictions.
 * When all restrictions are obeyed, true is returned.
 * Argument useThreads tells if any MPI rank doing PME uses more than 1 threads.
 * If at calling useThreads is unknown, pass true for conservative checking.
 *
 * The PME GPU restrictions are checked separately during pme_gpu_init().
 */
bool gmx_pme_check_restrictions(int pme_order,
                                int nkx, int nky, int nkz,
                                int nnodes_major,
                                bool useThreads,
                                bool errorsAreFatal);

/*! \brief Construct PME data
 *
 * \throws   gmx::InconsistentInputError if input grid sizes/PME order are inconsistent.
 * \returns  Pointer to newly allocated and initialized PME data.
 */
gmx_pme_t *gmx_pme_init(const t_commrec *cr,
                        int nnodes_major, int nnodes_minor,
                        const t_inputrec *ir, int homenr,
                        gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                        gmx_bool bReproducible,
                        real ewaldcoeff_q, real ewaldcoeff_lj,
                        int nthread,
                        PmeRunMode runMode,
                        PmeGpu *pmeGpu,
                        gmx_device_info_t *gpuInfo,
                        const gmx::MDLogger &mdlog);

/*! \brief Destroys the PME data structure.*/
void gmx_pme_destroy(gmx_pme_t *pme);

//@{
/*! \brief Flag values that control what gmx_pme_do() will calculate
 *
 * These can be combined with bitwise-OR if more than one thing is required.
 */
#define GMX_PME_SPREAD        (1<<0)
#define GMX_PME_SOLVE         (1<<1)
#define GMX_PME_CALC_F        (1<<2)
#define GMX_PME_CALC_ENER_VIR (1<<3)
/* This forces the grid to be backtransformed even without GMX_PME_CALC_F */
#define GMX_PME_CALC_POT      (1<<4)

#define GMX_PME_DO_ALL_F  (GMX_PME_SPREAD | GMX_PME_SOLVE | GMX_PME_CALC_F)
//@}

/*! \brief Do a PME calculation on a CPU for the long range electrostatics and/or LJ.
 *
 * The meaning of \p flags is defined above, and determines which
 * parts of the calculation are performed.
 *
 * \return 0 indicates all well, non zero is an error code.
 */
int gmx_pme_do(struct gmx_pme_t *pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real chargeA[],  real chargeB[],
               real c6A[],      real c6B[],
               real sigmaA[],   real sigmaB[],
               matrix box,      t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix vir_q,    matrix vir_lj,
               real *energy_q,  real *energy_lj,
               real lambda_q,   real lambda_lj,
               real *dvdlambda_q, real *dvdlambda_lj,
               int flags);

/*! \brief Called on the nodes that do PME exclusively (as slaves) */
int gmx_pmeonly(struct gmx_pme_t *pme,
                struct t_commrec *cr,     t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                t_inputrec *ir, PmeRunMode runMode);

/*! \brief Calculate the PME grid energy V for n charges.
 *
 * The potential (found in \p pme) must have been found already with a
 * call to gmx_pme_do() with at least GMX_PME_SPREAD and GMX_PME_SOLVE
 * specified. Note that the charges are not spread on the grid in the
 * pme struct. Currently does not work in parallel or with free
 * energy.
 */
void gmx_pme_calc_energy(struct gmx_pme_t *pme, int n, rvec *x, real *q, real *V);

/*! \brief Send the charges and maxshift to out PME-only node. */
void gmx_pme_send_parameters(struct t_commrec *cr,
                             const interaction_const_t *ic,
                             gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                             real *chargeA, real *chargeB,
                             real *sqrt_c6A, real *sqrt_c6B,
                             real *sigmaA, real *sigmaB,
                             int maxshift_x, int maxshift_y);

/*! \brief Send the coordinates to our PME-only node and request a PME calculation */
void gmx_pme_send_coordinates(struct t_commrec *cr, matrix box, rvec *x,
                              real lambda_q, real lambda_lj,
                              gmx_bool bEnerVir,
                              gmx_int64_t step);

/*! \brief Tell our PME-only node to finish */
void gmx_pme_send_finish(struct t_commrec *cr);

/*! \brief Tell our PME-only node to reset all cycle and flop counters */
void gmx_pme_send_resetcounters(struct t_commrec *cr, gmx_int64_t step);

/*! \brief PP nodes receive the long range forces from the PME nodes */
void gmx_pme_receive_f(struct t_commrec *cr,
                       gmx::ForceWithVirial *forceWithVirial,
                       real *energy_q, real *energy_lj,
                       real *dvdlambda_q, real *dvdlambda_lj,
                       float *pme_cycles);

/*! \brief
 * This function updates the local atom data on GPU after DD (charges, coordinates, etc.).
 * TODO: it should update the PME CPU atom data as well.
 * (currently PME CPU call gmx_pme_do() gets passed the input pointers for each computation).
 *
 * \param[in] pme            The PME structure.
 * \param[in] nAtoms         The number of particles.
 * \param[in] charges        The pointer to the array of particle charges.
 */
void gmx_pme_reinit_atoms(const gmx_pme_t *pme, const int nAtoms, const real *charges);

/* A block of PME GPU functions */

/*! \brief Checks whether the input system allows to run PME on GPU.
 * TODO: this mostly duplicates an internal PME assert function
 * pme_gpu_check_restrictions(), except that works with a
 * formed gmx_pme_t structure. Should that one go away/work with inputrec?
 *
 * \param[in]  ir     Input system.
 * \param[out] error  The error message if the input is not supported on GPU.
 *
 * \returns true if PME can run on GPU with this input, false otherwise.
 */
bool pme_gpu_supports_input(const t_inputrec *ir, std::string *error);

/*! \brief
 * Returns the active PME codepath (CPU, GPU, mixed).
 * \todo This is a rather static data that should be managed by the higher level task scheduler.
 *
 * \param[in]  pme            The PME data structure.
 * \returns active PME codepath.
 */
PmeRunMode pme_run_mode(const gmx_pme_t *pme);

/*! \brief
 * Tells if PME is enabled to run on GPU (not necessarily active at the moment).
 * \todo This is a rather static data that should be managed by the hardware assignment manager.
 * For now, it is synonymous with the active PME codepath (in the absence of dynamic switching).
 *
 * \param[in]  pme            The PME data structure.
 * \returns true if PME can run on GPU, false otherwise.
 */
inline bool pme_gpu_task_enabled(const gmx_pme_t *pme)
{
    return (pme != nullptr) && (pme_run_mode(pme) != PmeRunMode::CPU);
}

/*! \brief
 * Resets the PME GPU timings. To be called at the reset step.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_reset_timings(const gmx_pme_t *pme);

/*! \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_pme_t structure (for log output). To be called at the run end.
 *
 * \param[in] pme               The PME structure.
 * \param[in] timings           The gmx_wallclock_gpu_pme_t structure.
 */
void pme_gpu_get_timings(const gmx_pme_t         *pme,
                         gmx_wallclock_gpu_pme_t *timings);

/* The main PME GPU functions */

/*! \brief
 * Prepares PME on GPU computation (updating the box if needed)
 * \param[in] pme               The PME data structure.
 * \param[in] needToUpdateBox   Tells if the stored unit cell parameters should be updated from \p box.
 * \param[in] box               The unit cell box.
 * \param[in] wcycle            The wallclock counter.
 * \param[in] flags             The combination of flags to affect this PME computation.
 *                              The flags are the GMX_PME_ flags from pme.h.
 */
void pme_gpu_prepare_computation(gmx_pme_t      *pme,
                                 bool            needToUpdateBox,
                                 const matrix    box,
                                 gmx_wallcycle_t wcycle,
                                 int             flags);

/*! \brief
 * Launches first stage of PME on GPU - H2D input transfers, spreading kernel, and D2H grid transfer if needed.
 *
 * \param[in] pme               The PME data structure.
 * \param[in] x                 The array of local atoms' coordinates.
 * \param[in] wcycle            The wallclock counter.
 */
void pme_gpu_launch_spread(gmx_pme_t      *pme,
                           const rvec     *x,
                           gmx_wallcycle_t wcycle);

/*! \brief
 * Launches middle stages of PME (FFT R2C, solving, FFT C2R) either on GPU or on CPU, depending on the run mode.
 *
 * \param[in] pme               The PME data structure.
 * \param[in] wcycle            The wallclock counter.
 */
void pme_gpu_launch_complex_transforms(gmx_pme_t       *pme,
                                       gmx_wallcycle_t  wcycle);

/*! \brief
 * Launches last stage of PME on GPU - force gathering and D2H force transfer.
 *
 * \param[in]  pme               The PME data structure.
 * \param[in]  wcycle            The wallclock counter.
 * \param[in]  forceTreatment    Tells how data should be treated. The gathering kernel either stores
 *                               the output reciprocal forces into the host array, or copies its contents to the GPU first
 *                               and accumulates. The reduction is non-atomic.
 */
void pme_gpu_launch_gather(const gmx_pme_t        *pme,
                           gmx_wallcycle_t         wcycle,
                           PmeForceOutputHandling  forceTreatment);

/*! \brief
 * Blocks until PME GPU tasks are completed, and gets the output forces and virial/energy
 * (if they were to be computed).
 *
 * \param[in]  pme            The PME data structure.
 * \param[out] wcycle         The wallclock counter.
 * \param[out] forces         The output forces.
 * \param[out] virial         The output virial matrix.
 * \param[out] energy         The output energy.
 */
void pme_gpu_wait_finish_task(const gmx_pme_t                *pme,
                              gmx_wallcycle_t                 wcycle,
                              gmx::ArrayRef<const gmx::RVec> *forces,
                              matrix                          virial,
                              real                           *energy);
/*! \brief
 * Attempts to complete PME GPU tasks.
 *
 * The \p completionKind argument controls whether the function blocks until all
 * PME GPU tasks enqueued completed (as pme_gpu_wait_finish_task() does) or only
 * checks and returns immediately if they did not.
 * When blocking or the tasks have completed it also gets the output forces
 * by assigning the ArrayRef to the \p forces pointer passed in.
 * Virial/energy are also outputs if they were to be computed.
 *
 * Note: also launches the reinitalization of the PME output buffers.
 * TODO: this should be moved out to avoid miscounting its wall-time (as wait iso launch).
 *
 * \param[in]  pme            The PME data structure.
 * \param[in]  wcycle         The wallclock counter.
 * \param[out] forces         The output forces.
 * \param[out] virial         The output virial matrix.
 * \param[out] energy         The output energy.
 * \param[in]  completionKind  Indicates whether PME task completion should only be checked rather than waited for
 * \returns                   True if the PME GPU tasks have completed
 */
bool pme_gpu_try_finish_task(const gmx_pme_t                *pme,
                             gmx_wallcycle_t                 wcycle,
                             gmx::ArrayRef<const gmx::RVec> *forces,
                             matrix                          virial,
                             real                           *energy,
                             GpuTaskCompletion               completionKind);


#endif
