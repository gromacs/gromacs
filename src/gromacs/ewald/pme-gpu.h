/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * \brief This file contains function definitions for performing the PME calculations on GPU.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef PMEGPU_H
#define PMEGPU_H

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/timing/gpu_timing.h"
#include "pme-gpu-types.h"
#include "pme-internal.h"

struct gmx_hw_info_t;
struct gmx_gpu_opt_t;

/* Internal data handling */

/*! \brief \internal
 * Copies the forces from the CPU buffer (pme->gpu->forcesHost) to the GPU
 * (to reduce them with the PME GPU gathered forces).
 * To be called after the bonded calculations.
 * Does nothing on non-CUDA builds.
 *
 * \param[in] pme            The PME structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_forces(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme)) CUDA_FUNC_TERM


/*! \brief \internal
 * Copies the input coordinates from the CPU buffer (pme->gpu->coordinatesHost) onto the GPU.
 *
 * \param[in] pme            The PME structure.
 *
 * Needs to be called every MD step. The coordinates are then used in the spline calculation.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_coordinates(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme));


// nice external functions

/*! \brief \internal
 * Finds out if PME is set to run on GPU.
 *
 * \param[in] pme  The PME structure.
 * \returns        TRUE if PME runs on GPU, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_enabled(const gmx_pme_t *pme)
{
    /* Something to think about: should this function be called from all the CUDA_FUNC_QUALIFIER functions?
     * In other words, should we plan for dynamic toggling of the PME GPU?
     */
    return (pme != NULL) && pme->bGPU;
}

/*! \brief \internal
 * (Re-)initializes the PME GPU data at the beginning of the run or on DLB. Does nothing on non-CUDA builds.
 *
 * \param[in] pme     The PME structure.
 * \param[in] hwinfo  The hardware information structure.
 * \param[in] gpu_opt The GPU information structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reinit(gmx_pme_t           *CUDA_FUNC_ARGUMENT(pme),
                                        const gmx_hw_info_t *CUDA_FUNC_ARGUMENT(hwinfo),
                                        const gmx_gpu_opt_t *CUDA_FUNC_ARGUMENT(gpu_opt)) CUDA_FUNC_TERM

/*! \brief \internal
 * Destroys the PME GPU data at the end of the run. Does nothing on non-CUDA builds.
 *
 * \param[in] pme     The PME structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy(gmx_pme_t *CUDA_FUNC_ARGUMENT(pme)) CUDA_FUNC_TERM

/*! \brief
 * Starts the PME GPU step (copies coordinates onto GPU, possibly sets the unit cell parameters).
 * Does nothing if PME on GPU is disabled.
 *
 * \param[in] pme     The PME structure.
 * \param[in] box     The unit cell box which does not necessarily change every step (only with pressure coupling enabled).
 *                    Currently it is simply compared with the previous one to determine if it needs to be updated.
 */
void pme_gpu_start_step(const gmx_pme_t *pme, const matrix box);

/*! \brief \internal
 * Finishes the PME GPU step, waiting for the output forces and/or energy/virial to be copied to the host. Does nothing on non-CUDA builds.
 *
 * \param[in] pme            The PME structure.
 * \param[in] bCalcForces    The left-over flag from the CPU code which tells the function to copy the forces to the CPU side. Should be passed to the launch call instead.
 * \param[in] bCalcEnerVir   The left-over flag from the CPU code which tells the function to copy the energy/virial to the CPU side. Should be passed to the launch call instead.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_finish_step(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme),
                                             const gmx_bool   CUDA_FUNC_ARGUMENT(bCalcForces),
                                             const gmx_bool   CUDA_FUNC_ARGUMENT(bCalcEnerVir)) CUDA_FUNC_TERM

/*! \brief \internal
 * Gets the PME GPU output virial/energy. Should be called after pme_gpu_finish_step. Does nothing on non-CUDA builds.
 *
 * \param[in]  pme  The PME structure.
 * \param[out] energy  The output energy pointer.
 * \param[out] virial  The output virial matrix.
 *
 * Should thsi be merged with pme_gpu_finish_step?
 */
void pme_gpu_get_energy_virial(const gmx_pme_t *pme, real *energy, matrix virial);

/*! \brief \internal
 * Reallocates the local atoms data (charges, coordinates, etc.). Copies the charges. Does nothing on non-CUDA builds.
 *
 * \param[in] pme            The PME structure.
 * \param[in] nAtoms         The number of particles.
 * \param[in] coefficients   The pointer to the host-side array of particle charges.
 *
 * This is a function that should only be called in the beginning of the run and on domain decomposition.
 * Should be called before the pme_gpu_set_io_ranges.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reinit_atoms(const gmx_pme_t  *CUDA_FUNC_ARGUMENT(pme),
                                              const int         CUDA_FUNC_ARGUMENT(nAtoms),
                                              float            *CUDA_FUNC_ARGUMENT(coefficients)) CUDA_FUNC_TERM

/*! \brief \internal
 * Sets the host-side I/O buffers in the PME GPU.
 * Does nothing if PME on GPU is disabled.
 *
 * \param[in] pme            The PME structure.
 * \param[in] coordinates    The pointer to the host-side array of particle coordinates in rvec format.
 * \param[in] forces         The pointer to the host-side array of particle forces.
 *                           It will be used for output, but can also be used for input,
 *                           if bClearForces is passed as false to the pme_gpu_launch_gather.
 */
void pme_gpu_set_io_ranges(const gmx_pme_t *pme, rvec *coordinates, rvec *forces);

/*! \brief \internal
 * Resets the PME GPU timings. To be called at the reset step. Does nothing on non-CUDA builds.
 *
 * \param[in] pme            The PME structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reset_timings(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme)) CUDA_FUNC_TERM

/*! \brief \internal
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end. Does nothing on non-CUDA builds.
 *
 * \param[in] pme               The PME structure
 * \param[in] timings           The gmx_wallclock_gpu_t structure (with some shamelessly duplicated fields for the PME GPU timings).
 */
CUDA_FUNC_QUALIFIER void pme_gpu_get_timings(const gmx_pme_t      *CUDA_FUNC_ARGUMENT(pme),
                                             gmx_wallclock_gpu_t **CUDA_FUNC_ARGUMENT(timings)) CUDA_FUNC_TERM


/* GPU framework-agnostic functions follow */

/* The convenience PME GPU status getters */

/*! \brief \internal
 * Tells if PME performs the gathering stage on GPU.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if the gathering is performed on GPU, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_performs_gather(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && pme->gpu->bGPUGather;
}

/*! \brief \internal
 * Tells if PME performs the FFT stages on GPU.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if FFT is performed on GPU, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_performs_FFT(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && pme->gpu->bGPUFFT;
}

/*! \brief \internal
 * Tells if PME performs the grid (un-)wrapping on GPU.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if (un-)wrapping is performed on GPU, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_performs_wrapping(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && pme->gpu->bGPUSingle;
}

/*! \brief \internal
 * Tells if PME performs the grid solving on GPU.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if solving is performed on GPU, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_performs_solve(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && pme->gpu->bGPUSolve;
}

/*! \brief \internal
 * Tells if PME runs on multiple GPUs.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if PME runs on multiple GPUs, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_uses_dd(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && !pme->gpu->bGPUSingle;
}

/* The main PME GPU functions (separate stages and the whole loop) will live here*/

/*! \brief \internal
 * Gets the output forces and virial/energy if corresponding flags are (were?) passed in.
 *
 * \param[in]  pme            The PME data structure.
 * \param[in]  wcycle         The wallclock counter.
 * \param[out] vir_q          The output virial matrix.
 * \param[out] energy_q       The output energy.
 * \param[in]  flags          The combination of flags to affect the output.
 *                            Pass GMX_PME_CALC_ENER_VIR to get the virial and energy.
 *                            GMX_PME_CALC_F should be affecting the force output,
 *                            but likely will not as the force copy has already been scheduled before.
 *                            TODO: rethink the flag handling.
 */
void pme_gpu_get_results(const gmx_pme_t *pme,
                         gmx_wallcycle_t  wcycle,
                         matrix           vir_q,
                         real            *energy_q,
                         int              flags);

#endif // PMEGPU_H
