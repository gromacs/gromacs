/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

// FIXME: remove the "__" prefix in front of the group def when we move the
//        nonbonded code into separate dir.

/*! \libinternal \defgroup __module_nb_verlet Short-range non-bonded interaction module
 * \ingroup group_mdrun
 *
 * \brief Computes forces and energies for short-range pair-interactions
 * based on the Verlet algorithm. The algorithm uses pair-lists generated
 * at fixed intervals as well as various flavors of pair interaction kernels
 * implemented for a wide range of CPU and GPU architectures.
 *
 * The module includes support for flavors of Coulomb and Lennard-Jones interaction
 * treatment implemented for a large range of SIMD instruction sets for CPU
 * architectures as well as in CUDA and OpenCL for GPU architectures.
 * Additionally there is a reference CPU non-SIMD and a reference CPU
 * for GPU pair-list setup interaction kernel.
 *
 * The implementation of the kernels is based on the cluster non-bonded algorithm
 * which in the code is referred to as the NxN algorithms ("nbnxn_" prefix);
 * for details of the algorithm see DOI:10.1016/j.cpc.2013.06.003.
 *
 * Algorithmically, the non-bonded computation has two different modes:
 * A "classical" mode: generate a list every nstlist steps containing at least
 * all atom pairs up to a distance of rlistOuter and compute pair interactions
 * for all pairs that are within the interaction cut-off.
 * A "dynamic pruning" mode: generate an "outer-list" up to cut-off rlistOuter
 * every nstlist steps and prune the outer-list using a cut-off of rlistInner
 * every nstlistPrune steps to obtain a, smaller, "inner-list". This
 * results in fewer interaction computations and allows for a larger nstlist.
 * On a GPU, this dynamic pruning is performed in a rolling fashion, pruning
 * only a sub-part of the list each (second) step. This way it can often
 * overlap with integration and constraints on the CPU.
 * Currently a simple heuristic determines which mode will be used.
 *
 * TODO: add a summary list and brief descriptions of the different submodules:
 * search, CPU kernels, GPU glue code + kernels.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Anca Hamuraru <anca@streamcomputing.eu>
 * \author Teemu Virolainen <teemu@streamcomputing.eu>
 * \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *
 * TODO: add more authors!
 */

/*! \libinternal \file
 *
 * \brief This file contains the public interface of the non-bonded Verlet module
 * that implements the NxN cluster non-bonded algorithm to efficiently compute
 * pair forces.
 *
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 *
 * \inlibraryapi
 * \ingroup __module_nb_verlet
 */


#ifndef NB_VERLET_H
#define NB_VERLET_H

#include <memory>

#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

//! Help pass GPU-emulation parameters with type safety.
enum class EmulateGpuNonbonded : bool
{
    //! Do not emulate GPUs.
    No,
    //! Do emulate GPUs.
    Yes
};

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Nonbonded NxN kernel types: plain C, CPU SIMD, GPU, GPU emulation */
typedef enum
{
    nbnxnkNotSet = 0,
    nbnxnk4x4_PlainC,
    nbnxnk4xN_SIMD_4xN,
    nbnxnk4xN_SIMD_2xNN,
    nbnxnk8x8x8_GPU,
    nbnxnk8x8x8_PlainC,
    nbnxnkNR
} nbnxn_kernel_type;

/*! \brief Return a string identifying the kernel type.
 *
 * \param [in] kernel_type   nonbonded kernel types, takes values from the nbnxn_kernel_type enum
 * \returns                  a string identifying the kernel corresponding to the type passed as argument
 */
const char *lookup_nbnxn_kernel_name(int kernel_type);

/*! \brief Ewald exclusion types */
enum {
    ewaldexclTable, ewaldexclAnalytical
};

/*! \brief Atom locality indicator: local, non-local, all.
 *
 * Used for calls to:
 * gridding, pair-search, force calculation, x/f buffer operations
 * */
enum {
    eatLocal = 0, eatNonlocal = 1, eatAll
};

/*! \brief Tests for local atom range */
#define LOCAL_A(x)               ((x) == eatLocal)
/*! \brief Tests for non-local atom range */
#define NONLOCAL_A(x)            ((x) == eatNonlocal)
/*! \brief Tests for either local or non-local atom range */
#define LOCAL_OR_NONLOCAL_A(x)   (LOCAL_A(x) || NONLOCAL_A(x))

/*! \brief Interaction locality indicator
 *
 * Used in pair-list search/calculations in the following manner:
 *  - local interactions require local atom data and affect local output only;
 *  - non-local interactions require both local and non-local atom data and
 *    affect both local- and non-local output.
 */
enum {
    eintLocal = 0, eintNonlocal = 1
};

/*! \brief Tests for local interaction indicator */
#define LOCAL_I(x)               ((x) == eintLocal)
/*! \brief Tests for non-local interaction indicator */
#define NONLOCAL_I(x)            ((x) == eintNonlocal)

/*! \brief Flag to tell the nonbonded kernels whether to clear the force output buffers */
enum {
    enbvClearFNo, enbvClearFYes
};

/*! \libinternal
 *  \brief Non-bonded interaction group data structure. */
typedef struct nonbonded_verlet_group_t {
    nbnxn_pairlist_set_t  nbl_lists;   /**< pair list(s)                       */
    nbnxn_atomdata_t     *nbat;        /**< atom data                          */
    int                   kernel_type; /**< non-bonded kernel - see enum above */
    int                   ewald_excl;  /**< Ewald exclusion - see enum above   */
} nonbonded_verlet_group_t;

/*! \libinternal
 *  \brief Top-level non-bonded data structure for the Verlet-type cut-off scheme. */
typedef struct nonbonded_verlet_t {
    std::unique_ptr<NbnxnListParameters> listParams;      /**< Parameters for the search and list pruning setup */
    nbnxn_search_t                       nbs;             /**< n vs n atom pair searching data       */
    int                                  ngrp;            /**< number of interaction groups          */
    nonbonded_verlet_group_t             grp[2];          /**< local and non-local interaction group */

    gmx_bool                             bUseGPU;         /**< TRUE when non-bonded interactions are computed on a physical GPU */
    EmulateGpuNonbonded                  emulateGpu;      /**< true when non-bonded interactions are computed on the CPU using GPU-style pair lists */
    gmx_nbnxn_gpu_t                     *gpu_nbv;         /**< pointer to GPU nb verlet data     */
    int                                  min_ci_balanced; /**< pair list balancing parameter
                                                               used for the 8x8x8 GPU kernels    */
} nonbonded_verlet_t;

/*! \brief Getter for bUseGPU */
gmx_bool
usingGpu(nonbonded_verlet_t *nbv);

#ifdef __cplusplus
}
#endif

#endif /* NB_VERLET_H */
