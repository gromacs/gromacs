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
#ifndef GMX_MDTYPES_TYPES_FORCEREC_H
#define GMX_MDTYPES_TYPES_FORCEREC_H

#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "locality.h"

/* Abstract type for PME that is defined only in the routine that use them. */
struct gmx_pme_t;
struct bonded_threading_t;
class DispersionCorrection;
class ListedForces;
class CpuPpLongRangeNonbondeds;
struct t_fcdata;
struct t_forcetable;
struct interaction_const_t;

namespace gmx
{
struct nonbonded_verlet_t;
class DeviceStreamManager;
class ListedForcesGpu;
class GpuForceReduction;
class ForceProviders;
class MdGpuGraph;
class StatePropagatorDataGpu;
class PmePpCommGpu;
class WholeMoleculeTransform;
} // namespace gmx

/* Value to be used in mdrun for an infinite cut-off.
 * Since we need to compare with the cut-off squared,
 * this value should be slighlty smaller than sqrt(GMX_FLOAT_MAX).
 */
#define GMX_CUTOFF_INF 1E+18
//! Check the cuttoff
real cutoff_inf(real cutoff);

/* Forward declaration of type for managing Ewald tables */
struct gmx_ewald_tab_t;

/*! \brief Helper force buffers for ForceOutputs
 *
 * This class stores intermediate force buffers that are used
 * internally in the force calculation and which are reduced into
 * the output force buffer passed to the force calculation.
 */
class ForceHelperBuffers
{
public:
    /*! \brief Constructs helper buffers
     *
     * When the forces that will be accumulated with help of these buffers
     * have direct virial contributions, set the parameter to true, so
     * an extra force buffer is available for these forces to enable
     * correct virial computation.
     */
    ForceHelperBuffers(bool haveDirectVirialContributions);

    //! Returns whether we have a direct virial contribution force buffer
    bool haveDirectVirialContributions() const { return haveDirectVirialContributions_; }

    //! Returns the buffer for direct virial contributions
    gmx::ArrayRef<gmx::RVec> forceBufferForDirectVirialContributions()
    {
        GMX_ASSERT(haveDirectVirialContributions_, "Buffer can only be requested when present");
        return forceBufferForDirectVirialContributions_;
    }

    //! Returns the buffer for shift forces, size c_numShiftVectors
    gmx::ArrayRef<gmx::RVec> shiftForces() { return shiftForces_; }

    //! Resizes the direct virial contribution buffer, when present
    void resize(int numAtoms);

private:
    //! True when we have contributions that are directly added to the virial
    bool haveDirectVirialContributions_ = false;
    //! Force buffer for force computation with direct virial contributions
    std::vector<gmx::RVec> forceBufferForDirectVirialContributions_;
    //! Shift force array for computing the virial, size c_numShiftVectors
    std::vector<gmx::RVec> shiftForces_;
};
// NOLINTNEXTLINE (clang-analyzer-optin.performance.Padding)
struct t_forcerec
{
    // Declare an explicit constructor and destructor, so they can be
    // implemented in a single source file, so that not every source
    // file that includes this one needs to understand how to find the
    // destructors of the objects pointed to by unique_ptr members.
    t_forcerec();
    ~t_forcerec();

    std::unique_ptr<interaction_const_t> ic;

    /* PBC stuff */
    PbcType pbcType = PbcType::Xyz;
    //! Tells whether atoms inside a molecule can be in different periodic images,
    //  i.e. whether we need to take into account PBC when computing distances inside molecules.
    //  This determines whether PBC must be considered for e.g. bonded interactions.
    bool            bMolPBC     = false;
    RefCoordScaling rc_scaling  = RefCoordScaling::No;
    gmx::RVec       posres_com  = { 0, 0, 0 };
    gmx::RVec       posres_comB = { 0, 0, 0 };

    // Tells whether the box is continuosly deformed
    bool haveBoxDeformation = false;

    bool use_simd_kernels = false;

    /* Interaction for calculated in kernels. In many cases this is similar to
     * the electrostatics settings in the inputrecord, but the difference is that
     * these variables always specify the actual interaction in the kernel - if
     * we are tabulating reaction-field the inputrec will say reaction-field, but
     * the kernel interaction will say cubic-spline-table. To be safe we also
     * have a kernel-specific setting for the modifiers - if the interaction is
     * tabulated we already included the inputrec modification there, so the kernel
     * modification setting will say 'none' in that case.
     */
    NbkernelElecType     nbkernel_elec_interaction = NbkernelElecType::None;
    NbkernelVdwType      nbkernel_vdw_interaction  = NbkernelVdwType::None;
    InteractionModifiers nbkernel_elec_modifier    = InteractionModifiers::None;
    InteractionModifiers nbkernel_vdw_modifier     = InteractionModifiers::None;

    /* Cut-Off stuff.
     * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
     */
    real rlist = 0;

    /* Charge sum for topology A/B ([0]/[1]) for Ewald corrections */
    std::array<double, 2> qsum  = { 0 };
    std::array<double, 2> q2sum = { 0 };
    std::array<double, 2> c6sum = { 0 };

    /* Dispersion correction stuff */
    std::unique_ptr<DispersionCorrection> dispersionCorrection;

    /* Fudge factors */
    real fudgeQQ = 0;

    std::unique_ptr<t_forcetable> pairsTable; /* for 1-4 interactions, [pairs] and [pairs_nb] */

    /* Free energy */
    FreeEnergyPerturbationType efep = FreeEnergyPerturbationType::No;

    /* Information about atom properties for the molecule blocks in the global topology */
    std::vector<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock;
    /* Information about atom properties for local and non-local atoms */
    std::vector<int32_t> atomInfo;

    std::vector<gmx::RVec> shift_vec;

    std::unique_ptr<gmx::WholeMoleculeTransform> wholeMoleculeTransform;

    /* The Nbnxm Verlet non-bonded machinery */
    std::unique_ptr<gmx::nonbonded_verlet_t> nbv;

    /* The wall tables (if used) */
    int                                                     nwall = 0;
    std::vector<std::vector<std::unique_ptr<t_forcetable>>> wall_tab;

    /* The number of atoms participating in do_force_lowlevel */
    int natoms_force = 0;
    /* The number of atoms participating in force calculation and constraints */
    int natoms_force_constr = 0;

    /* List of helper buffers for ForceOutputs, one for each time step with MTS */
    std::vector<ForceHelperBuffers> forceHelperBuffers;

    /* Data for PPPM/PME/Ewald */
    gmx_pme_t*   pmedata                = nullptr;
    LongRangeVdW ljpme_combination_rule = LongRangeVdW::Geom;

    /* Non bonded Parameter lists */
    int               ntype          = 0; /* Number of atom types */
    bool              haveBuckingham = false;
    std::vector<real> nbfp;
    std::vector<real> ljpme_c6grid; /* C6-values used on grid in LJPME */

    /* Energy group pair flags */
    int* egp_flags = nullptr;

    /* Shell molecular dynamics flexible constraints */
    real fc_stepsize = 0;

    /* If > 0 signals Test Particle Insertion,
     * the value is the number of atoms of the molecule to insert
     * Only the energy difference due to the addition of the last molecule
     * should be calculated.
     */
    int n_tpi = 0;

    /* Limit for printing large forces, negative is don't print */
    real print_force = 0;

    /* User determined parameters, copied from the inputrec */
    int  userint1  = 0;
    int  userint2  = 0;
    int  userint3  = 0;
    int  userint4  = 0;
    real userreal1 = 0;
    real userreal2 = 0;
    real userreal3 = 0;
    real userreal4 = 0;

    /* Data for special listed force calculations */
    std::unique_ptr<t_fcdata> fcdata;

    // The listed forces calculation data, 1 entry or multiple entries with multiple time stepping
    std::vector<ListedForces> listedForces;

    // The listed forces calculation data for GPU
    std::unique_ptr<gmx::ListedForcesGpu> listedForcesGpu;

    // The long range non-bonded forces
    std::unique_ptr<CpuPpLongRangeNonbondeds> longRangeNonbondeds;

    gmx::ForceProviders* forceProviders = nullptr;

    // The stateGpu object is created in runner, forcerec just keeps the copy of the pointer.
    // TODO: This is not supposed to be here. StatePropagatorDataGpu should be a part of
    //       general StatePropagatorData object that is passed around
    gmx::StatePropagatorDataGpu* stateGpu = nullptr;
    // TODO: Should not be here. This is here only to pass the pointer around.
    gmx::DeviceStreamManager* deviceStreamManager = nullptr;

    /* For PME-PP GPU communication */
    std::unique_ptr<gmx::PmePpCommGpu> pmePpCommGpu;

    /* For GPU force reduction (on both local and non-local atoms) */
    gmx::EnumerationArray<gmx::AtomLocality, std::unique_ptr<gmx::GpuForceReduction>> gpuForceReduction;

    gmx::EnumerationArray<MdGraphEvenOrOddStep, std::unique_ptr<gmx::MdGpuGraph>> mdGraph;
};

/* Important: Starting with Gromacs-4.6, the values of c6 and c12 in the nbfp array have
 * been scaled by 6.0 or 12.0 to save flops in the kernels. We have corrected this everywhere
 * in the code, but beware if you are using these macros externally.
 */
#define C6(nbfp, ntp, ai, aj) (nbfp)[2 * ((ntp) * (ai) + (aj))]
#define C12(nbfp, ntp, ai, aj) (nbfp)[2 * ((ntp) * (ai) + (aj)) + 1]
#define BHAMC(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj))]
#define BHAMA(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj)) + 1]
#define BHAMB(nbfp, ntp, ai, aj) (nbfp)[3 * ((ntp) * (ai) + (aj)) + 2]

#endif
