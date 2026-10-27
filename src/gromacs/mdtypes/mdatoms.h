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
#ifndef GMX_MDTYPES_MDATOMS_H
#define GMX_MDTYPES_MDATOMS_H

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/booltype.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_inputrec;

enum class ParticleType : int;

namespace gmx
{
template<typename T>
class ArrayRef;
class DeviceStreamManager;

/*! \libinternal
 * \brief Contains a C-style t_mdatoms while managing some of its
 * memory with C++ vectors with allocators.
 *
 * \todo The group-scheme kernels needed a plain C-style t_mdatoms, so
 * this type combines that with the memory management needed for
 * efficient PME on GPU transfers. The mdAtoms_ member should be
 * removed. */
class MDAtoms
{
private:
    //! Molecular topology
    const gmx_mtop_t& mtop_;
    //! Input record
    const t_inputrec& inputrec_;

public:
    //! Total mass in state A
    real tmassA;
    //! Total mass in state B
    real tmassB;
    //! Total mass
    real tmass;
    // Are there center-of-mass-motion removal groups in the whole system?
    bool bVCMgrps = false;
    //! Number of perturbed atoms
    int nPerturbed = 0;
    //! Number of atoms for which the mass is perturbed
    int nMassPerturbed = 0;
    //! Number of atoms for which the charge is perturbed
    int nChargePerturbed = 0;
    //! Number of atoms for which the type is perturbed
    int nTypePerturbed = 0;
    //! Number of home atoms in this domain
    int numHomeAtoms;
    //! Atomic mass in A state
    std::vector<real> massA;
    //! Atomic mass in B state
    std::vector<real> massB;
    //! Atomic mass in present state
    std::vector<real> massT;
    //! Inverse atomic mass per atom, 0 for vsites and shells
    PaddedHostVector<real> invmass;
    //! Inverse atomic mass per atom and dimension, 0 for vsites, shells and frozen dimensions
    std::vector<RVec> invMassPerDim;
    //! Atomic charges in the A state
    PaddedHostVector<real> chargeA;
    //! Atomic charges in the B state
    PaddedHostVector<real> chargeB;
    //! Dispersion constant C6 in A state
    std::vector<real> sqrt_c6A;
    //! Dispersion constant C6 in A state
    std::vector<real> sqrt_c6B;
    //! Van der Waals radius sigma in the A state
    std::vector<real> sigmaA;
    //! Van der Waals radius sigma in the B state
    std::vector<real> sigmaB;
    //! Van der Waals radius sigma^3 in the A state
    std::vector<real> sigma3A;
    //! Van der Waals radius sigma^3 in the B state
    std::vector<real> sigma3B;
    //! Is this atom perturbed?
    std::vector<BoolType> bPerturbed;
    //! Type of atom in the A state
    std::vector<int> typeA;
    //! Type of atom in the B state
    std::vector<int> typeB;
    //! Particle type
    std::vector<ParticleType> ptype;
    //! Group index for temperature coupling
    HostVector<unsigned short> cTC;
    //! Group index for energy matrix
    std::vector<unsigned short> cENER;
    //! Group index for acceleration
    std::vector<unsigned short> cACC;
    //! Group index for freezing
    std::vector<unsigned short> cFREEZE;
    //! Group index for center of mass motion removal
    std::vector<unsigned short> cVCM;
    //! Group index for user 1
    std::vector<unsigned short> cU1;
    //! Group index for user 2
    std::vector<unsigned short> cU2;

    MDAtoms(const gmx_mtop_t&          mtop,
            const t_inputrec&          inputrec,
            bool                       rankHasPmeGpuTask,
            bool                       useGpuForUpdate,
            const DeviceStreamManager* deviceStreamManager);

    /*! \brief Write out per-atom arrays for this domain
     *
     * Force, update, and constraint modules need various per-atom
     * fields that are efficient to look up when fully written out
     * for the whole domain. As particles may have moved between
     * domain repartitionings, these vectors must be re-filled at
     * the start of every domain lifetime.
     *
     * \param[in]    numTotalAtoms  The number of atoms used on this domain
     * \param[in]    index          Lookup table for global atom index.
     * \param[in]    numHomeAtoms   Number of home atoms for this domain
     * \param[in]    massLambda     The FEP lambda value relevant to mass perturbation
     *
     * When using DD, numTotalAtoms indicates the size of the map from
     * local to global atom indices.
     *
     * If \c index is non-empty only the indexed atoms are copied.
     */
    void updateAfterPartition(int numTotalAtoms, ArrayRef<const int> index, int numHomeAtoms, real massLambda);

    /*! \brief This routine updates the mass-related values in \c mdAtoms to
     * respect \c massLambda
     *
     * Intended to cater for possible per-step changes to the
     * lambda value for masses with slow-growth FEP, but also
     * called by \c atoms2md() during reinit after DD repartitioning.
     *
     * \param[in]    massLambda    The FEP lambda value relevant to mass perturbation
     */
    void updateForMassLambda(real massLambda);
};

} // namespace gmx

#endif
