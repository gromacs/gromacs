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
public:
    //! Total mass in state A
    real tmassA;
    //! Total mass in state B
    real tmassB;
    //! Total mass
    real tmass;
    //! Do we have multiple center of mass motion removal groups
    bool bVCMgrps;
    //! Do we have any virtual sites?
    bool haveVsites;
    //! Do we have atoms that are frozen along 1 or 2 (not 3) dimensions?
    bool havePartiallyFrozenAtoms;
    //! Number of perturbed atoms
    int nPerturbed;
    //! Number of atoms for which the mass is perturbed
    int nMassPerturbed;
    //! Number of atoms for which the charge is perturbed
    int nChargePerturbed;
    //! Number of atoms for which the type is perturbed
    int nTypePerturbed;
    //! Do we have orientation restraints
    bool bOrires;
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
    //! Group index for orientation restraints
    std::vector<unsigned short> cORF;
    //! The mass lambda value used to update the contents of the struct
    real massLambda;

    // TODO make this private
    MDAtoms(bool rankHasPmeGpuTask, bool useGpuForUpdate, const DeviceStreamManager* deviceStreamManager);
    //! Builder function.
    friend std::unique_ptr<MDAtoms> makeMDAtoms(FILE*                      fp,
                                                const gmx_mtop_t&          mtop,
                                                const t_inputrec&          ir,
                                                bool                       rankHasPmeGpuTask,
                                                const DeviceStreamManager* deviceStreamManager);
};

//! Builder function for MdAtomsWrapper.
std::unique_ptr<MDAtoms> makeMDAtoms(FILE*                      fp,
                                     const gmx_mtop_t&          mtop,
                                     const t_inputrec&          ir,
                                     bool                       useGpuForPme,
                                     bool                       useGpuForUpdate,
                                     const DeviceStreamManager* deviceStreamManager);

} // namespace gmx

/*! \brief This routine copies the atoms->atom struct into mdAtoms.
 *
 * \param[in]    mtop          The molecular topology.
 * \param[in]    inputrec      The input record.
 * \param[in]    nindex        If nindex>=0 we are doing DD.
 * \param[in]    index         Lookup table for global atom index.
 * \param[in]    numHomeAtoms  Number of home atoms on this rank.
 * \param[inout] mdAtoms       Data set up by this routine.
 *
 * If index!=NULL only the indexed atoms are copied.
 * For the masses the A-state (lambda=0) mass is used.
 * Sets mdAtoms->massLambda = 0.
 * In free-energy runs, update_mdatoms() should be called after atoms2md()
 * to set the masses corresponding to the value of the mass lambda at each step.
 */
void atoms2md(const gmx_mtop_t&        mtop,
              const t_inputrec&        inputrec,
              int                      nindex,
              gmx::ArrayRef<const int> index,
              int                      numHomeAtoms,
              gmx::MDAtoms*            mdAtoms);

void update_mdatoms(gmx::MDAtoms* mdAtoms, real massLambda);
/* When necessary, sets all the mass parameters to values corresponding
 * to the (mass) free-energy parameter lambda.
 * Sets mdAtoms->massLambda = massLambda.
 */

#endif
