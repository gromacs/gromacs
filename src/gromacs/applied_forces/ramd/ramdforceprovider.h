/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * \brief
 * Declares force provider for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H

#include <string>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

#include "ramdoutputprovider.h"
#include "ramdparameters.h"
#include "randomsphericaldirectiongenerator.h"

namespace gmx
{

/*! \internal \brief
 * Declares IForceProvider for RAMD.
 */
class RAMDForceProvider final : public IForceProvider
{
public:
    RAMDForceProvider(const RAMDParameters&                             parameters,
                      const std::vector<std::unique_ptr<LocalAtomSet>>& localAtoms,
                      const gmx_mtop_t&                                 topology,
                      PbcType                                           pbcType,
                      const MDLogger&                                   logger,
                      RAMDOutputProvider&                               ramdOutputProvider);

    //! Destruct force provider for RAMD
    ~RAMDForceProvider();

    /*!\brief Calculate forces of RAMD.
     * \param[in] fInput input for force provider
     * \param[out] fOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput) override;

private:
    DVec calc_com(ArrayRef<const RVec> x, const std::vector<Index>& indices)
    {
        DVec com        = DVec(0.0, 0.0, 0.0);
        real total_mass = 0.0;
        for (auto idx : indices)
        {
            const real mass = mTopLookUp_.getAtomParameters(idx).m;
            for (int j = 0; j < DIM; ++j)
            {
                com[j] += mass * x[idx][j];
            }
            total_mass += mass;
        }
        if (total_mass > 0.0)
        {
            com /= total_mass;
        }
        return com;
    }

    //! The parameters for RAMD
    const RAMDParameters& parameters_;

    //! Reference to local atom sets
    const std::vector<std::unique_ptr<LocalAtomSet>>& localAtoms_;

    //! Reference to topology
    const gmx_mtop_t& topology_;

    const PbcType       pbcType_;
    const MDLogger&     logger_;
    RAMDOutputProvider& ramdOutputProvider_;

    //! Random pull direction
    RandomSphericalDirectionGenerator random_spherical_direction_generator;

    //! Current pull direction
    std::vector<DVec> direction_;

    //! COM of receptor of last RAMD evaluation step
    std::vector<DVec> com_rec_prev_;

    //! COM of ligand of last RAMD evaluation step
    std::vector<DVec> com_lig_prev_;

    //! Has the ligand left his binding site?
    std::vector<int> ligand_exited_;

    //! Control trajectory output
    gmx_bool write_trajectory_;

    //! Lookup for molecule topology information
    MTopLookUp mTopLookUp_;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_RAMDFORCEPROVIDER_H
