/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements nblib tpr reading
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "nblib/tpr.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#include "nblib/box.h"

namespace nblib
{

class TprReader::Impl
{
public:
    gmx_mtop_t molecularTopology_;
    t_forcerec forceRecord_;
    t_inputrec inputRecord_;
    t_state    globalState_;

    std::unique_ptr<gmx::MDAtoms>   mdAtoms_;
    std::unique_ptr<gmx_localtop_t> localtopology_;
};


TprReader::TprReader(std::string filename) : impl_(new Impl)
{
    gmx::SimulationWorkload simulationWorkload;

    // If the file does not exist, this function will throw
    PartialDeserializedTprFile partialDeserializedTpr = read_tpx_state(
            filename.c_str(), &impl_->inputRecord_, &impl_->globalState_, &impl_->molecularTopology_);

    // Currently dispersion correction not supported by nblib
    impl_->inputRecord_.eDispCorr = DispersionCorrectionType::No;

    // init forcerec
    t_commrec           commrec{};
    gmx::ForceProviders forceProviders;
    impl_->forceRecord_.forceProviders = &forceProviders;
    init_forcerec(nullptr,
                  gmx::MDLogger(),
                  simulationWorkload,
                  &impl_->forceRecord_,
                  impl_->inputRecord_,
                  impl_->molecularTopology_,
                  &commrec,
                  impl_->globalState_.box,
                  nullptr,
                  nullptr,
                  gmx::ArrayRef<const std::string>{},
                  -1);

    nonbondedParameters_ = makeNonBondedParameterLists(impl_->molecularTopology_.ffparams.atnr,
                                                       impl_->molecularTopology_.ffparams.iparams,
                                                       impl_->forceRecord_.haveBuckingham);

    impl_->localtopology_ = std::make_unique<gmx_localtop_t>(impl_->molecularTopology_.ffparams);
    gmx_mtop_generate_local_top(impl_->molecularTopology_, impl_->localtopology_.get(), false);
    exclusionListElements_ = std::vector<int>(impl_->localtopology_->excls.elementsView().begin(),
                                              impl_->localtopology_->excls.elementsView().end());
    exclusionListRanges_   = std::vector<int>(impl_->localtopology_->excls.listRangesView().begin(),
                                            impl_->localtopology_->excls.listRangesView().end());

    int ntopatoms = impl_->molecularTopology_.natoms;
    impl_->mdAtoms_ = gmx::makeMDAtoms(nullptr, impl_->molecularTopology_, impl_->inputRecord_, false);
    atoms2md(impl_->molecularTopology_, impl_->inputRecord_, -1, {}, ntopatoms, impl_->mdAtoms_.get());
    update_mdatoms(impl_->mdAtoms_->mdatoms(), impl_->inputRecord_.fepvals->init_lambda);
    auto numParticles = impl_->mdAtoms_->mdatoms()->nr;
    charges_.resize(numParticles);
    particleTypeIdOfAllParticles_.resize(numParticles);
    inverseMasses_.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
    {
        charges_[i]                      = impl_->mdAtoms_->mdatoms()->chargeA[i];
        particleTypeIdOfAllParticles_[i] = impl_->mdAtoms_->mdatoms()->typeA[i];
        inverseMasses_[i]                = impl_->mdAtoms_->mdatoms()->invmass[i];
    }
    particleInteractionFlags_ = impl_->forceRecord_.atomInfo;

    for (int i = 0; i < dimSize; ++i)
    {
        for (int j = 0; j < dimSize; ++j)
        {
            boxMatrix_[j + i * dimSize] = impl_->globalState_.box[i][j];
        }
    }

    coordinates_.assign(impl_->globalState_.x.begin(), impl_->globalState_.x.end());
    velocities_.assign(impl_->globalState_.v.begin(), impl_->globalState_.v.end());

    // Copy listed interactions data
    interactionDefinitions_ = std::make_unique<InteractionDefinitions>(impl_->localtopology_->idef);
    ffparams_               = std::make_unique<gmx_ffparams_t>(impl_->molecularTopology_.ffparams);
}

Box TprReader::getBox() const
{
    return Box(boxMatrix_);
}

TprReader::~TprReader() {}

} // namespace nblib
