/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Implementation of the Plumed force provider class
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */

#include "plumedforceprovider.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/units.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/utility/exceptions.h"

#include "plumedOptions.h"

#define __PLUMED_WRAPPER_FORTRAN 0 // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_LINK_RUNTIME 1 // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_EXTERN 0       // NOLINT(bugprone-reserved-identifier)

#define __PLUMED_WRAPPER_CXX 1            // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX11 1       // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_LIBCXX17 1       // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_WRAPPER_IMPLEMENTATION 1 // NOLINT(bugprone-reserved-identifier)
#define __PLUMED_HAS_DLOPEN               // NOLINT(bugprone-reserved-identifier)

#include "external/plumed/Plumed.h"

namespace gmx
{

PlumedForceProvider::~PlumedForceProvider() = default;
PlumedForceProvider::PlumedForceProvider(const PlumedOptions& options)
try : plumed_(std::make_unique<PLMD::Plumed>())
{
    // I prefer to pass a struct with data because it stops the coupling
    // at the implementation and not at the function signature:
    // less code to edit when adding new options :)
#if GMX_THREAD_MPI
    if (options.cr_->nnodes > 1)
    {
        GMX_THROW(InvalidInputError(
                "plumed MPI interface is not compatible with THREAD_MPI when uses more than one "
                "rank"));
    }
#endif
    int real_precision = sizeof(real);
    plumed_->cmd("setRealPrecision", &real_precision);
    real energyUnits = 1.0;
    plumed_->cmd("setMDEnergyUnits", &energyUnits);
    real lengthUnits = 1.0;
    plumed_->cmd("setMDLengthUnits", &lengthUnits);
    real timeUnits = 1.0;
    plumed_->cmd("setMDTimeUnits", &timeUnits);

    plumed_->cmd("setPlumedDat", options.plumedFile_.c_str());

    plumed_->cmd("getApiVersion", &plumedAPIversion_);
    if (plumedAPIversion_ > 1)
    {
        /* setting kbT is only implemented with api>1) */
        if (options.ensembleTemperature_.has_value())
        {
            real kbT = options.ensembleTemperature_.value() * gmx::c_boltz;
            plumed_->cmd("setKbT", &kbT);
        }
    }

    if (plumedAPIversion_ > 2)
    {
        if ((options.startingBehavior_ != StartingBehavior::NewSimulation))
        {
            int res = 1;
            plumed_->cmd("setRestart", &res);
        }
    }

    if (havePPDomainDecomposition(options.cr_))
    {
        plumed_->cmd("setMPIComm", &options.cr_->mpi_comm_mygroup);
    }

    plumed_->cmd("setNatoms", options.natoms_);
    plumed_->cmd("setMDEngine", "gromacs");
    plumed_->cmd("setLogFile", "PLUMED.OUT");

    plumed_->cmd("setTimestep", &options.simulationTimeStep_);
    plumed_->cmd("init", nullptr);
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while initializing the PLUMED force provider:\n") + ex.what()));
}

void PlumedForceProvider::writeCheckpointData()
try
{
    if (plumedAPIversion_ > 3)
    {
        int checkp = 1;
        plumed_->cmd("doCheckPoint", &checkp);
    }
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was writing the checkpoint data\n:") + ex.what()));
}

void PlumedForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                          ForceProviderOutput*      forceProviderOutput)
try
{
    // setup: these instructions in the original patch are BEFORE do_force()
    // now this is called within do_force(), but this does not impact the results
    const t_commrec* cr    = &(forceProviderInput.cr_);
    long int         lstep = forceProviderInput.step_;
    plumed_->cmd("setStepLong", &lstep);
    if (haveDDAtomOrdering(*cr))
    {
        int nat_home = dd_numHomeAtoms(*cr->dd);
        plumed_->cmd("setAtomsNlocal", &nat_home);
        plumed_->cmd("setAtomsGatindex", cr->dd->globalAtomIndices.data());
    }

    plumed_->cmd("setPositions", &(forceProviderInput.x_.data()->as_vec()[0]));
    plumed_->cmd("setMasses", forceProviderInput.massT_.data());
    plumed_->cmd("setCharges", forceProviderInput.chargeA_.data());
    plumed_->cmd("setBox", &forceProviderInput.box_[0][0]);

    plumed_->cmd("prepareCalc", nullptr);

    int plumedWantsToStop = 0;
    plumed_->cmd("setStopFlag", &plumedWantsToStop);

    real* fOut = &(forceProviderOutput->forceWithVirial_.force_.data()->as_vec()[0]);
    plumed_->cmd("setForces", fOut);

    matrix plumed_vir;
    clear_mat(plumed_vir);
    plumed_->cmd("setVirial", &plumed_vir[0][0]);

    // end setup: these instructions in the original patch are BEFORE do_force()
    // in the original patch do_force() was called HERE

    // Do the work
    plumed_->cmd("performCalc", nullptr);

    msmul(plumed_vir, 0.5, plumed_vir);
    forceProviderOutput->forceWithVirial_.addVirialContribution(plumed_vir);
}
catch (const std::exception& ex)
{
    GMX_THROW(InternalError(
            std::string("An error occurred while PLUMED was calculating the forces\n:") + ex.what()));
}
} // namespace gmx
