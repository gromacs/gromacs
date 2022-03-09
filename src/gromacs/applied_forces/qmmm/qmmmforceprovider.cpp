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
/*! \internal \file
 * \brief
 * Implements force provider for QMMM
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "qmmmforceprovider.h"

#include <libcp2k.h>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \brief Helper function that dumps string into the file.
 *
 * \param[in] filename Name of the file to write.
 * \param[in] str String to write into the file.
 * \throws    std::bad_alloc if out of memory.
 * \throws    FileIOError on any I/O error.
 */
void writeStringToFile(const std::string& filename, const std::string& str)
{

    TextOutputFile fOut(filename);
    fOut.write(str.c_str());
    fOut.close();
}

} // namespace

QMMMForceProvider::QMMMForceProvider(const QMMMParameters& parameters,
                                     const LocalAtomSet&   localQMAtomSet,
                                     const LocalAtomSet&   localMMAtomSet,
                                     PbcType               pbcType,
                                     const MDLogger&       logger) :
    parameters_(parameters),
    qmAtoms_(localQMAtomSet),
    mmAtoms_(localMMAtomSet),
    pbcType_(pbcType),
    logger_(logger),
    box_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }
{
}

QMMMForceProvider::~QMMMForceProvider()
{
    if (force_env_ != -1)
    {
        cp2k_destroy_force_env(force_env_);
        if (GMX_LIB_MPI)
        {
            cp2k_finalize_without_mpi();
        }
        else
        {
            cp2k_finalize();
        }
    }
}

bool QMMMForceProvider::isQMAtom(index globalAtomIndex)
{
    return std::find(qmAtoms_.globalIndex().begin(), qmAtoms_.globalIndex().end(), globalAtomIndex)
           != qmAtoms_.globalIndex().end();
}

void QMMMForceProvider::appendLog(const std::string& msg)
{
    GMX_LOG(logger_.info).asParagraph().appendText(msg);
}

void QMMMForceProvider::initCP2KForceEnvironment(const t_commrec& cr)
{
    // Check that we have filename either defined in KVT or deduced from *.tpr name
    GMX_RELEASE_ASSERT(!parameters_.qmFileNameBase_.empty(),
                       "Names of CP2K input/output files is not defined in QMMMParameters");


    // Names of CP2K Input and Output files
    const std::string cp2kInputName  = parameters_.qmFileNameBase_ + ".inp";
    const std::string cp2kPdbName    = parameters_.qmFileNameBase_ + ".pdb";
    const std::string cp2kOutputName = parameters_.qmFileNameBase_ + ".out";

    // Write CP2K input if we are Master
    if (MASTER(&cr))
    {
        // In the CP2K Input we need to substitute placeholder with the actuall *.pdb file name
        writeStringToFile(cp2kInputName, formatString(parameters_.qmInput_.c_str(), cp2kPdbName.c_str()));

        // Write *.pdb with point charges for CP2K
        writeStringToFile(cp2kPdbName, parameters_.qmPdb_);
    }

    // Check if we have external MPI library
    if (GMX_LIB_MPI)
    {
        /* Attempt to init CP2K and create force environment in case we have an external MPI library.
         * _without_mpi - means that libcp2k will not call MPI_Init() itself,
         * but rather expects it to be called already
         */
        cp2k_init_without_mpi();

        // Here #if (GMX_LIB_MPI) needed because of MPI_Comm_c2f() usage
        // TODO: Probably there should be more elegant solution
#if GMX_LIB_MPI
        cp2k_create_force_env_comm(
                &force_env_, cp2kInputName.c_str(), cp2kOutputName.c_str(), MPI_Comm_c2f(cr.mpi_comm_mysim));
#endif
    }
    else
    {

        // If we have thread-MPI or no-MPI then we should initialize CP2P differently
        if (cr.nnodes > 1)
        {
            // CP2K could not use thread-MPI parallelization. In case -ntmpi > 1 throw an error.
            std::string msg =
                    "Using GROMACS with CP2K requires the GROMACS build to be configured "
                    "with an MPI library or to use a single thread-MPI rank (-ntmpi 1). "
                    "In the latter case, manual use of -ntomp is also advisable when "
                    "the node has many cores to fill them with threads.";
            GMX_THROW(NotImplementedError(msg.c_str()));
        }

        // Attempt to init CP2K and create force environment
        cp2k_init();
        cp2k_create_force_env(&force_env_, cp2kInputName.c_str(), cp2kOutputName.c_str());
    }

    // Set flag of successful initialization
    isCp2kLibraryInitialized_ = true;

} // namespace gmx

void QMMMForceProvider::calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput)
{
    // Total number of atoms in the system
    size_t numAtoms = qmAtoms_.numAtomsGlobal() + mmAtoms_.numAtomsGlobal();

    // Save box
    copy_mat(fInput.box_, box_);

    // Initialize PBC
    t_pbc pbc;
    set_pbc(&pbc, pbcType_, box_);

    /*
     * 1) If calculateForce called first time, then we need to init CP2K,
     *    as it was not possible during QMMMForceProvider constructor
     *    due to absence of full communication record at that point.
     */
    if (!isCp2kLibraryInitialized_)
    {
        try
        {
            initCP2KForceEnvironment(fInput.cr_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    /*
     * 2) We need to gather fInput.x_ in case of MPI / DD setup
     */

    // x - coordinates (gathered across nodes in case of DD)
    std::vector<RVec> x(numAtoms, RVec({ 0.0, 0.0, 0.0 }));

    // Fill cordinates of local QM atoms and add translation
    for (size_t i = 0; i < qmAtoms_.numAtomsLocal(); i++)
    {
        x[qmAtoms_.globalIndex()[qmAtoms_.collectiveIndex()[i]]] =
                fInput.x_[qmAtoms_.localIndex()[i]] + parameters_.qmTrans_;
    }

    // Fill cordinates of local MM atoms and add translation
    for (size_t i = 0; i < mmAtoms_.numAtomsLocal(); i++)
    {
        x[mmAtoms_.globalIndex()[mmAtoms_.collectiveIndex()[i]]] =
                fInput.x_[mmAtoms_.localIndex()[i]] + parameters_.qmTrans_;
    }

    // If we are in MPI / DD conditions then gather coordinates over nodes
    if (havePPDomainDecomposition(&fInput.cr_))
    {
        gmx_sum(3 * numAtoms, x.data()->as_vec(), &fInput.cr_);
    }

    // Put all atoms into the central box (they might be shifted out of it because of the translation)
    put_atoms_in_box(pbcType_, fInput.box_, ArrayRef<RVec>(x));

    /*
     * 3) Cast data to double format of libcp2k
     *    update coordinates and box in CP2K and perform QM calculation
     */
    // x_d - coordinates casted to linear dobule vector for CP2K with parameters_.qmTrans_ added
    std::vector<double> x_d(3 * numAtoms, 0.0);
    for (size_t i = 0; i < numAtoms; i++)
    {
        x_d[3 * i]     = static_cast<double>((x[i][XX]) / c_bohr2Nm);
        x_d[3 * i + 1] = static_cast<double>((x[i][YY]) / c_bohr2Nm);
        x_d[3 * i + 2] = static_cast<double>((x[i][ZZ]) / c_bohr2Nm);
    }

    // box_d - box_ casted to linear dobule vector for CP2K
    std::vector<double> box_d(9);
    for (size_t i = 0; i < DIM; i++)
    {
        box_d[3 * i]     = static_cast<double>(box_[0][i] / c_bohr2Nm);
        box_d[3 * i + 1] = static_cast<double>(box_[1][i] / c_bohr2Nm);
        box_d[3 * i + 2] = static_cast<double>(box_[2][i] / c_bohr2Nm);
    }

    // Update coordinates and box in CP2K
    cp2k_set_positions(force_env_, x_d.data(), 3 * numAtoms);
    cp2k_set_cell(force_env_, box_d.data());

    // Run CP2K calculation
    cp2k_calc_energy_force(force_env_);

    /*
     * 4) Get output data
     * We need to fill only local part into fOutput
     */

    // Only master process should add QM + QMMM energy
    if (MASTER(&fInput.cr_))
    {
        double qmEner = 0.0;
        cp2k_get_potential_energy(force_env_, &qmEner);
        fOutput->enerd_.term[F_EQM] += qmEner * c_hartree2Kj * c_avogadro;
    }

    // Get Forces they are in Hartree/Bohr and will be converted to kJ/mol/nm
    std::vector<double> cp2kForce(3 * numAtoms, 0.0);
    cp2k_get_forces(force_env_, cp2kForce.data(), 3 * numAtoms);

    // Fill forces on QM atoms first
    for (size_t i = 0; i < qmAtoms_.numAtomsLocal(); i++)
    {
        fOutput->forceWithVirial_.force_[qmAtoms_.localIndex()[i]][XX] +=
                static_cast<real>(cp2kForce[3 * qmAtoms_.globalIndex()[qmAtoms_.collectiveIndex()[i]]])
                * c_hartreeBohr2Md;

        fOutput->forceWithVirial_.force_[qmAtoms_.localIndex()[i]][YY] +=
                static_cast<real>(cp2kForce[3 * qmAtoms_.globalIndex()[qmAtoms_.collectiveIndex()[i]] + 1])
                * c_hartreeBohr2Md;

        fOutput->forceWithVirial_.force_[qmAtoms_.localIndex()[i]][ZZ] +=
                static_cast<real>(cp2kForce[3 * qmAtoms_.globalIndex()[qmAtoms_.collectiveIndex()[i]] + 2])
                * c_hartreeBohr2Md;
    }

    // Filll forces on MM atoms then
    for (size_t i = 0; i < mmAtoms_.numAtomsLocal(); i++)
    {
        fOutput->forceWithVirial_.force_[mmAtoms_.localIndex()[i]][XX] +=
                static_cast<real>(cp2kForce[3 * mmAtoms_.globalIndex()[mmAtoms_.collectiveIndex()[i]]])
                * c_hartreeBohr2Md;

        fOutput->forceWithVirial_.force_[mmAtoms_.localIndex()[i]][YY] +=
                static_cast<real>(cp2kForce[3 * mmAtoms_.globalIndex()[mmAtoms_.collectiveIndex()[i]] + 1])
                * c_hartreeBohr2Md;

        fOutput->forceWithVirial_.force_[mmAtoms_.localIndex()[i]][ZZ] +=
                static_cast<real>(cp2kForce[3 * mmAtoms_.globalIndex()[mmAtoms_.collectiveIndex()[i]] + 2])
                * c_hartreeBohr2Md;
    }
};

} // namespace gmx
