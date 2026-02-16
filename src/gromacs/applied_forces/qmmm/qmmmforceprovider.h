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
 * Declares force provider for QMMM
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_QMMMFORCEPROVIDER_H

#include <string>
#include <string_view>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/logger.h"

#include "qmmmforceproviderstate.h"
#include "qmmmtypes.h"

namespace gmx
{

CLANG_DIAGNOSTIC_IGNORE("-Wunused-private-field")

struct MDModulesWriteCheckpointData;
class MpiComm;

//! Type for CP2K force environment handle
typedef int force_env_t;

/*! \internal \brief
 * Implements IForceProvider for QM/MM.
 */
class QMMMForceProvider final : public IForceProvider
{
public:
    QMMMForceProvider(const QMMMParameters&         parameters,
                      const LocalAtomSet&           localQMAtomSet,
                      const LocalAtomSet&           localMMAtomSet,
                      PbcType                       pbcType,
                      const MDLogger&               logger,
                      const MpiComm&                mpiComm,
                      const QMMMForceProviderState& state);

    //! Destruct force provider for QMMM and finalize libcp2k
    ~QMMMForceProvider();

    /*!\brief Calculate forces of QMMM.
     * \param[in] fInput input for force provider
     * \param[out] fOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& fInput, ForceProviderOutput* fOutput) override;

    /*! \brief Write internal QMMM data to checkpoint file.
     * \param[in] checkpointWriting enables writing to the Key-Value-Tree
     *                              that is used for storing the checkpoint
     *                              information
     * \param[in] moduleName names the module that is checkpointing this force-provider
     *
     */
    void writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting, std::string_view moduleName);

private:
    //! Write message to the log
    void appendLog(const std::string& msg);

    /*!\brief Check if atom belongs to the global index of qmAtoms_
     * \param[in] globalAtomIndex global index of the atom to check
     */
    bool isQMAtom(Index globalAtomIndex);

    /*!\brief Initialization of QM program.
     * \param[in] mpiComm communication object
     */
    void initCP2KForceEnvironment(const MpiComm& mpiComm);

    const QMMMParameters& parameters_;
    const LocalAtomSet&   qmAtoms_;
    const LocalAtomSet&   mmAtoms_;
    const PbcType         pbcType_;
    const MDLogger&       logger_;

    //! Struct holding the information stored in the checkpoint file
    QMMMForceProviderState state_;

    //! Internal copy of PBC box
    matrix box_;

    //! CP2K force environment handle
    force_env_t force_env_ = -1;
};

CLANG_DIAGNOSTIC_RESET

//! Returns information for describing the CP2K QM/MM support
std::string qmmmDescription();

} // namespace gmx

#endif // GMX_APPLIED_FORCES_QMMMFORCEPROVIDER_H
