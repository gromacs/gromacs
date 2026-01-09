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
 * Declares QMMM force provider state for checkpointing
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMFORCEPROVIDERSTATE_H
#define GMX_APPLIED_FORCES_QMMMFORCEPROVIDERSTATE_H

#include <string_view>

#include "gromacs/utility/vectypes.h"

namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;

/*! \internal
 * \brief Class defining the internal QMMM force provider state.
 */
class QMMMForceProviderState
{

public:
    //! Default constructor
    QMMMForceProviderState() = default;

    /*! \brief Write internal QMMM data into a key value tree.
     * The entries to the kvt are identified with identifier, so that a variable
     * is indentified with the key "identifier-variablename"
     *
     * \param[in] kvtBuilder enables writing to the Key-Value-Tree
     *                              the state is written to
     *
     * \param[in] identifier denotes the module that is checkpointing the data
     */
    void writeState(KeyValueTreeObjectBuilder kvtBuilder, std::string_view identifier) const;

    /*! \brief Read the internal parameters from the checkpoint file on master
     * \param[in] kvtData holding the checkpoint information
     * \param[in] identifier identifies the data in a key-value-tree
     */
    void readState(const KeyValueTreeObject& kvtData, std::string_view identifier);

    /*! \brief Returns whether a QMMM state was read from checkpoint.
     * \return true if state was read, false otherwise
     */
    bool isStateRead() const { return stateRead_; }

    /*! \brief Returns QMMM atoms translation vector.
     * \return stored QM translation vector
     */
    const RVec& qmTrans() const { return qmTrans_; }

    /*! \brief Sets QM translation vector.
     * \param[in] qmTrans QM translation vector to set
     */
    void setQMTrans(const RVec& qmTrans) { qmTrans_ = qmTrans; }

private:
    /*! \brief Indicate if a QMMM state was read.
     */
    bool stateRead_ = false;

    /*! \brief Information about QMMM atoms translation vector during simulation.
     *
     * This is saved to the checkpoint. */
    RVec qmTrans_ = { 0.0, 0.0, 0.0 };
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_QMMMFORCEPROVIDERSTATE_H
