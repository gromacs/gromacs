/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#ifndef GROMACS_WORKFLOW_IMPL_H
#define GROMACS_WORKFLOW_IMPL_H

/*! \internal \file
 * \brief Implementation details for Workflow infrastructure.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <memory>
#include <string>

#include "gmxapi/exceptions.h"

// Local module internal headers.
#include "workflow.h"

namespace gmxapi
{

class WorkflowKeyError : public BasicException<WorkflowKeyError>
{
public:
    using BasicException::BasicException;
};

/*!
 * \brief Work graph node for MD simulation.
 */
class MDNodeSpecification : public NodeSpecification
{
public:
    //! Uses parameter type of base class.
    using NodeSpecification::paramsType;

    /*!
     * \brief Simulation node from file input
     *
     * \param filename TPR input filename.
     */
    explicit MDNodeSpecification(const std::string& filename);

    /*
     * \brief Implement NodeSpecification::clone()
     *
     * \returns a node to launch a simulation from the same input as this
     *
     * Returns nullptr if clone is not possible.
     */
    std::unique_ptr<NodeSpecification> clone() override;

    /*! \brief Implement NodeSpecification::params()
     *
     * \return Copy of internal params value.
     */
    paramsType params() const noexcept override;

private:
    //! The TPR input filename, set during construction
    paramsType tprfilename_;
};


} // end namespace gmxapi

#endif // GROMACS_WORKFLOW_IMPL_H
