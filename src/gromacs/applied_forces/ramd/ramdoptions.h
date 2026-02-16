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
 * Declares options for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDOPTIONS_H
#define GMX_APPLIED_FORCES_RAMDOPTIONS_H

#include <string>

#include "gromacs/mdtypes/imdpoptionprovider.h"

#include "ramdparameters.h"

namespace gmx
{

class IndexGroupsAndNames;
class KeyValueTreeObject;
class KeyValueTreeBuilder;
class MDLogger;

/*! \internal
 * \brief Input data storage for RAMD
 */
class RAMDOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for RAMD to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     * \note This should be symmetrical to option initialization without
     *       employing manual prefixing with the section name string once
     *       the legacy code blocking this design is removed.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Report if this set of options is active
    bool active() const;

    //! Get parameters_ instance
    const RAMDParameters& parameters();

    //! Set the MDLogger instance
    void setLogger(const MDLogger& logger);

    //! Get the logger instance
    const MDLogger& logger() const;
    //! Set atom groups
    void setInputGroupIndices(const IndexGroupsAndNames&);

    //! Store the parameters that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

private:
    void readConfigString();

    //! Parameter values for force & energy evaluation
    RAMDParameters parameters_;

    /*! \brief MDLogger during preprocessing
     *
     * This is a pointer only because we need an "optional reference"
     * to a const MDLogger before the notification always provides the
     * actual reference. */
    const MDLogger* logger_ = nullptr;

    //! RAMD groups input file
    std::string groupsFile_;

    //! Content of the RAMD groups file
    std::string groupsString_;
};

} // namespace gmx

#endif
