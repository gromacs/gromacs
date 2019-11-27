/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements nblib Topology and TopologyBuilder
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/topology/block.h"

#include "molecules.h"

struct t_blocka;

namespace nblib {

class Topology {
public:

    const std::vector<int>& getAtoms() const;

    const std::vector<real>& getCharges() const;

    const std::vector<real>& getMasses() const;

    const std::vector<real>& getNonbondedParameters() const;

    const std::vector<int>& getAtomInfoAllVdw() const;

    // TODO: This function is only needed for testing. Need
    //       another way for testing exclusion correctness
    const t_blocka& getGMXexclusions() const
    {
        return excls;
    }

private:
    Topology() = default;

    friend class TopologyBuilder;

    //! Storage for parameters for short range interactions.
    std::vector<real>      nonbondedParameters;
    //! Storage for atom type parameters.
    std::vector<int>       atomTypes;
    //! Storage for atom partial charges.
    std::vector<real>      charges;
    //! Atom masses
    std::vector<real>      masses;
    //! Atom info where all atoms are marked to have Van der Waals interactions
    std::vector<int>       atomInfoAllVdw;
    //! Information about exclusions.
    t_blocka               excls;
};

class TopologyBuilder {
public:
    TopologyBuilder();

    Topology buildTopology();

    TopologyBuilder& addMolecule(Molecule moleculeType, int nMolecules);

private:
    Topology topology_;

    int numAtoms_;
    std::vector<std::tuple<Molecule, int>> molecules_;

    t_blocka createExclusionsList() const;

    template <class Extractor>
    std::vector<real> extractAtomTypeQuantity(Extractor extractor);

    std::vector<real> extractCharge();
};

} // namespace nblib

#endif //GROMACS_TOPOLOGY_H
