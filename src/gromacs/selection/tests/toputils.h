/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Helper routines for constructing topologies for tests.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_TESTS_TOPUTILS_H
#define GMX_SELECTION_TESTS_TOPUTILS_H

#include <vector>

struct t_topology;
struct t_trxframe;

namespace gmx
{
namespace test
{

class TopologyManager
{
    public:
        TopologyManager();
        ~TopologyManager();

        void requestFrame();
        void requestVelocities();
        void requestForces();

        void loadTopology(const char *filename);
        void initAtoms(int count);
        void initAtomTypes(int count, const char *const types[]);
        template <int count>
        void initAtomTypes(const char *const (&types)[count])
        {
            initAtomTypes(count, types);
        }
        void initUniformResidues(int residueSize);
        void initUniformMolecules(int moleculeSize);

        t_topology *topology() { return top_; }
        t_trxframe *frame() { return frame_; }

    private:
        t_topology             *top_;
        t_trxframe             *frame_;
        std::vector<char *>     atomtypes_;
};

} // namespace test
} // namespace gmx

#endif
