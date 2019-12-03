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
 * This implements basic nblib tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/box.h"
#include "gromacs/nblib/coords.h"
#include "gromacs/nblib/interactions.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/nbkerneloptions.h"
#include "gromacs/nblib/nbkernelsystem.h"
#include "gromacs/nblib/topology.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, BasicEndToEndTest)
{
    int             sizeFactor = 1;
    NBKernelSystem  kernelSystem(sizeFactor);
    NBKernelOptions kernelOptions;
    kernelOptions.numIterations = 2;
    EXPECT_NO_THROW(nbKernel(kernelSystem, kernelOptions, false));
}

// TEST(NBlibTest, CoordinatesChange)
// {
//     int             sizeFactor = 1;
//     NBKernelSystem  kernelSystem(sizeFactor);
//     real            atom1InitialX = kernelSystem.coordinates[0][0];
//     NBKernelOptions kernelOptions;
//     kernelOptions.numIterations = 2;
//     nbKernel(kernelSystem, kernelOptions, false);
//     real atom1FinalX = kernelSystem.coordinates[0][0];
//     EXPECT_NE(atom1InitialX, atom1FinalX);
// }

TEST(NBlibTest, BasicArgonSetupTest)
{
    constexpr int NArgonAtoms = 100;

    nblib::Atom argonAtom("AR", 39.94800, 0.0, 0.0062647225, 9.847044e-06);

    Molecule argonMolecule("AR");
    argonMolecule.addAtom("AR", argonAtom);

    TopologyBuilder topBuilder;

    topBuilder.add(argonMolecule, NArgonAtoms);
    Topology top = topBuilder.buildTopology();

    Box box(7.73950);

    std::vector<gmx::RVec> coords = {
        { 5.158, 6.923, 3.413 },
        { 2.891, 6.634, 0.759 },
        { 4.356, 2.932, 1.414 },
        { 0.902, 5.872, 1.938 },
        { 0.074, 1.594, 7.358 },
        { 6.174, 1.962, 3.891 },
        { 0.444, 7.325, 4.901 },
        { 3.431, 3.521, 1.711 },
        { 2.761, 6.224, 4.937 },
        { 5.247, 3.657, 2.433 },
        { 3.561, 4.047, 0.460 },
        { 1.140, 6.853, 4.265 },
        { 1.336, 0.150, 7.228 },
        { 0.423, 2.200, 5.378 },
        { 4.364, 3.673, 1.201 },
        { 5.739, 1.297, 7.467 },
        { 6.919, 2.499, 3.090 },
        { 0.885, 2.419, 2.972 },
        { 3.895, 1.101, 7.121 },
        { 2.204, 0.213, 4.006 },
        { 0.556, 6.842, 7.265 },
        { 3.853, 0.558, 6.386 },
        { 7.174, 3.164, 5.165 },
        { 6.371, 0.064, 0.045 },
        { 7.400, 3.669, 6.568 },
        { 2.838, 3.889, 2.602 },
        { 3.497, 0.210, 2.784 },
        { 6.672, 4.529, 7.095 },
        { 1.219, 5.564, 0.465 },
        { 1.977, 2.186, 4.723 },
        { 5.414, 7.329, 3.934 },
        { 0.812, 3.466, 1.740 },
        { 1.683, 3.403, 4.855 },
        { 1.339, 7.445, 0.650 },
        { 7.174, 0.334, 3.986 },
        { 3.055, 2.214, 3.627 },
        { 0.667, 4.898, 3.606 },
        { 1.145, 2.904, 1.329 },
        { 2.606, 4.578, 1.721 },
        { 0.401, 5.614, 3.721 },
        { 1.227, 3.986, 6.081 },
        { 2.644, 2.238, 5.904 },
        { 1.401, 0.689, 6.317 },
        { 3.728, 5.849, 7.589 },
        { 3.476, 2.661, 2.815 },
        { 5.531, 4.391, 7.161 },
        { 6.887, 7.395, 4.518 },
        { 5.041, 7.280, 4.935 },
        { 3.676, 7.469, 3.273 },
        { 0.753, 6.868, 6.012 },
        { 2.346, 5.137, 0.896 },
        { 7.459, 7.701, 0.700 },
        { 3.381, 3.645, 5.993 },
        { 6.203, 4.601, 7.663 },
        { 2.189, 4.085, 1.251 },
        { 0.292, 7.132, 0.637 },
        { 0.837, 1.907, 6.153 },
        { 5.962, 4.325, 3.423 },
        { 1.707, 0.973, 2.160 },
        { 6.883, 6.808, 6.886 },
        { 6.899, 2.441, 0.083 },
        { 7.027, 6.668, 1.093 },
        { 6.033, 4.569, 5.681 },
        { 2.178, 2.946, 1.847 },
        { 1.298, 4.898, 1.441 },
        { 2.784, 2.323, 7.119 },
        { 6.814, 4.721, 2.273 },
        { 6.457, 5.465, 6.555 },
        { 7.063, 0.978, 0.708 },
        { 3.987, 3.185, 0.452 },
        { 5.052, 1.942, 4.802 },
        { 5.272, 1.323, 1.031 },
        { 2.890, 5.216, 6.061 },
        { 6.226, 4.075, 5.097 },
        { 5.526, 2.754, 3.738 },
        { 3.416, 1.297, 1.701 },
        { 6.658, 1.703, 6.706 },
        { 7.040, 6.544, 4.783 },
        { 0.455, 1.126, 1.969 },
        { 2.224, 0.683, 6.267 },
        { 5.864, 6.840, 5.273 },
        { 5.981, 1.779, 0.599 },
        { 5.271, 7.664, 1.537 },
        { 6.163, 4.723, 4.759 },
        { 6.354, 7.006, 6.059 },
        { 4.931, 1.359, 7.123 },
        { 1.426, 0.773, 0.122 },
        { 1.889, 5.085, 2.465 },
        { 5.432, 0.349, 3.569 },
        { 5.767, 1.071, 1.590 },
        { 4.643, 0.077, 2.502 },
        { 2.945, 3.368, 3.480 },
        { 5.821, 2.840, 2.285 },
        { 7.418, 1.903, 2.961 },
        { 0.481, 5.155, 1.421 },
        { 2.010, 7.343, 4.368 },
        { 1.278, 0.490, 5.614 },
        { 4.012, 1.630, 6.454 },
        { 5.515, 4.780, 4.228 },
        { 5.468, 3.128, 5.651 },
    };

    SimState simState(coords, box, top);

    EXPECT_EQ(top.getMasses().size(), NArgonAtoms);
    // EXPECT_EQ(top.getCharges().size(), NATOMS);
    // EXPECT_EQ(top.getMasses().size(), NATOMS);

    // ...
}

// TEST(NBlibTest, BasicWaterSetupTest)
// {
//     constexpr int NWaterMolecules = 100;

//     Atom oxygenAtom("OW", 16, -0.6, 1.0, 1.0);
//     Atom hydrogenAtom("HW", 1, 0.3, 0, 0);

//     Molecule waterMolecule("HOH");
//     waterMolecule.addAtom("O", oxygenAtom);
//     waterMolecule.addAtom("H1", hydrogenAtom);
//     waterMolecule.addAtom("H2", hydrogenAtom);

//     TopologyBuilder topBuilder();

//     top.add(waterMolecule, NWaterMolecules);
//     Topology top.buildTopology();

//     Box box(3.0);

//     // TODO: generate coords
//     std::vector<gmx::RVec> coords = {

//     };
//     EXPECT_EQ(coords().size(), top.numAtoms);

//     SimState simState(coord, box, topo);

//     EXPECT_EQ(top.simState().size(), NWaterMolecules);
//     // EXPECT_EQ(top.getCharges().size(), NATOMS);
//     // EXPECT_EQ(top.getMasses().size(), NATOMS);

//     // ...
// }

}  // namespace
}  // namespace test
}  // namespace nblib
