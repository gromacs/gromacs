/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_DIMETHYL_BUTENE_DATA_HPP
#define NBLIB_DIMETHYL_BUTENE_DATA_HPP

#include "nblib/listed_forces/traits.h"
#include "nblib/topologyhelpers.h"

namespace nblib
{

/*! \brief sets up an interaction tuple for a poly-dimethyl-butene chain
 *
 *                 ( gamma_i   delta_i            ) gamma_N   delta_N
 *                 (        \ /                   )        \ /
 *     alpha_-1 -- (     -- beta_i -- alpha_i --  )     -- beta_N -- alpha_N
 *                 (                   / \        )                   / \
 *                 (          epsilon_i   zeta_i  )__N       epsilon_N   zeta_N
 *
 */
class PolyDimethylButene
{
public:
    //! \brief creates a polymer of length <length>, where <length> is N in the diagram above
    explicit PolyDimethylButene(int length) :
        x0(0.0), y0(0.5), z0(0.5), phi(109.5), bondLength(0.01)
    {
        real planarAngle = 30; // ~(180 - phi) / 2
        a                = bondLength * std::cos(planarAngle * DEG2RAD);
        b                = bondLength * std::sin(planarAngle * DEG2RAD);
        setInteractionParameters();
        makePolymer(length);

        box = Box(3 * a * length, 1, 1);
    }

    //! \brief if this is called, listed interactions will be represented as aggregates where possible
    void createAggregates()
    { /* aggregateTransformations(interactions); */
    }

    std::vector<gmx::RVec> x;

    ListedInteractionData interactions;

    Box box{ 0 };

private:
    // x0,y0,z0: position of first atom
    real x0, y0, z0;
    // phi: ~tetrahedron angle
    real phi;

    real bondLength;

    // a: monomer length in X
    // b: monomer width in Y and height in Z
    real a, b;

    int paramIndex(int i) const { return i % 2; }

    void setInteractionParameters()
    {
        HarmonicBondType                bond1{ 37600, bondLength };
        HarmonicBondType                bond2{ 31380, bondLength };
        std::array<HarmonicBondType, 2> bondsp{ bond1, bond2 };

        HarmonicAngle                angle1{ 100, Degrees(phi) };
        HarmonicAngle                angle2{ 200, Degrees(phi) };
        std::array<HarmonicAngle, 2> anglesp{ angle1, angle2 };

        ProperDihedral                dihedral1{ Degrees(0), 40, 1 };
        ProperDihedral                dihedral2{ Degrees(0), 50, 1 };
        std::array<ProperDihedral, 2> dihedralsp{ dihedral1, dihedral2 };

        {
            auto& p = pickType<HarmonicBondType>(interactions).parametersA;
            p.insert(begin(p), begin(bondsp), end(bondsp));
        }
        {
            auto& p = pickType<HarmonicAngle>(interactions).parametersA;
            p.insert(begin(p), begin(anglesp), end(anglesp));
        }
        {
            auto& p = pickType<ProperDihedral>(interactions).parametersA;
            p.insert(begin(p), begin(dihedralsp), end(dihedralsp));
        }
    }

    void makePolymer(int length)
    {
        addHead();

        for (int i = 0; i < length; ++i)
        {
            addBulk(i);
        }

        addTail(length);
    }

    void addHead() { x.push_back(alpha(0)); }

    // add i-th monomer
    void addBulk(int i)
    {
        auto& basicInteractions = interactions;

        x.push_back(beta(i));
        x.push_back(gamma(i));
        x.push_back(delta(i));
        x.push_back(epsilon(i));
        x.push_back(zeta(i));
        x.push_back(alpha(i));

        auto& bonds = pickType<HarmonicBondType>(basicInteractions);

        // connect to previous monomer
        bonds.indices.push_back({ seqAlpha(i - 1), seqBeta(i), paramIndex(i) });
        // back-bone link
        bonds.indices.push_back({ seqBeta(i), seqAlpha(i), paramIndex(i) });
        // beta (upward) fins
        bonds.indices.push_back({ seqBeta(i), seqGamma(i), paramIndex(i) });
        bonds.indices.push_back({ seqBeta(i), seqDelta(i), paramIndex(i) });
        // alpha (downward) fins
        bonds.indices.push_back({ seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        bonds.indices.push_back({ seqAlpha(i), seqZeta(i), paramIndex(i) });

        auto& angles = pickType<HarmonicAngle>(basicInteractions);

        // beta angles
        angles.indices.push_back({ seqAlpha(i - 1), seqBeta(i), seqAlpha(i), paramIndex(i) });
        angles.indices.push_back({ seqAlpha(i - 1), seqBeta(i), seqGamma(i), paramIndex(i) });
        angles.indices.push_back({ seqAlpha(i - 1), seqBeta(i), seqDelta(i), paramIndex(i) });
        angles.indices.push_back({ seqGamma(i), seqBeta(i), seqDelta(i), paramIndex(i) });
        angles.indices.push_back({ seqAlpha(i), seqBeta(i), seqGamma(i), paramIndex(i) });
        angles.indices.push_back({ seqAlpha(i), seqBeta(i), seqDelta(i), paramIndex(i) });

        // alpha angles
        angles.indices.push_back({ seqBeta(i), seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        angles.indices.push_back({ seqBeta(i), seqAlpha(i), seqZeta(i), paramIndex(i) });
        angles.indices.push_back({ seqBeta(i), seqAlpha(i), seqBeta(i + 1), paramIndex(i) });
        angles.indices.push_back({ seqBeta(i + 1), seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        angles.indices.push_back({ seqBeta(i + 1), seqAlpha(i), seqZeta(i), paramIndex(i) });
        angles.indices.push_back({ seqEpsilon(i), seqAlpha(i), seqZeta(i), paramIndex(i) });

        auto& dihedrals = pickType<ProperDihedral>(basicInteractions);

        // beta_i -- alpha_i axis
        // alpha_-1 -- beta leg
        dihedrals.indices.push_back(
                { seqAlpha(i - 1), seqBeta(i), seqAlpha(i), seqBeta(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqAlpha(i - 1), seqBeta(i), seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqAlpha(i - 1), seqBeta(i), seqAlpha(i), seqZeta(i), paramIndex(i) });
        // gamma -- beta leg
        dihedrals.indices.push_back(
                { seqGamma(i), seqBeta(i), seqAlpha(i), seqBeta(i + 1), paramIndex(i) });
        dihedrals.indices.push_back({ seqGamma(i), seqBeta(i), seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        dihedrals.indices.push_back({ seqGamma(i), seqBeta(i), seqAlpha(i), seqZeta(i), paramIndex(i) });
        // delta -- beta leg
        dihedrals.indices.push_back(
                { seqDelta(i), seqBeta(i), seqAlpha(i), seqBeta(i + 1), paramIndex(i) });
        dihedrals.indices.push_back({ seqDelta(i), seqBeta(i), seqAlpha(i), seqEpsilon(i), paramIndex(i) });
        dihedrals.indices.push_back({ seqDelta(i), seqBeta(i), seqAlpha(i), seqZeta(i), paramIndex(i) });

        // alpha_i -- beta_i+1 axis
        // beta_i -- alpha_i leg
        dihedrals.indices.push_back(
                { seqBeta(i), seqAlpha(i), seqBeta(i + 1), seqAlpha(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqBeta(i), seqAlpha(i), seqBeta(i + 1), seqGamma(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqBeta(i), seqAlpha(i), seqBeta(i + 1), seqDelta(i + 1), paramIndex(i) });
        // epsilon_i -- alpha_i leg
        dihedrals.indices.push_back(
                { seqEpsilon(i), seqAlpha(i), seqBeta(i + 1), seqAlpha(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqEpsilon(i), seqAlpha(i), seqBeta(i + 1), seqGamma(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqEpsilon(i), seqAlpha(i), seqBeta(i + 1), seqDelta(i + 1), paramIndex(i) });
        // zeta_i -- alpha_i leg
        dihedrals.indices.push_back(
                { seqZeta(i), seqAlpha(i), seqBeta(i + 1), seqAlpha(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqZeta(i), seqAlpha(i), seqBeta(i + 1), seqGamma(i + 1), paramIndex(i) });
        dihedrals.indices.push_back(
                { seqZeta(i), seqAlpha(i), seqBeta(i + 1), seqDelta(i + 1), paramIndex(i) });
    }

    void addTail(int n)
    {
        x.push_back(beta(n));
        x.push_back(gamma(n));
        x.push_back(delta(n));
        x.push_back(epsilon(n));
        x.push_back(zeta(n));
        x.push_back(alpha(n));
    }


    [[nodiscard]] gmx::RVec alpha(int n) const { return gmx::RVec{ x0 + a * (2 * n + 2), y0, z0 }; }

    [[nodiscard]] gmx::RVec beta(int n) const
    {
        return gmx::RVec{ x0 + a * (2 * n + 1), y0, z0 + b };
    }

    [[nodiscard]] gmx::RVec gamma(int n) const
    {
        return gmx::RVec{ x0 + a * (2 * n + 1), y0 + a, z0 + 2 * b };
    }

    [[nodiscard]] gmx::RVec delta(int n) const
    {
        return gmx::RVec{ x0 + a * (2 * n + 1), y0 - a, z0 + 2 * b };
    }

    [[nodiscard]] gmx::RVec epsilon(int n) const
    {
        return gmx::RVec{ x0 + a * (2 * n + 2), y0 + a, z0 - b };
    }

    [[nodiscard]] gmx::RVec zeta(int n) const
    {
        return gmx::RVec{ x0 + a * (2 * n + 2), y0 - a, z0 - b };
    }

    static int seqAlpha(int n) { return 6 * (n + 1); }

    static int seqBeta(int n) { return 6 * n + 1; }

    static int seqGamma(int n) { return 6 * n + 2; }

    static int seqDelta(int n) { return 6 * n + 3; }

    static int seqEpsilon(int n) { return 6 * n + 4; }

    static int seqZeta(int n) { return 6 * n + 5; }
};

} // namespace nblib

#endif // NBLIB_DIMETHYL_BUTENE_HPP
