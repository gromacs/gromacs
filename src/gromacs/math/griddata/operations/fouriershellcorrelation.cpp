/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * Implements helper class for autocorrelation tests
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "fouriershellcorrelation.h"
#include "fouriertransform.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/griddata/griddata.h"

#include <set>
#include <functional>
#include <numeric>
#include <algorithm>

namespace gmx
{

FourierShellCorrelation::FourierShellCorrelation(const GridWithTranslation<DIM> &realGrid)
{
    auto        reciprocalGrid = convertGridToReciprocalSpace(realGrid);
    const auto &unitcell       =  reciprocalGrid->unitCell();
    auto        spacing        = 2*std::max({unitcell.basisVectorLength(XX), unitcell.basisVectorLength(YY), unitcell.basisVectorLength(ZZ)});
    const auto &cell           = realGrid.cell();
    auto        highestK       = norm(RVec {cell.basisVectorLength(XX), cell.basisVectorLength(YY), cell.basisVectorLength(ZZ)});
    for (real binEdge = 0; binEdge < highestK +  spacing; binEdge += spacing)
    {
        binEdges_.insert(binEdge);
    }
    allocateShellDataContainersFromBins_(binEdges_);
};

void
FourierShellCorrelation::allocateShellDataContainersFromBins_(const std::set<real> &binEdges)
{
    for (const auto binEdge : binEdges)
    {
        referenceShells_[binEdge] = std::vector<t_complex>();
        otherShells_[binEdge]     = std::vector<t_complex>();
    }
}
/*! \brief Calculate fourier shells with custom binning. */
FourierShellCorrelation::FourierShellCorrelation(const std::set<real> &binEdges)
{
    allocateShellDataContainersFromBins_(binEdges);
};

const std::set<real> &
FourierShellCorrelation::getBinEdges() const
{
    return binEdges_;
};

class
FourierShellCorrelation::BinShells_
{
    public:
        BinShells_(fourierShell &shell) : shell_ {shell}
        {};
        void operator() (const t_complex &value, RVec k)
        {
            auto currentBin = shell_.lower_bound(norm(k));
            if (currentBin != shell_.end())
            {
                currentBin->second.push_back(value);
            }
        };
    private:
        fourierShell &shell_;
};

std::vector<real>
FourierShellCorrelation::getFscCurve(const GridDataReal3D &reference, const GridDataReal3D &other)
{
    for (auto &shell : referenceShells_)
    {
        shell.second.clear();
    }
    for (auto &shell : otherShells_)
    {
        shell.second.clear();
    }
    auto referenceFT = FourierTransformRealToComplex3D(reference).normalize().result();
    auto otherFT     = FourierTransformRealToComplex3D(other).normalize().result();

    ApplyToUnshiftedFourierTransform(*referenceFT).apply(BinShells_(referenceShells_));
    ApplyToUnshiftedFourierTransform(*otherFT).apply(BinShells_(otherShells_));

    std::vector<real> fscCurve;
    auto              otherShellIterator = std::begin(otherShells_);
    for (auto &referenceShell : referenceShells_)
    {
        fscCurve.push_back(correlateComplex_(referenceShell.second, otherShellIterator->second));
        ++otherShellIterator;
    }
    return fscCurve;
};

real
FourierShellCorrelation::correlateComplex_(const std::vector<t_complex> &a, const std::vector<t_complex> &b) const
{
    auto sumAbsoluteValues    =  [](const real &accumulat, const t_complex &value){
            return accumulat + value.re*value.re + value.im*value.im;
        };
    auto normASquare          = std::accumulate(a.begin(), a.end(), 0., sumAbsoluteValues);
    auto normBSquare          = std::accumulate(b.begin(), b.end(), 0., sumAbsoluteValues);
    real mulSum               = 0.;
    auto bIterator            = std::begin(b);
    for (auto aValue : a)
    {
        mulSum += aValue.re * bIterator->re + aValue.im * bIterator->im;
        ++bIterator;
    }
    return mulSum/sqrt(normASquare*normBSquare);
};


} /* gmx */
