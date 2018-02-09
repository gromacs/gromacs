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
#ifndef GMX_MATH_FOURIERSHELLCORRELATION_H_
#define GMX_MATH_FOURIERSHELLCORRELATION_H_

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/gmxcomplex.h"
#include <set>
#include <map>
#include <vector>

namespace gmx
{

class FourierShellCorrelation
{
    public:
        typedef std::map < real, std::vector < t_complex>> fourierShell;
        FourierShellCorrelation() = default;
        /*! \brief Set bins from real-space grid guaranteeing six datapoints per shell.
         *
         */
        FourierShellCorrelation(const GridWithTranslation<DIM> &RealGrid);
        /*! \brief Calculate fourier shells with custom binning. */
        FourierShellCorrelation(const std::set<real> &binEdges);
        const std::set<real> &getBinEdges() const;
        std::vector<real> getFscCurve(const GridDataReal3D &reference, const GridDataReal3D &other);

    private:
        class BinShells_;
        void allocateShellDataContainersFromBins_(const std::set<real> &binEdges);
        real correlateComplex_(const std::vector<t_complex> &a, const std::vector<t_complex> &b) const;
        std::set<real> binEdges_;
        fourierShell   referenceShells_;
        fourierShell   otherShells_;

};

}
#endif /* end of include guard: GMX_MATH_FOURIERSHELLCORRELATION_H_ */
