/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/gmxana/gmx_ana.h"

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionfiletype.h"

#include "gromacs/math/griddata/griddata.h"

#include "gromacs/fileio/griddataview.h"
#include "gromacs/fileio/griddataio.h"

#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/container/mask.h"
#include "gromacs/math/griddata/operations/gridinterpolator.h"

namespace gmx
{

namespace
{

//! \brief TODO
float maskedNormalizedInnerProduct(const GridDataReal3D::container_type &vec1,
                                   const GridDataReal3D::container_type &vec2,
                                   GridDataReal3D::value_type            mean1,
                                   GridDataReal3D::value_type            mean2,
                                   const Mask                           &mask)
{
    double     innerProduct = 0;
    const auto minMapSize   = std::min(vec1.size(), vec2.size());
#pragma omp parallel shared(innerProduct)
    {

#pragma omp for reduction(+:innerProduct)
        for (size_t l = 0; l < minMapSize; l++)
        {
            if (mask[l])
            {
                innerProduct += (vec1[l]-mean1) * (vec2[l]-mean2);
            }
        }
    }

    return innerProduct;
}

GridDataReal3D::value_type maskedMean(const GridDataReal3D::container_type &v, const Mask &mask)
{
    double sum = 0.;
#pragma omp parallel shared(sum)
    {
#pragma omp for reduction(+:sum)
        for (size_t l = 0; l < v.size(); ++l)
        {
            if (mask[l])
            {
                sum += v[l];
            }
        }
    }
    return sum/v.size();
}

//! \brief TODO
float correlationCoefficient(const GridDataReal3D::container_type &map1,
                             const GridDataReal3D::container_type &map2,
                             const Mask                           &mask)
{
    const auto mean1 = maskedMean(map1, mask);
    const auto mean2 = maskedMean(map2, mask);
    return maskedNormalizedInnerProduct(map1, map2, mean1, mean2, mask) / sqrt(maskedNormalizedInnerProduct(map2, map2, mean2, mean2, mask) *
                                                                               maskedNormalizedInnerProduct(map1, map1, mean1, mean1, mask));
}


//! \brief Calculate the local correlation
GridDataReal3D localCorrelationMap(const GridDataReal3D &map1, const GridDataReal3D &map2, const Mask &mask, float range)
{
    auto localCorrelationMap = map1;
    // Initialize the output localcorrelationmap map, setting the borders and masked points,
    // with non-assigned NAN values
    auto maxDistX = static_cast<int>(range / map1.getGrid().unitCell().basisVectorLength(XX)) + 1;
    auto maxDistY = static_cast<int>(range / map1.getGrid().unitCell().basisVectorLength(YY)) + 1;
    auto maxDistZ = static_cast<int>(range / map1.getGrid().unitCell().basisVectorLength(ZZ)) + 1;

    std::fill(std::begin(localCorrelationMap), std::end(localCorrelationMap), nanf(""));
    const auto &lattice   = map1.getGrid().lattice();
    for (size_t i = 0; i < map1.size(); ++i)
    {
        auto  latticeIndex = lattice.vectoriseLinearIndex(i);
        float sum12        = 0.0;
        float sum11        = 0.0;
        float sum22        = 0.0;
        float sum1         = 0.;
        float sum2         = 0.;

        for (int iZ = latticeIndex[ZZ] - maxDistZ; iZ <= latticeIndex[ZZ] + maxDistZ; iZ++)
        {
            for (int iY = latticeIndex[YY] - maxDistY; iY <= latticeIndex[YY] + maxDistY; iY++)
            {
                for (int iX = latticeIndex[XX] - maxDistX; iX <= latticeIndex[XX] + maxDistX; iX++)
                {
                    ColumnMajorLattice<DIM>::MultiIndex localSurroundingIndex({{iX, iY, iZ}});
                    if (lattice.inLattice(localSurroundingIndex) && mask[i])
                    {
                        const auto v1 = *(map1.iteratorAtMultiIndex(localSurroundingIndex));
                        const auto v2 = *(map2.iteratorAtMultiIndex(localSurroundingIndex));
                        sum12 += v1 * v2;
                        sum11 += v1 * v1;
                        sum22 += v2 * v2;
                        sum1  += v1;
                        sum2  += v2;
                    }
                }
            }
        }

        *(localCorrelationMap.iteratorAtMultiIndex(latticeIndex)) = (sum12 - sum1 * sum2) / sqrt((sum22-sum2*sum2) * (sum11-sum1*sum1));
    }
    return localCorrelationMap;
}
}
int gmx_mapcompare(int argc, char *argv[])
{
    std::string desc =
        "[THISMODULE] calculates the difference between two density maps.";

    std::vector<const char *> bugs =
    {"This tool is under construction.", ""};
    Options                   options;

    std::string               fnInput;
    options.addOption(StringOption("refmap").store(&fnInput));
    std::string               fnCompare;
    options.addOption(StringOption("compare").store(&fnCompare));

    std::string               fnLocalCorrelation = "localcorrelation.ccp4";
    options.addOption(FileNameOption("localcorrelation").filetype(eftCCP4).outputFile().defaultBasename("localcorrelation").description("Local correlation between computed and reference map").store(&fnLocalCorrelation));

    float rLocal = 0.3;
    options.addOption(FloatOption("rlocal").store(&rLocal)
                          .description("Local correlation range"));

    std::string               fnMask;
    options.addOption(StringOption("mask").store(&fnMask).required(false));
    float                     maskCutOff = 0.01;
    options.addOption(FloatOption("maskcutoff").store(&maskCutOff).required(false));

    gmx::CommandLineParser(&options).parse(&argc, argv);
    options.finish();

    const CommandLineHelpContext *context = GlobalCommandLineHelpContext::get();
    if (context != nullptr)
    {
        CommandLineHelpWriter(options).setHelpText(desc)
            .setKnownIssues(arrayRefFromArray(bugs.data(), bugs.size()))
            .writeHelp(*context);
    }
    Mask mask;
    if (!fnMask.empty())
    {
        mask = Mask(MrcFile().read(fnMask), maskCutOff);
    }
    const auto map        = MapConverter(fnInput).map();
    auto       compareMap = MapConverter(fnCompare).map();

    // If maps have different grids, interpolate the compare map onto the grid of the input map
    if (!(map.getGrid() == compareMap.getGrid()))
    {
        const auto &interpolatedMap = interpolateLinearly(compareMap, map.getGrid());
        compareMap = interpolatedMap;
        fprintf(stderr, "\nWARNING: Map grids don't match. Performing linear interpolation between grids.\n");
    }

    fprintf(stdout, "Map correlation coefficient = %g .", correlationCoefficient(map, compareMap, mask));

    if (!fnLocalCorrelation.empty())
    {
        auto lcMap = localCorrelationMap(map, compareMap, mask, rLocal);
        mask.maskToNan(&lcMap);
        MrcFile().write(fnLocalCorrelation, lcMap);
    }
    return 0;
}
}
