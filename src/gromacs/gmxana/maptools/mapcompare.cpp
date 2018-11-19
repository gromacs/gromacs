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

#include "mapcompare.h"

#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/griddataview.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/container/mask.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/operations/gridinterpolator.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"

namespace gmx
{

class MapCompare final : public ICommandLineOptionsModule
{
    public:
        void init(CommandLineModuleSettings *settings) override;
        void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings) override;
        void optionsFinished() override;
        int run() override;
    private:
        static const char *const helpText_[];
        std::string              fnInput_;
        std::string              fnCompare_;
        std::string              fnLocalCorrelation_;
        std::string              fnMask_;
        float                    rLocal_;
        float                    maskCutOff_;
};

namespace
{

//! \brief TODO
float maskedNormalizedInnerProduct(const GridDataFloat3D::container_type &vec1,
                                   const GridDataFloat3D::container_type &vec2,
                                   GridDataFloat3D::value_type            mean1,
                                   GridDataFloat3D::value_type            mean2,
                                   const Mask                            &mask)
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

GridDataFloat3D::value_type maskedMean(const GridDataFloat3D::container_type &v, const Mask &mask)
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
float correlationCoefficient(const GridDataFloat3D::container_type &map1,
                             const GridDataFloat3D::container_type &map2,
                             const Mask                            &mask)
{
    const auto mean1 = maskedMean(map1, mask);
    const auto mean2 = maskedMean(map2, mask);
    return maskedNormalizedInnerProduct(map1, map2, mean1, mean2, mask) / sqrt(maskedNormalizedInnerProduct(map2, map2, mean2, mean2, mask) *
                                                                               maskedNormalizedInnerProduct(map1, map1, mean1, mean1, mask));
}


//! \brief Calculate the local correlation
GridDataFloat3D localCorrelationMap(const GridDataFloat3D &map1, const GridDataFloat3D &map2, const Mask &mask, float range)
{
    constexpr int minNumberGridPointsForLocalCorrelation = 0;
    auto          localCorrelationMap                    = map1;
    containeroperation::setZero(&localCorrelationMap);
    // Initialize the output localcorrelationmap map, setting the borders and masked points,
    // with non-assigned NAN values
    auto        maxDistX = static_cast<int>(range / map1.getGrid().unitCell()[XX]) + 1;
    auto        maxDistY = static_cast<int>(range / map1.getGrid().unitCell()[YY]) + 1;
    auto        maxDistZ = static_cast<int>(range / map1.getGrid().unitCell()[ZZ]) + 1;

    const auto &lattice      = map1.getGrid().lattice();
    auto        latticeIndex = lattice.begin();
    for (size_t i = 0; i < map1.size(); ++i, ++latticeIndex)
    {
        double sum12        = 0.0;
        double sum11        = 0.0;
        double sum22        = 0.0;
        double sum1         = 0.;
        double sum2         = 0.;
        int    n            = 0;

        for (int iZ = (*latticeIndex)[ZZ] - maxDistZ; iZ <= (*latticeIndex)[ZZ] + maxDistZ; iZ++)
        {
            for (int iY = (*latticeIndex)[YY] - maxDistY; iY <= (*latticeIndex)[YY] + maxDistY; iY++)
            {
                for (int iX = (*latticeIndex)[XX] - maxDistX; iX <= (*latticeIndex)[XX] + maxDistX; iX++)
                {
                    const offset<DIM> localSurroundingIndex = {iX, iY, iZ};
                    if (lattice.contains(localSurroundingIndex) && mask[i])
                    {
                        const auto v1 = map1[localSurroundingIndex];
                        const auto v2 = map2[localSurroundingIndex];
                        sum11 += v1 * v1;
                        sum22 += v2 * v2;

                        sum12 += v1 * v2;
                        sum1  += v1;
                        sum2  += v2;
                        ++n;
                    }
                }
            }
        }
        if (n > minNumberGridPointsForLocalCorrelation)
        {
            localCorrelationMap[*latticeIndex] = (n * sum12 - sum1 * sum2) / sqrt( (n * sum11 - sum1 * sum1) * (n * sum22 - sum2 * sum2) );
        }
    }
    return localCorrelationMap;
}

}   // anomymous namespace

void MapCompare::init(CommandLineModuleSettings * /*settings*/)
{}

void MapCompare::optionsFinished()
{}

const char *const MapCompare::helpText_[] = {
    "[THISMODULE] compares two density maps by calculating their"
    " cross-correlation coefficient. Map regions may be excluded"
    " from comparison using a [TT]-mask[tt] map, wich may be the"
    " same as one of the input maps. All values in the mask below"
    " [TT]-maskcutoff[tt] are excluded from comparison. If the"
    " grid of the reference and the comparison map don't match,"
    " the comparison map will be interpolated linearly onto the"
    " reference map grid. The [TT]-localcorrelation[tt] outputs"
    " a local correlation map on the same grid as the reference"
    " map where each voxels reports the correlation along a sliding"
    " box of length [TT]-rlocal[tt], which is centered at the voxel"
    " and cut off at the map boudaries."
};

void MapCompare::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings)
{

    settings->setHelpText(helpText_);

    options->addOption(StringOption("refmap")
                           .store(&fnInput_)
                           .defaultValue("reference.ccp4")
                           .required());

    options->addOption(StringOption("compare")
                           .store(&fnCompare_)
                           .defaultValue("compare.ccp4")
                           .required());

    options->addOption(FileNameOption("localcorrelation")
                           .filetype(eftCCP4)
                           .outputFile()
                           .defaultBasename("localcorrelation")
                           .description("Local correlation between computed and reference map")
                           .store(&fnLocalCorrelation_));

    options->addOption(FloatOption("rlocal")
                           .defaultValue(0.3)
                           .store(&rLocal_)
                           .description("Local correlation range"));

    options->addOption(StringOption("mask")
                           .store(&fnMask_) );

    options->addOption(FloatOption("maskcutoff")
                           .defaultValue(1e-5)
                           .store(&maskCutOff_) );

}

int MapCompare::run()
{
    Mask mask;
    if (!fnMask_.empty())
    {
        mask = Mask(MrcFile().read(fnMask_), maskCutOff_);
    }
    const auto map        = MapConverter(fnInput_).map();
    auto       compareMap = MapConverter(fnCompare_).map();

    // If maps have different grids, interpolate the compare map onto the grid of the input map
    if (!(map.getGrid() == compareMap.getGrid()))
    {
        const auto &interpolatedMap = interpolateLinearly(compareMap, map.getGrid());
        compareMap = interpolatedMap;
        fprintf(stderr, "\nWARNING: Map grids don't match. Performing linear interpolation between grids.\n");
    }

    fprintf(stdout, "Map correlation coefficient = %g .", correlationCoefficient(map, compareMap, mask));

    if (!fnLocalCorrelation_.empty())
    {
        auto lcMap = localCorrelationMap(map, compareMap, mask, rLocal_);
        mask.maskToValue(&lcMap, 0.);
        MrcFile().write(fnLocalCorrelation_, lcMap);
    }
    return EXIT_SUCCESS;
}

const char mapcompareInfo::name[]             = "mapcompare";
const char mapcompareInfo::shortDescription[] = "Compare density maps";
ICommandLineOptionsModulePointer mapcompareInfo::create()
{
    return compat::make_unique<MapCompare>();
}


}
