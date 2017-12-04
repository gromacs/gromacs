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

}

int gmx_mapcompare(int argc, char *argv[])
{
    std::string desc =
        "[THISMODULE] calculates the difference between two density maps.";

    std::vector<const char *> bugs =
    {"This tool is under construction.", ""};
    Options                   options;

    std::string               fnInput;
    options.addOption(StringOption("mi").store(&fnInput));
    std::string               fnCompare;
    options.addOption(StringOption("compare").store(&fnCompare));
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
    auto map        = MapConverter(fnInput).map();
    auto compareMap = MapConverter(fnCompare).map();

    fprintf(stdout, "Map correlation coefficient = %g .", correlationCoefficient(map, compareMap, mask));

    return 0;
}

}
