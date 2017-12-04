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

#include "densfitdata.h"

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/mrcmetadata.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/txtdump.h"

namespace gmx
{



namespace
{
const char *densityPotentialNames[static_cast<int>(DensityPotential::NumberOfPotentials) +1] = { "cross-correlation", "cross-entropy", nullptr };
}

bool TimeDependentDensfitParameters::vary() const
{
    return timePoints_.size() > 1;
}

int TimeDependentDensfitParameters::str_nelem(const char *str, int maxptr, char *ptr[])
{
    int   np = 0;
    char *copy0, *copy;

    copy0 = gmx_strdup(str);
    copy  = copy0;
    ltrim(copy);
    while (*copy != '\0')
    {
        if (np >= maxptr)
        {
            gmx_fatal(FARGS, "Too many groups on line: '%s' (max is %d)",
                      str, maxptr);
        }
        if (ptr)
        {
            ptr[np] = copy;
        }
        np++;
        while ((*copy != '\0') && !isspace(*copy))
        {
            copy++;
        }
        if (*copy != '\0')
        {
            *copy = '\0';
            copy++;
        }
        ltrim(copy);
    }
    if (ptr == nullptr)
    {
        sfree(copy0);
    }

    return np;

}

TimeDependentDensfitParameters::TimeDependentDensfitParameters(
        const char timePointString[STRLEN], const char forceConstantString[STRLEN],
        const char sigmaString[STRLEN], warninp *wi)
{
    #define MAXPTR 254

    /* Get values for the refinement points from input strings */
    /* Read time values from input string*/
    char *ptr1[MAXPTR];

    int   nRead = str_nelem(timePointString, MAXPTR, ptr1);
    for (int i = 0; i < nRead; i++)
    {
        timePoints_.push_back(strtod(ptr1[i], nullptr));
    }
    if (nRead == 0)
    {
        timePoints_.push_back(0.);
    }
    if (!std::is_sorted(std::begin(timePoints_), std::end(timePoints_)))
    {
        gmx_fatal(
                FARGS,
                "Density fitting time points must be given in increasing order!\n");
    }

    /* Read sigma values from input string */

    nRead = str_nelem(sigmaString, MAXPTR, ptr1);
    for (int i = 0; i < nRead; i++)
    {
        sigmas_.push_back(strtod(ptr1[i], nullptr));
    }

    /* Read k values from input string */
    nRead = str_nelem(forceConstantString, MAXPTR, ptr1);
    for (int i = 0; i < nRead; i++)
    {
        forceConstants_.push_back(strtod(ptr1[i], nullptr));
        if (forceConstants_.back() < 0)
        {
            warning(wi, "Negative density fitting force constants expel atoms from the map.");
        }
    }
}

real TimeDependentDensfitParameters::currentForceConstant(real time) const
{
    return interpolateLinearly(time, timePoints_, forceConstants_);
}
real TimeDependentDensfitParameters::currentSigma(real time) const
{
    return interpolateLinearly(time, timePoints_, sigmas_);
}

//! \brief TODO
real TimeDependentDensfitParameters::interpolateLinearly(
        real x0, const std::vector<real> &x, const std::vector<real> &y) const
{
    // Assuming sorted x, find the first element larger or equal x0
    const auto xLowerBound     = std::lower_bound(std::begin(x), std::end(x), x0);
    const auto lowerBoundIndex = std::distance(std::begin(x), xLowerBound);

    // If x is larger than any x_value or its x's lower boundary lacks a
    // corresponding y value,
    // interpolate with a constant that is the last y value
    if ((xLowerBound == std::end(x)) ||
        (lowerBoundIndex >= static_cast<decltype(lowerBoundIndex)>(y.size())))
    {
        return y.back();
    }

    // If x is smaller than any x_value, interpolate with a constant that is the
    // first y value
    if (xLowerBound == std::begin(x))
    {
        return y.front();
    }

    /* for x0 between smallest and largest x, interpolate y linearly
     * as y = r * y(lowerboundary(x0)) + (r-1) * (y(beforelowerboundary(x0)))
     * r is the relative distance in the interval from next smaller
     * x value from x0s the lower boundary to x0's lower boundary */
    const auto r =
        (x0 - *(xLowerBound - 1)) / (*xLowerBound - *(xLowerBound - 1));
    return r * y[lowerBoundIndex] + (1 - r) * y[lowerBoundIndex - 1];
}

void TimeDependentDensfitParameters::print(FILE *fp, int indent) const
{
    pr_reals(fp, indent, "densfit-time ", timePoints_.data(), timePoints_.size());
    pr_reals(fp, indent, "densfit-sigma", sigmas_.data(), sigmas_.size());
    pr_reals(fp, indent, "densfit-k    ", forceConstants_.data(),
             forceConstants_.size());
}

void TimeDependentDensfitParameters::serialize(ISerializer *  serializer)
{

    int npoints = 0;
    if (!serializer->reading())
    {
        npoints = timePoints_.size();
    }
    serializer->doInt(&npoints);

    if (serializer->reading())
    {
        timePoints_.resize(npoints);
        forceConstants_.resize(npoints);
        sigmas_.resize(npoints);
    }
    for (auto &time : timePoints_)
    {
        serializer->doFloat(&time);
    }
    for (auto &forceConstant : forceConstants_)
    {
        serializer->doFloat(&forceConstant);
    }
    for (auto &sigma : sigmas_)
    {
        serializer->doFloat(&sigma);
    }

}

const TimeDependentDensfitParameters &DensfitData::timeDependent() const
{
    return timeDependent_;
}

void DensfitData::printToOut(FILE *fp, const gmx_output_env_t *oenv) const
{
    fprintf(fp, "# Density fitting forces for %lu atoms are updated every %d "
            "time step%s.\n",
            group_.ind_.size(), nstfit_, nstfit_ != 1 ? "s" : "");
    fprintf(fp, "# Output is written in intervals of %d time step%s.\n", nstout_,
            nstout_ > 1 ? "s" : "");
    fprintf(fp, "# Parameter sigma_dist = %g\n", dist_);
    if (nstmapout_ != 0)
    {
        fprintf(fp, "# Calculated map will be saved every %d steps\n", nstmapout_);
    }
    std::vector<std::string> setname;
    std::stringstream        LegendStr;
    /* Print a nice legend */
    LegendStr << "#";
    LegendStr << std::setw(12)<< "time";

    setname.emplace_back("temperature (K)");
    LegendStr << std::setw(12)<< "T";
    setname.emplace_back("k (kJ/mol)");
    LegendStr << std::setw(12)<< "k";
    setname.emplace_back("sigma (nm)");
    LegendStr << std::setw(12) << "sigma";
    setname.emplace_back("correlation coefficient");
    LegendStr << std::setw(12) << "c.c.";
    setname.emplace_back("V_dens (kJ/mol)");
    LegendStr << std::setw(12)  << "V_dens";

    xvgrLegend(fp, setname, oenv);

    fprintf(fp, "#\n# Legend for the following data columns:\n");
    fprintf(fp, "%s\n", LegendStr.str().c_str());
}

void DensityFittingGroup::serialize(ISerializer * serializer)
{
    int nat = ind_.size();
    serializer->doInt(&nat);
    if (serializer->reading())
    {
        ind_.resize(nat);
    }
    for (auto &i : ind_)
    {
        serializer->doInt(&i);
    }
    serializer->doString(&name_);
}

const DensityFittingGroup &DensfitData::fittingGroup() const
{
    return group_;
}

void DensfitData::serialize(ISerializer *serializer)
{
    map_ref_.readwrite(serializer);
    timeDependent_.serialize(serializer);
    group_.serialize(serializer);

    serializer->doFloat(&dist_);
    serializer->doInt(&nstfit_);
    serializer->doInt(&nstout_);
    serializer->doInt(&nstmapout_);

    // Encode the densitypotential enmum class as int for serializing
    int densityPotentialEnumAsInt;
    if (!serializer->reading())
    {
        densityPotentialEnumAsInt = static_cast<int>(densityPotential_);
    }
    serializer->doInt(&densityPotentialEnumAsInt);
    if (serializer->reading())
    {
        densityPotential_ = DensityPotential(densityPotentialEnumAsInt);
    }

    serializer->doBool(&bKeepAndNumberMaps_);
    serializer->doBool(&normalizeMapsToUnity_);

}

bool DensfitData::keepAndNumberMaps() const { return bKeepAndNumberMaps_; };

bool DensfitData::outputMapThisStep(gmx_int64_t step) const
{
    return do_per_step(step, nstmapout_);
}

bool DensfitData::fitThisStep(gmx_int64_t step) const
{
    return do_per_step(step, nstfit_);
}
void DensfitData::print(FILE *fp, int indent) const
{
    fprintf(fp, "density fitting options:\n");

    timeDependent_.print(fp, indent);
    const auto densiftRefMapString =  MrcMetaData().fromGrid(map_ref_.getGrid()).set_grid_stats(map_ref_).to_string();
    pr_str(fp, indent, "densfit-referenceMap", densiftRefMapString.c_str());

    pr_real(fp, indent, "densfit-dist", dist_);
    pr_int(fp, indent, "densfit-nstfit", nstfit_);
    pr_int(fp, indent, "densfit-nstout", nstout_);
    pr_int(fp, indent, "densfit-nstmapout", nstmapout_);
    pr_str(fp, indent, "densfit-potential", densityPotentialNames[static_cast<int>(densityPotential_)]);
    pr_str(fp, indent, "densfit-bKeepMaps", boolToString(bKeepAndNumberMaps_));
    pr_str(fp, indent, "densfit-normalizemapstounity", boolToString(normalizeMapsToUnity_));
    pr_ivec_block(fp, indent, "atom", group_.ind_.data(), group_.ind_.size(), true);

}

const GridDataReal3D &DensfitData::referenceMap() const
{
    return map_ref_;
}

void DensfitData::makeGroups(t_blocka *grps, char **gnames)
{
    /* Group index from group name */
    auto ig = 0;

    for (auto i = 0; i < grps->nr; ++i)
    {
        if (gmx_strcasecmp(group_.name_.c_str(), gnames[i]) == 0)
        {
            ig = i;
            break;
        }
    }

    for (auto i = grps->index[ig]; i < grps->index[ig + 1]; ++i)
    {
        group_.ind_.push_back(grps->a[i]);
    }
    fprintf(stderr, "Density fitting group '%s' has %lu atoms\n", group_.name_.c_str(), group_.ind_.size());
}

int DensfitData::nStFit() const
{return nstfit_; }

void DensfitData::read_densparams(int *ninp_p, t_inpfile **inp_p, warninp *wi)
{

    int          ninp = *ninp_p;
    t_inpfile   *inp  = *inp_p;
    CTYPE("Density Potential Type: 'cross-correlation' or 'cross-entropy' ");
    int          densityPotentialEnumInt_;
    EETYPE("densfit-potential", densityPotentialEnumInt_, densityPotentialNames);
    densityPotential_ = DensityPotential(densityPotentialEnumInt_);

    CTYPE("Index group of atoms used for spreading");
    char         groupCharArray[STRLEN];
    const char * tmp;
    STYPE("spread-group", groupCharArray, "");
    group_.name_ = std::string(groupCharArray);
    char         referenceDensityCharArray[STRLEN];
    STYPE("reference-density", referenceDensityCharArray, "");
    std::string  referenceDensityFn(referenceDensityCharArray);
    map_ref_ = MrcFile().read(referenceDensityFn);
    /* We need at least one value for k and sigma. But if the user wishes
     * to do multiple refinement steps in an automated way, additional
     * values can be added. Interpolation between these values will be
     * done like in the simulated annealing case. */
    CTYPE("List of times (ps) at which density fitting / refinement parameters "
          "are given and exactly matched,");
    CTYPE("interpolate linearly in between");
    char timePointsString[STRLEN];
    STYPE("densfit-time", timePointsString, NULL);

    CTYPE("Gaussian width(s) used to transform the discrete atomic positions "
          "into a density (one sigma per refinement point)");
    char densfitSigmasString[STRLEN];
    STYPE("densfit-sigma", densfitSigmasString, NULL);

    CTYPE("Strength constant(s) k (kJ/mol) of map potential V_fit = k*(1 - "
          "corr.coeff.) (one k per refinement point)");
    char densfitForceConstantString[STRLEN];
    STYPE("densfit-k", densfitForceConstantString, NULL);

    CTYPE("Cutoff distance (in sigmas) for density spreading, spread only within "
          "+/- sigma*dist");
    RTYPE("dist", dist_, 4.0);
    CTYPE("Calculate V_fit every nstfit time steps");
    ITYPE("nstfit", nstfit_, 1);
    CTYPE("Write diagnostic output every nstout time steps");
    ITYPE("nstout", nstout_, 1);
    CTYPE("Write simulated density map to file every nstmapout time steps");
    ITYPE("nstmapout", nstmapout_, 1000);
    CTYPE("Keep and number simulated density maps (results in a bunch of map files)");
    EETYPE("KeepMaps", bKeepAndNumberMaps_, yesno_names);
    CTYPE("Normalize maps to unity.");
    EETYPE("densfit-normalizemapstounity", normalizeMapsToUnity_, yesno_names);

    if (normalizeMapsToUnity_)
    {
        containeroperation::divide(&map_ref_, containermeasure::norm(map_ref_));
    }

    timeDependent_ = TimeDependentDensfitParameters(
                timePointsString, densfitForceConstantString, densfitSigmasString, wi);


    *ninp_p = ninp;
    *inp_p  = inp;
}

DensityPotential DensfitData::densityPotential() const
{
    return densityPotential_;
};
bool DensfitData::normalizeMapsToUnity() const
{
    return normalizeMapsToUnity_;
};

real DensfitData::dist() const { return dist_; };

} // namespace gmx
