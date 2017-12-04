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

#include "densfitdata.h"

#include <cassert>
#include <cstdio>
#include <vector>

#include "gromacs/utility/iserializer.h"

#include "gromacs/fileio/griddataio.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/inmemoryserializer.h"

#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdlib/sim_util.h"

#include "gromacs/mdlib/groupcoord.h"

#include "gromacs/domdec/domdec_struct.h"

#include "gromacs/timing/cyclecounter.h"

#include "gromacs/domdec/domdec.h"

#include "gromacs/fileio/griddataio.h"

#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/pbcutil/pbc.h"

#include "gromacs/fileio/mrcmetadata.h"

#include "gromacs/math/vecdump.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/txtdump.h"

#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"

static const char *FitStr = ""; //!< Used in gmx map output
//! \brief Environment variable for setting nstfit
static const char *NSTFIT_ENVVAR = "GMX_NSTFIT";

#define MAXPTR 254

namespace gmx
{

namespace
{
//! \brief Adds 'buf' to 'str'
static void add_to_string(char **str, char *buf)
{
    int len;

    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}

//! \brief TODO
static void add_to_string_aligned(char **str, const char *buf)
{
    char buf_aligned[STRLEN];

    sprintf(buf_aligned, "%12s", buf);

    add_to_string(str, buf_aligned);
}
}

TimeDependentDensfitParameters::TimeDependentDensfitParameters(
        const std::vector<real> &times, const std::vector<real> &forceConstants,
        const std::vector<real> &sigmas)
    : timePoints_(times), forceConstants_(forceConstants), sigmas_(sigmas) {}

bool TimeDependentDensfitParameters::vary() const
{
    return timePoints_.size() > 1;
}

TimeDependentDensfitParameters::TimeDependentDensfitParameters(
        const char timePointString[STRLEN], const char forceConstantString[STRLEN],
        const char sigmaString[STRLEN], warninp *wi)
{

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
            warning(wi, "Got a negative value for a density fitting k value. This "
                    "will expel atoms from the map!");
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
    fprintf(fp, "#   Note that the force update frequency can be overwritten "
            "by the environment variable %s.\n",
            NSTFIT_ENVVAR);
    fprintf(fp, "# Output is written in intervals of %d time step%s.\n", nstout_,
            nstout_ > 1 ? "s" : "");
    fprintf(fp, "# Parameter sigma_dist = %g\n", dist_);
    if (0 != nstmapout_)
    {
        fprintf(fp, "# Calculated map will be saved every %d steps\n", nstmapout_);
    }
    int          nsets;
    const char **setname;
    char         buf[50];
    char        *LegendStr = NULL;
    /* Print a nice legend */
    snew(LegendStr, 1);
    LegendStr[0] = '\0';
    sprintf(buf, "#     %6s", "time");
    add_to_string_aligned(&LegendStr, buf);

    nsets = 0;
    snew(setname, 5);
    setname[nsets++] = strdup("temperature (K)");
    add_to_string_aligned(&LegendStr, "T");
    setname[nsets++] = strdup("k (kJ/mol)");
    add_to_string_aligned(&LegendStr, "k");
    setname[nsets++] = strdup("sigma (nm)");
    add_to_string_aligned(&LegendStr, "sigma");
    setname[nsets++] = strdup("correlation coefficient");
    add_to_string_aligned(&LegendStr, "c.c.");
    setname[nsets++] = strdup("V_dens (kJ/mol)");
    add_to_string_aligned(&LegendStr, "V_dens");
    xvgr_legend(fp, nsets, setname, oenv);
    sfree(setname);

    fprintf(fp, "#\n# Legend for the following data columns:\n");
    fprintf(fp, "%s\n", LegendStr);
    sfree(LegendStr);
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

    serializer->doBool(&bKeepAndNumberMaps_);

}

DensfitData::DensfitData(const DensfitData &other)
{
    map_ref_            = other.map_ref_;
    timeDependent_      = other.timeDependent_;
    dist_               = other.dist_;
    nstfit_             = other.nstfit_;
    nstout_             = other.nstout_;
    nstmapout_          = other.nstmapout_;
    bKeepAndNumberMaps_ = other.bKeepAndNumberMaps_;
    group_              = other.group_;
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
    pr_str(fp, indent, "densfit-bKeepMaps", boolToString(bKeepAndNumberMaps_));
    pr_ivec_block(fp, indent, "atom", group_.ind_.data(), group_.ind_.size(), true);

}

DensfitData::DensfitData(const std::vector<real> &sigma,
                         const std::vector<real> &k, real sigma_dist,
                         int nAtoms, int *index, bool bVerbose)
{

    dist_        = sigma_dist;

    TimeDependentDensfitParameters({0.}, sigma, k);
    if (index == nullptr)
    {
        if (bVerbose)
        {
            fprintf(stdout, "%sYou did not provide an index group - will use all "
                    "atoms for spreading.\n",
                    FitStr);
        }

        group_.ind_.resize(nAtoms);
        std::iota(std::begin(group_.ind_), std::end(group_.ind_), 0);
    }
    else
    {
        group_.ind_.assign(index, index+nAtoms);
    }
}

const GridDataReal3D &DensfitData::referenceMap() const
{
    return map_ref_;
}

void DensfitData::makeGroups(t_blocka *grps, char **gnames)
{
    /* Group index from group name */
    int ig = search_string(group_.name_.c_str(), grps->nr, gnames);

    fprintf(stderr, "Density fitting group '%s' has %lu atoms\n", group_.name_.c_str(), group_.ind_.size());
    for (int i = grps->index[ig]; i < grps->index[ig + 1]; i++)
    {
        group_.ind_.push_back(grps->a[i]);
    }
}

void DensfitData::read_densparams(int *ninp_p, t_inpfile **inp_p, warninp *wi)
{

    int          ninp = *ninp_p;
    t_inpfile   *inp  = *inp_p;
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
    CTYPE("Keep and number simulated density maps (results in a bunch of map "
          "files)");
    EETYPE("KeepMaps", bKeepAndNumberMaps_, yesno_names);

    timeDependent_ = TimeDependentDensfitParameters(
                timePointsString, densfitForceConstantString, densfitSigmasString, wi);

    *ninp_p = ninp;
    *inp_p  = inp;
}

real DensfitData::dist() const { return dist_; };

} // namespace gmx
