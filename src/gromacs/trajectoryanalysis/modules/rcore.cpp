/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017 by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::RCore.
 *
 * \author Anatoly Titov <toluk@omrb.pnpi.spb.ry>
 * \author Alexey Shvetsov <alexxy@omrb.pnpi.spb.ru>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "rcore.h"

#include <cmath>

#include <vector>
#include <cstdlib>
#include <cstdio>

#include "gromacs/fileio/trxio.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

namespace analysismodules
{

namespace
{


/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class RCore : public TrajectoryAnalysisModule
{
    public:
        RCore();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        Selection                                           sel_;
        std::string                                         fnNdx_;
        std::vector< std::vector < RVec > >                 trajectory_;
        std::vector< double >                               rmsd_;
        std::vector< int >                                  ndx_;
        double                                              rmsdRef_;
        int                                                 frRef_;
};


RCore::RCore()
    : rmsdRef_(0.1),
      frRef_(0)
{
}


void
RCore::initOptions(IOptionsContainer          *options,
                   TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates rigid core of macromolecule,",
        "and write index file with atoms corresponding to it.",
        "The group in resulting index file can be used as reference",
        "group for other tools e.g. gmx rms."
        "",
        "Option [TT]-select[tt] sets subset of atoms that will be uses",
        "to find rigid core from",
        "",
        "Option [TT]-on[tt] sets name for output index file",
        "",
        "Option [TT]-rmsdref[tt] sets how rigid core should be:",
        "only atoms which rmsd is lower then [TT]-rmsdref[tt] will be",
        "selected as rigid ones",
        "",
        "Option [TT]-frref[tt] used to set against which frame RMSD",
        "should be calculated"
    };

    // Add the descriptive text (program help text) to the options
    settings->setHelpText(desc);
    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select")
                           .store(&sel_).required()
                           .description("Atoms that are considered as part of the excluded volume"));
    // Add option for rcore index file name
    options->addOption(FileNameOption("on").filetype(eftIndex)
                           .outputFile().required()
                           .store(&fnNdx_).defaultBasename("rcore")
                           .description("Index file from the rcore"));
    // Add option for rmsdRef constant
    options->addOption(DoubleOption("rmsdref")
                           .store(&rmsdRef_)
                           .description("Minimal RMSD for RCore"));
    // Add option for frRef constant
    options->addOption(gmx::IntegerOption("frref")
                           .store(&frRef_)
                           .description("Reference frame"));
}

void
RCore::initAnalysis(const TrajectoryAnalysisSettings  & /*settings*/,
                    const TopologyInformation         & /*top*/)
{
    ndx_.resize(0);
    ConstArrayRef< int > atomind = sel_.atomIndices();
    for (ConstArrayRef< int >::iterator ai = atomind.begin(); (ai < atomind.end()); ++ai)
    {
        ndx_.push_back(*ai);
    }
}


void
RCore::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/,
                    TrajectoryAnalysisModuleData * /*pdata*/)
{
    trajectory_.resize(frnr+1);
    trajectory_[frnr].resize(ndx_.size());
    for (size_t i = 0; i < ndx_.size(); i++)
    {
        trajectory_[frnr][i] = fr.x[ndx_[i]];
    }
}


void
RCore::finishAnalysis(int nframes)
{
    real               *w_rls;
    rvec               *x1, *x2, temp;
    std::vector< bool > flag;
    double              noise_mid = 9999;
    int                 count     = 0;
    FILE               *fpNdx_;

    snew(x1, ndx_.size());
    snew(x2, ndx_.size());
    snew(w_rls, ndx_.size());
    flag.resize(ndx_.size(), true);
    rmsd_.resize(ndx_.size(), 0);

    while (noise_mid > rmsdRef_)
    {
        noise_mid = 0;
        count     = 0;
        rmsd_.assign(ndx_.size(), 0);

        for (size_t i = 0; i < ndx_.size(); i++)
        {
            if (flag[i])
            {
                w_rls[i] = 1;
            }
            else
            {
                w_rls[i] = 0;
            }
        }

        for (int i = 0; i < nframes; i++)
        {
            for (size_t j = 0; j < trajectory_[i].size(); j++)
            {
                copy_rvec(trajectory_[i][j].as_vec(), x2[j]);
                copy_rvec(trajectory_[frRef_][j].as_vec(), x1[j]); // i think its nessesary to do because of lots of reset_x
            }
            reset_x(ndx_.size(), NULL, ndx_.size(), NULL, x1, w_rls);
            reset_x(ndx_.size(), NULL, ndx_.size(), NULL, x2, w_rls);
            do_fit(ndx_.size(), w_rls, x1, x2);

            for (size_t j = 0; j < ndx_.size(); j++)
            {
                rvec_sub(x1[j], x2[j], temp);
                rmsd_[j] += norm(temp);
            }
        }

        for (size_t i = 0; i < ndx_.size(); i++)
        {
            rmsd_[i] /= nframes;
        }

        for (size_t i = 0; i < ndx_.size(); i++)
        {
            if (flag[i])
            {
                count++;
                noise_mid += rmsd_[i];
            }
        }
        noise_mid /= count;

        for (size_t i = 0; i < ndx_.size(); i++)
        {
            if (rmsd_[i] > noise_mid)
            {
                flag[i] = false;
            }
        }
    }

    int write_count = 0;

    fpNdx_ = std::fopen(fnNdx_.c_str(), "w+");
    std::fprintf(fpNdx_, "[ rcore ]\n");
    for (size_t i = 0; i < rmsd_.size(); i++)
    {
        if (rmsd_[i] <= rmsdRef_)
        {
            write_count++;
            if (write_count > 20)
            {
                write_count -= 20;
                std::fprintf(fpNdx_, "\n");
            }
            std::fprintf(fpNdx_, "%5d ", ndx_[i] + 1);
        }
    }
    std::fprintf(fpNdx_, "\n");
    std::fclose(fpNdx_);
    sfree(x1);
    sfree(x2);
    sfree(w_rls);
}


void
RCore::writeOutput()
{
}

}       // namespace

const char RCoreInfo::name[]             = "rcore";
const char RCoreInfo::shortDescription[] =
    "Calculate Rigid Core of Macromolecule";

TrajectoryAnalysisModulePointer RCoreInfo::create()
{
    return TrajectoryAnalysisModulePointer(new RCore);
}

} // namespace analysismodules

} // namespace gmx
