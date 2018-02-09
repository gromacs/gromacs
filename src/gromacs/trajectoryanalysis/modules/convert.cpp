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
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Convert.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "convert.h"

#include <algorithm>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/unique_cptr.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * Convert
 */

class Convert : public TrajectoryAnalysisModule
{
    public:
        Convert();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    protected:
        void modifyFrame(t_trxframe *newFrame, const t_trxframe *oldFrame);
        void initOutput();
        void writeFrame(const t_trxframe *coord) const;
        void closeFile();
        void setLegacyInformation(t_atoms *local);
        void trjOpenTng(const char *filename, const char *mode, t_trxstatus *output) const;
        t_trxstatus *trjOpenTrr(const char *filename, const char *mode) const;
        void setConnections();




    private:
        Selection                           sel_;
        std::string                         name_;


        t_trxframe                localFrame_;
        bool                      bVel_;
        bool                      bForce_;
        bool                      bgenCon_;
        bool                      bPrec_;
        bool                      bTime_;
        double                    prec_;
        double                    time_;

        t_trxstatus              *trr_;

        char                      filemode_[5];
        t_trxstatus              *infile_;
        gmx_conect                connections_;
        t_atoms                   atoms_;
        bool                      bAtoms_;
        int filetype_;
        const gmx_mtop_t         *mtop_;

};

Convert::Convert() : sel_(nullptr), bVel_(false), bForce_(false), bgenCon_(false), bPrec_(false), bTime_(false),
                     prec_(0.0), time_(0.0), trr_(nullptr), infile_(nullptr), bAtoms_(false), filetype_(efNR), mtop_(nullptr)
{
    strcpy(filemode_, "w");
}


void
Convert::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] converts trajectory files between different formats."
    };

    settings->setHelpText(desc);

    options->addOption(SelectionOption("Output").store(&sel_).dynamicMask()
                           .required()
                           .description("Selection to write out"));

    options->addOption(FileNameOption("o").filetype(eftTrajectory).outputFile()
                           .store(&name_).defaultBasename("trajout")
                           .required()
                           .description("Output trajectory after conversion"));
    options->addOption(BooleanOption("velocity").store(&bVel_)
                           .description("Write velocities to file if possible"));
    options->addOption(BooleanOption("force").store(&bForce_)
                           .description("Write forces to file if possible"));
    options->addOption(BooleanOption("conect").store(&bgenCon_)
                           .description("Write connection information to file"));
    options->addOption(DoubleOption("prec").store(&prec_)
                           .description("Output precision for compressed files")
                           .storeIsSet(&bPrec_));
    options->addOption(DoubleOption("t0").store(&time_)
                           .description("Set custom start time")
                           .storeIsSet(&bTime_));



    // set correct flags to indicate we need a proper topology for the analysis
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
Convert::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
}


void
Convert::initAnalysis(const TrajectoryAnalysisSettings   & /*settings*/,
                      const TopologyInformation         &top)
{
    mtop_ = top.mtop();
    initOutput();
}

void
Convert::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /* pbc */,
                      TrajectoryAnalysisModuleData * /*pdata*/)
{
    clear_trxframe(&localFrame_, true);

    modifyFrame(&localFrame_, &fr);
    writeFrame(&localFrame_);
}

void
Convert::finishAnalysis(int /*nframes*/)
{
    closeFile();
}



void
Convert::writeOutput()
{
}

void
Convert::modifyFrame(t_trxframe *newFrame, const t_trxframe *oldFrame)
{
    const Selection *sel    = &sel_;
    int              natoms = sel->atomCount();

    *newFrame = *oldFrame;

    newFrame->time   = time_;
    newFrame->bV     = (oldFrame->bV && bVel_);
    newFrame->bF     = (oldFrame->bF && bForce_);
    newFrame->natoms = natoms;
    newFrame->bPrec  = (oldFrame->bPrec && bPrec_);
    newFrame->prec   = prec_;
    newFrame->atoms  = &atoms_;
    newFrame->bAtoms = bAtoms_;

    rvec *xmem = nullptr;
    rvec *vmem = nullptr;
    rvec *fmem = nullptr;
    snew(xmem, natoms);
    if (newFrame->bV)
    {
        snew(vmem, natoms);
    }
    if (newFrame->bF)
    {
        snew(fmem, natoms);
    }
    newFrame->x = xmem;
    newFrame->v = vmem;
    newFrame->f = fmem;


    for (int i = 0; i < natoms; i++)
    {
        int pos = sel->position(i).refId();
        copy_rvec(oldFrame->x[pos], newFrame->x[i]);
        if (newFrame->bV)
        {
            copy_rvec(oldFrame->v[pos], newFrame->v[i]);
        }
        if (newFrame->bF)
        {
            copy_rvec(oldFrame->f[pos], newFrame->f[i]);
        }
    }
}

void
Convert::closeFile()
{
    if (trr_ != nullptr)
    {
        close_trx(trr_);
    }
    trr_ = nullptr;
}

void
Convert::setLegacyInformation(t_atoms *local)
{
    t_atoms inputAtoms = gmx_mtop_global_atoms(mtop_);
    init_t_atoms(local, inputAtoms.nr, inputAtoms.havePdbInfo);
    sfree(local->resinfo);
    local->resinfo = inputAtoms.resinfo;
    int natoms = sel_.atomCount();
    for (int i = 0; (i < natoms); i++)
    {
        local->atomname[i] = inputAtoms.atomname[sel_.position(i).refId()];
        local->atom[i]     = inputAtoms.atom[sel_.position(i).refId()];
        if (inputAtoms.havePdbInfo)
        {
            local->pdbinfo[i] = inputAtoms.pdbinfo[sel_.position(i).refId()];
        }
        local->nres = std::max(local->nres, local->atom[i].resind+1);
    }
    local->nr = natoms;
    bAtoms_   = true;
}

void
Convert::trjOpenTng(const char *filename, const char *mode, t_trxstatus *output) const
{
    ArrayRef<const int> index     = sel_.atomIndices();
    int                 natoms    = sel_.atomCount();
    const char         *indexname = sel_.name();

    int                 size = index.size();
    int                *localindex;
    snew(localindex, size);
    int                 runner = 0;

    for (const int *ai = index.begin(); ai != index.end(); ++ai)
    {
        localindex[runner++] = *ai;
    }

    gmx_mtop_t mtop = (*const_cast<gmx_mtop_t*>(mtop_));

//    t_trxstatus &localInput = impl_->infile_;
    trjtools_gmx_prepare_tng_writing(filename,
                                     mode[0],
                                     nullptr, //infile_, //how to get the input file here?
                                     &output,
                                     nullptr,
                                     natoms,
                                     &mtop,
                                     localindex,
                                     indexname);
}

t_trxstatus *
Convert::trjOpenTrr(const char *filename, const char *mode) const
{
    return open_trx(filename, mode);
}

void
Convert::setConnections()
{
    gmx_mtop_t *mtop = const_cast<gmx_mtop_t*>(mtop_);
    t_topology  top  = gmx_mtop_t_to_t_topology(mtop, false);
    connections_ = gmx_conect(&top);
}


void
Convert::initOutput()
{
    if (!name_.empty())
    {
        filetype_ = fn2ftp(name_.c_str());
        switch (filetype_)
        {
            case (efTNG):
                trjOpenTng(name_.c_str(), filemode_, trr_);
                break;
            case (efTRR):
            case (efXTC):
            case (efPDB):
            case (efGRO):
            case (efG96):
                trr_ = trjOpenTrr(name_.c_str(), filemode_);
                setLegacyInformation(&atoms_);
                break;
            default:
                gmx_incons("Invalid file type");
                // handle this error
        }
    }
    if (bgenCon_ || (filetype_ == efPDB))
    {
        setConnections();
    }
}

/*! \cond libapi */

void
Convert::writeFrame(const t_trxframe *coord) const
{
    const t_trxframe *plocal = coord;
    t_trxframe        local  = (*const_cast<t_trxframe*>(plocal));

    switch (filetype_)
    {
        case (efTNG):
            write_tng_frame(trr_, &local); //local);
            break;
        case (efTRR):
        case (efXTC):
        case (efPDB):
        case (efGRO):
        case (efG96):
            write_trxframe(trr_, &local, connections_);
            break;
        default:
            gmx_incons("Illegal output file format");
    }
}

}       // namespace

const char ConvertInfo::name[]             = "convert";
const char ConvertInfo::shortDescription[] =
    "Converts between different trajectory types";

TrajectoryAnalysisModulePointer ConvertInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Convert);
}

} // namespace analysismodules

} // namespace gmx
