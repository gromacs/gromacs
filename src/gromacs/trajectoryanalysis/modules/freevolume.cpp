/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements gmx::analysismodules::Freevolume.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "freevolume.h"

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/atomprop.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

const char Freevolume::name[]             = "freevolume";
const char Freevolume::shortDescription[] =
    "Calculate free volume";

Freevolume::Freevolume()
    : TrajectoryAnalysisModule(name, shortDescription),
      adata_(new AnalysisDataFrameAverageModule())
{
    data_.setColumnCount(1);
    registerAnalysisDataset(&data_, "freevolume");
}


Freevolume::~Freevolume()
{
    gmx_rng_destroy(rng_);
    gmx_ana_nbsearch_free(nbsearch_);
}


void
Freevolume::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "g_freevol can calculate the free volume in a box as",
        "a function of time. The free volume is",
        "plotted as a fraction of the total volume.",
        "The group specified by the selection is considered to",
        "delineate non-free volume."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                       .store(&fnFreevol_).defaultBasename("freevolume")
                       .description("Computed free volume"));
    options->addOption(SelectionOption("select").required().valueCount(1)
                       .store(sel_));
    probe_radius_ = 0;
    options->addOption(DoubleOption("radius").store(&probe_radius_)
                       .description("Radius of the probe to be inserted (nm, 0 yields the true free volume)"));
    seed_ = -1;
    options->addOption(IntegerOption("seed").store(&seed_)
                       .description("Seed for random number generator."));
    ninsert_ = 1000;
    options->addOption(IntegerOption("ninsert").store(&ninsert_)
                       .description("Number of probe insertions per cubic nm to try for each frame in the trajectory."));
}


void
Freevolume::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation        &top)
{
    int i;
    int nnovdw;
    real value;
    gmx_atomprop_t aps;
    t_atoms *atoms;
    
    if (!top.hasTopology())
    {
        GMX_THROW(InvalidInputError("Oops, no topology"));
    }
    
    if (sel_[0].posCount() < 1)
    {
        GMX_THROW(InvalidInputError("The selection does not contain atoms"));
    }

    //data_.addModule(avem_);
    data_.addModule(adata_);
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(fnFreevol_);
    plotm_->setTitle("Free Volume");
    plotm_->setXAxisIsTime();
    plotm_->setYLabel("Free Volume (\%)");
    data_.addModule(plotm_);
    
    cutoff_ = 0;
    nnovdw  = 0;
    aps     = gmx_atomprop_init();
    atoms   = &(top.topology()->atoms);
    
    for(i=0; (i<atoms->nr); i++)
    {
        int resnr = atoms->atom[i].resind;
        if (TRUE == gmx_atomprop_query(aps,epropVDW,
                                       *(atoms->resinfo[resnr].name),
                                       *(atoms->atomname[i]),
                                       &value))
        {
            vdw_radius_.push_back(value);
            if (value > cutoff_)
            {
                cutoff_ = value;
            }
        }
        else
        {
            nnovdw++;
            if (nnovdw < 10)
            {
                fprintf(stderr,"Could not determine VDW radius for %s-%s. Set to zero.\n",
                        *(atoms->resinfo[resnr].name),
                        *(atoms->atomname[i]));
                }               
            vdw_radius_.push_back(0.0);
        }
    }
    gmx_atomprop_destroy(aps);
    cutoff_ += probe_radius_;
    if (nnovdw >= 10)
    {
        fprintf(stderr,"Could not determine VDW radius for %d particles. These were set to zero.\n",nnovdw);
    }

    printf("cutoff       = %g nm\n",cutoff_);
    printf("probe_radius = %g nm\n",probe_radius_);
    printf("seed         = %d\n",seed_);
    printf("ninsert      = %d probes per nm^3\n",ninsert_);
    rng_ = gmx_rng_init(seed_);
    nbsearch_ = gmx_ana_nbsearch_create(cutoff_,atoms->nr);
    nframes_ = 0;
    fractot_ = 0;
    frac2tot_ = 0;
}

static gmx_inline real det_c(const matrix a)
{
    return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
             -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
             +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}

static gmx_inline void mvmul_c(const matrix a, const rvec src, rvec dest)
{
    dest[XX] = a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
    dest[YY] = a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
    dest[ZZ] = a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}

void
Freevolume::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                         TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle       dh   = pdata->dataHandle(data_);
    const Selection         &sel1 = pdata->parallelSelection(sel_[0]);
    const SelectionPosition &p1   = sel1.position(0);
    rvec                     rand,ins,dx;
    double                   frac;
    real                     r,V,rins;
    int                      i,m,jp,Ninsert,NinsTot;
    bool                     bOverlap;
    gmx_ana_pos_t            allpos,pos;
    gmx_ana_index_t          g;
    
    dh.startFrame(frnr, fr.time);
    if (pbc == NULL)
    {
        exit(1);
    }
    // Compute volume and number of insertions to perform
    V = det_c(fr.box);
    Ninsert = ninsert_*V;
    
    // Use neighborsearching tools! This seems overly complex.
    gmx_ana_pos_clear(&pos);
    gmx_ana_pos_clear(&allpos);
    gmx_ana_index_clear(&g);
    gmx_ana_index_init_simple(&g,fr.natoms,strdup("Hoi"));
    gmx_ana_pos_reserve(&allpos,fr.natoms,fr.natoms);
    gmx_ana_pos_set_nr(&allpos,0);
    for(i=0; (i<fr.natoms); i++)
    {
        gmx_ana_pos_init_const(&pos,fr.x[i]);
        gmx_ana_pos_append_init(&allpos,&g,&pos,0);
        gmx_ana_pos_deinit(&pos);
    }
    gmx_ana_pos_append_finish(&allpos);
    gmx_ana_nbsearch_pos_init(nbsearch_,pbc,&allpos);
    gmx_ana_index_deinit(&g);
    
    // Then loop over insertions
    NinsTot = 0;
    for(i=0; (i<Ninsert); i++)
    {
        for(m=0; (m<DIM); m++)
        {
            // Generate random number between 0 and 1
            rand[m] = gmx_rng_uniform_real(rng_);
        }
        // Generate random 3D position within the box
        mvmul_c(fr.box,rand,ins);
        
        // Find the first reference position within the cutoff.
        bOverlap = false;
        if (true == gmx_ana_nbsearch_first_within(nbsearch_,ins,&jp))
        {
            pbc_dx(pbc,ins,fr.x[jp],dx);
            bOverlap = (norm(dx) < vdw_radius_[jp]);

            /** Finds the next reference position within the cutoff. */
            while (!bOverlap && 
                   (true == gmx_ana_nbsearch_next_within(nbsearch_,&jp)))
            {
                pbc_dx(pbc,ins,fr.x[jp],dx);
                bOverlap = (norm(dx) < vdw_radius_[jp]);
            }
        }
        
        if (!bOverlap)
        {
            NinsTot++;
        }
    }
    if (Ninsert > 0)
    {
        frac = (1.0*NinsTot)/Ninsert;
        fractot_  += frac;
        frac2tot_ += frac*frac;
        dh.setPoint(0,frac);
        nframes_++;
    }
    gmx_ana_pos_deinit(&allpos);
    dh.finishFrame();
}


void
Freevolume::finishAnalysis(int nframes)
{
}


void
Freevolume::writeOutput()
{
    double faver  = fractot_/nframes_;
    double fsigma = sqrt(frac2tot_/nframes_-faver*faver);
    printf("Free volume %.2f +/- %.2f %%\n",100*faver,100*fsigma);
    printf("Number of frames %d\n",data_.frameCount());
}

} // namespace analysismodules

} // namespace gmx
