/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares trajectory analysis module for RNA analysis (PDB structures and trajectories).
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders Gärdenäs <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_RNAANALYSIS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_RNAANALYSIS_H

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <list>
#include <sstream>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "basepairdb.h"

namespace gmx
{

namespace analysismodules
{

class RnaAnalysisInfo
{
    public:
        static const char name[];
        static const char shortDescription[];
        static TrajectoryAnalysisModulePointer create();
};

//! The type of iterator of pointer to BasePair
typedef std::vector<BasePair *>::iterator BasePairPtrIterator;

//! \internal \brief
//! Info about a base pair, base1 and base2 is the index of ResInfo
class PairInfo
{
    public:
        //! Constructor
        PairInfo(size_t              base1,
                 size_t              base2,
                 real                score,
                 real                RMSD,
                 BasePairPtrIterator bpi) : base1_(base1), base2_(base2),
                                            score_(score), RMSD_(RMSD),
                                            swap_(false), bpi_(bpi) { }

        //! Base index 1 in the structure
        size_t base1() const { return base1_; }

        //! Base index 2 in the structure
        size_t base2() const { return base2_; }

        //! The best score
        real score() const { return score_; }

        //! The best RMSD
        real RMSD() const { return RMSD_; }

        //! Whether should swap
        bool swap() const { return swap_; }

        //! Set the value of swap
        void setSwap(bool swap) { swap_ = swap; }

        //! An iterator of the pointer to the template base pair.
        BasePairPtrIterator bpi() const { return bpi_; }

        //! Set the value of score
        void setScore(real score) { score_ = score; }

        //! Set the value of RMSD
        void setRMSD(real RMSD) { RMSD_ = RMSD; }

        //! Set the pointer to entrance in the template data base
        void setBPI(BasePairPtrIterator bpi) { bpi_ = bpi; }

    private:
        size_t              base1_;
        size_t              base2_;
        real                score_;
        real                RMSD_;
        bool                swap_;
        BasePairPtrIterator bpi_;
};

//! \internal \brief
//! Holds the frame coordinates
struct frame
{
    //! The coordinates
    rvec * vec;
};

//! \internal \brief
//! Class used to analyze RNA structures and trajectories.
//! Inherits TrajectoryAnalysisModule and all functions from there.
//! \ingroup module_trajectoryanalysis
class RnaAnalysis : public TrajectoryAnalysisModule
{
    public:

        //! Constructor
        RnaAnalysis();

        //! Destructor
        virtual ~RnaAnalysis();

        //! Set the options and settings
        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);

        //! First routine called by the analysis frame work
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        //! Call for each frame of the trajectory
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        //! Last routine called by the analysis frame work
        virtual void finishAnalysis(int /* nframes*/);

        //! Routine to write output, that is additional over the built-in
        virtual void writeOutput();

    private:
        //! Check if it is a valid atom to save and get the atom type
        bool isValid(const char * name);

        size_t findSameAtom(const ResInfo &base, const char * name);

        //! Read the template data base
        void readTemplates();

        //! Read the topology
        void readTopology(const TopologyInformation &top);

        void printBasePair(const PairInfo &tempPair);

        real analyzeBasePair(BasePairPtrIterator bpi,
                             const ResInfo      &base1,
                             const ResInfo      &base2,
                             t_pbc              *pbc,
                             const Selection    &sel,
                             real               *RMSD,
                             rvec               *vec,
                             size_t              vecSize);
        Selection                        sel_;
        AnalysisData                     data_;
        AnalysisDataAverageModulePointer adata_;
        AnalysisNeighborhood             nb_;

        gmx_atomprop_t                   atomprop_;
        // Copy and assign disallowed by base.
        // A list list of a found pairs
        std::list < std::list < PairInfo>>   found;
        std::list<frame>                 coords;
        // Our storage of residue information
        std::vector<ResInfo *>           nucleotideInfo_;

        std::vector<BasePair *>          templateDB_;
        std::string                      DB;
        std::string                      outPath;
        std::string                      rnaAnalysis;
        std::string                      rnaLogFile;
        FILE                            *rnaLog;

        // Boolean options
        bool                             hydrogenRmsd_;
        bool                             sugarRmsd_;
        bool                             phosphateRmsd_;
        bool                             printAllBasePairs_;
        bool                             statistic;
        bool                             oneBPList;
        bool                             Xpm;

        double                           extraOffset;
        double                           addRMSD_;
        double                           bondDist_;
        bool                             datailedInfo_;
        t_atoms                         *iAtoms;
        int                              offsetAtom;
};
} // namespace analysismodules

} // namespace gmx

#endif
