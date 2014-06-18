/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::RnaAnalysis.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders Gärdenäs <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "rnaAnalysis.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <list>
#include <vector>

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#ifdef GMX_NATIVE_WINDOWS
#include <direct.h>
#include <io.h>
#endif

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/uniqueptr.h"

#include "basepairdb.h"

namespace gmx
{

namespace analysismodules
{

/*! \brief
 * Class used to analyze RNA structures and trajectories.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 *
 * \ingroup module_trajectoryanalysis
 */

//! Info about a base pair, base1 and base2 is the index of resInfo
struct pairInfo
{
    //! Base index 1
    int base1;
    //! Base index 2
    int base2;
    //! Index for the pair
    int pairIndex;
};

// Holds the frame coordinates
struct frame
{
    rvec * vec;
};

// Definition of class RnaAnalysis
class RnaAnalysis : public TrajectoryAnalysisModule
{
    public:

        //! Constructor
        RnaAnalysis();

        //! Destructor
        virtual ~RnaAnalysis();

        // Set the options and setting
        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);

        // First routine called by the analysis frame work
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        // Call for each frame of the trajectory
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        // Last routine called by the analysis frame work
        virtual void finishAnalysis(int /* nframes*/);

        // Routine to write output, that is additional over the built-in
        virtual void writeOutput();

    private:
        // Check if it is a valid atom to save and get the atom type
        bool getAtomType(const char * name, char * type);

        int  findSameAtom(int i, const char * name);

        //! Read the template data base
        void readTemplates();

        //! Read the topology
        void readTopology(const TopologyInformation &top);

        Selection                           sel_;
        AnalysisData                        data_;
        AnalysisDataAverageModulePointer    adata_;
        AnalysisNeighborhood                nb_;

        gmx_atomprop_t                      atomprop;
        // Copy and assign disallowed by base.
        // A list list of a found pairs
        std::list < std::list < pairInfo> > found;
        std::list <frame>                   coords;
        //! Our storage of residue information
        std::vector<resInfo *>              nucleotideInfo_;

        std::vector<BasePair *>             dataBase_;
        std::string                         DB;
        std::string                         outPath;
        std::string                         rnaAnalysis;
        std::string                         rnaLogFile;
        FILE                               *rnaLog;

        // Boolean options
        bool                               hydrogenRmsd_;
        bool                               sugarRmsd_;
        bool                               phosphateRmsd_;
        bool                               everyBPList;
        bool                               statistic;
        bool                               oneBPList;
        bool                               Xpm;

        double                             extraOffset;
        double                             maxRMSD;
        double                             bondDist;
        bool                               Info;
        t_atoms                          * iAtoms;
        int                                offsetAtom;
};

// Constructor
RnaAnalysis::RnaAnalysis()
    : TrajectoryAnalysisModule(RnaAnalysisInfo::name, RnaAnalysisInfo::shortDescription),
      adata_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(0, 1);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "rnaAnalysis");

    offsetAtom         = 0;
    hydrogenRmsd_      = false;
    sugarRmsd_         = false;
    phosphateRmsd_     = false;
    everyBPList        = false;
    statistic          = false;
    oneBPList          = true;
    Xpm                = false;
    extraOffset        = 0.0;
    maxRMSD            = 0.0;
    Info               = false;
    iAtoms             = NULL;
    bondDist           = 0.0;
    atomprop           = gmx_atomprop_init();
    rnaLog             = NULL;
}

// Destructor
RnaAnalysis::~RnaAnalysis()
{
    // Destroy C structures where there is no automatic memory release
    // C++ takes care of memory in classes (hopefully)
    gmx_atomprop_destroy(atomprop);
}

void
RnaAnalysis::initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Please note that we are still working on improving this tool!\n",
        "Reads and analyzes an RNA structure (-s [<.gro/..>]) or an MD simulation trajectory (-f [<.xtc/.trr/..>] -s [<.gro/..>])",
        "to detect all bases that are in contact (e.g., form base pairs).\n\n",
        "This is can be done in two ways:\n",
        "1. Check if hydrogen bonds between two RNA bases in the structure match a certain distance criteria (-templateDB hbonds).\n",
        "2. compare RNA structure bases to a certain base pair template database (-templateDB bps).\n\n",
        "The data is displayed on the terminal as raw data. More information about",
        "hydrogen bond distances and RMSD to template structures can be displayed(-moreInfo true).\n",
        "The data can also printed to an xpm file can be coverted to an eps file displaying a matrix",
        "describing all base pairs (-matrix true).\n",
        "If you want, you can also store all detected base pairs to pdb files (-pdbFile yes),",
        "you can also specify the directory in which the pdb files are stored (-outPath <string>)."

    };

    // Add the descriptive text (program help text) to the options
    settings->setHelpText(desc);

    // Add option for optional output file
    options->addOption(FileNameOption("g").filetype(eftGenericData).outputFile()
                           .store(&rnaLogFile).defaultBasename("result")
                           .description("log file" ).required());

    // Add option for optional output file
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&rnaAnalysis).defaultBasename("result")
                           .description("output file" ));

    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select").valueCount(1)
                           .store(&sel_).defaultSelectionText("all").onlyAtoms()
                           .description("Use default (all); select which atoms should be used for analysis"));

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);

    // Template databases in which base pairs are stored
    options->addOption(StringOption("templateDB").store(&DB).defaultValue("bps")
                           .description("define the template DB you want to use (currently just bps available)"));

    // Print all results about RNA base pairs during analysis
    options->addOption(BooleanOption("bpLists").store(&everyBPList).defaultValue(false)
                           .description("print all possible base pairs during analysis"));

    // Print the result about RNA base pairs in the end
    options->addOption(BooleanOption("bpListFinal").store(&oneBPList).defaultValue(true)
                           .description("print all identified base pairs after analysis"));

    // Include hydrogen atoms in RMSD calculation
    options->addOption(BooleanOption("hydro").store(&hydrogenRmsd_).defaultValue(true)
                           .description("include the hydrogen atoms of the nucleotide in RMSD calculations"));

    // Include  phosphate group atoms in RMSD calculation
    options->addOption(BooleanOption("phos").store(&phosphateRmsd_).defaultValue(true)
                           .description("include the phosphate atoms of the nucleotide in RMSD calculation"));

    // Include ONLY base atoms of nucleotides in RMSD calculations
    options->addOption(BooleanOption("sugar").store(&sugarRmsd_).defaultValue(true)
                           .description("include the sugar atoms of the nucleotide in RMSD calculation"));

    // Distance of hydrogen bonds
    options->addOption(DoubleOption("bondDist").store(&bondDist).defaultValue(0.3)
                           .description("distance between hydrogen bond donor and acceptor atom"));

    // Maximum RMSD for between identified base pair and template base pair
    options->addOption(DoubleOption("addRMSD").store(&maxRMSD).defaultValue(0.0)
                           .description("increase the maximum RMSD (cut-off value) by offset for identifying base pairs (default = 0.0)"));

    // More information about output data
    options->addOption(BooleanOption("moreInfo").store(&Info).defaultValue(false)
                           .description("details of H-bond distances and RMSD values between identified and template base pairs"));


    // Extra offset used when looking for base pairs in close proximity
    options->addOption(DoubleOption("addOffset").store(&extraOffset).defaultValue(0.5)
                           .description("increase distance by offset when searching for atoms within a certain range"));

    // Matrix output file
    options->addOption(BooleanOption("matrix").store(&Xpm).defaultValue(false)
                           .description("matrix xpm output file is printed, one file for each frame of the trajectory (so far only for templateDB bps)"));

    // Saves coordinates of base pairs in PDB files
    options->addOption(BooleanOption("pdbFiles").store(&statistic).defaultValue(false)
                           .description("save coordinates of all identified RNA base pairs in PDB files"));

    // Output path to directory in which PDB files are stored
    options->addOption(StringOption("outPath").store(&outPath).defaultValue("")
                           .description("output path to directory which already exists (default current directory)"));

}

void RnaAnalysis::readTemplates()
{
    // Store all file names of RNA base pair templates within input directory
    // in vector pdbTemplates
    std::string dbName = "RNA-" + DB + ".pdb";
    std::string str1   = "amber99sb-ildn.ff/" + dbName;

    // Set path and load file pointer
    real      maxDist = 0;

    // Initialize RNA template data base to compare given RNA structure to
    BasePair   *bp = NULL;
    TextReader  fp(gmxlibfn(str1.c_str()));
    std::string line;
    while (fp.readLine(&line))
    {
        if (NULL == bp)
        {
            bp = new BasePair(hydrogenRmsd_, sugarRmsd_, phosphateRmsd_);
        }
        std::vector<std::string> elements = splitString(line);
        if  (elements.size() == 0)
        {
            continue;
        }
        int         nElem    = 0;
        std::string pdbType  = elements[nElem++];

        if (pdbType.compare("MODEL") == 0)
        {
            continue;
        }
        else if ((pdbType.compare("REMARK") == 0) && (elements.size() >= 3))
        {
            bp->setNucleotides(elements[nElem++]);
            if ((DB.compare("bps") == 0) && (elements.size() >= 5))
            {
                bp->setBondType(elements[nElem++]);
                bp->setIsomerism(elements[nElem++]);
                bp->setTemplateRMSD(std::strtod(elements[nElem++].c_str(), NULL));
            }
        }
        else if ((pdbType.compare("CRYST1") == 0) && (elements.size() >= 7))
        {
            matrix box;
            int    ePBC;
            read_cryst1(line.c_str(), &ePBC, box);
            bp->setPBC(ePBC, box);
        }
        else if ( (pdbType.compare("ATOM") == 0) &&
                  (elements.size() >= 8) )
        {
            int         anum  = std::atoi(elements[nElem++].c_str());
            std::string aname = elements[nElem++];
            std::string rname = elements[nElem++];
            nElem++; // Residue number
            rvec        x;
            x[XX]            = std::strtod(elements[nElem++].c_str(), NULL);
            x[YY]            = std::strtod(elements[nElem++].c_str(), NULL);
            x[ZZ]            = std::strtod(elements[nElem++].c_str(), NULL);
            real        m;
            if (gmx_atomprop_query(atomprop, epropMass, "???",
                                   aname.c_str(), &m))
            {
                bp->addAtom(x, anum, aname, m);
            }
            else
            {
                GMX_THROW(InvalidInputError(line));
            }

        }
        else if ((pdbType.compare("CONECT") == 0) &&
                 (elements.size() >= 3))
        {
            int ai = std::atoi(elements[nElem++].c_str());
            int aj = std::atoi(elements[nElem++].c_str());
            bp->addBondAtom(ai, aj);
        }
        else if (pdbType.compare("ENDMDL") == 0)
        {
            // Assign the center of mass to origin of all atoms
            // in three dimensions.
            reset_x_ndim(3, bp->nrAtoms(),  NULL,
                         bp->nrAtoms(), NULL,
                         bp->atomPosition(),
                         bp->atomMass());

            // Set the distance between atoms
            bp->setAtomDist(offsetAtom);
            if (bp->maximumDistance() > maxDist)
            {
                maxDist = bp->maximumDistance();
            }
            dataBase_.push_back(bp);
            bp           = NULL;
        }
        else
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "Unknown pdbType %s", pdbType.c_str());
            GMX_THROW(InvalidInputError(buf));
        }
    }
    fp.close();

    fprintf(rnaLog, "There are %u RNA base pair types in template database %s.\n",
            static_cast<unsigned int>(dataBase_.size()),
            dbName.c_str());
    fprintf(rnaLog, "Will use %g nm for neighborsearching\n",
            maxDist+extraOffset);

    // Initiate the neighborsearching code
    nb_.setCutoff(maxDist + extraOffset);
    // TODO: Did not work before, but we don't know why....
    nb_.setMode(AnalysisNeighborhood::eSearchMode_Grid);
}

void RnaAnalysis::readTopology(const TopologyInformation &top)
{
    // Initialize variables
    iAtoms = &(top.topology()->atoms);
    std::vector<int>     rnaRes;
    char                 type;
    const char          *grpnames;

    // Set up the residue type
    gmx_residuetype_t *rt;
    gmx_residuetype_init(&rt);
    assert(rt);

    // Get all the RNA residues
    for (int i = 0; i < iAtoms->nres; i++)
    {
        gmx_residuetype_get_type(rt, *iAtoms->resinfo[i].name, &grpnames);
        if (strcmp("RNA", grpnames) == 0)
        {
            rnaRes.push_back(i);
        }
    }

    // Set all the values of the nucleotide
    int       cAtom = 0;
    bool      taken = false;

    for (unsigned int i = 0; i < rnaRes.size(); i++)
    {
        resInfo *ri = NULL;

        for (; cAtom < iAtoms->nr; cAtom++)
        {
            if (taken && (rnaRes[i] != iAtoms->atom[cAtom].resind))
            {
                // New nucleotide
                taken = false;
                break;
            }
            // Test if the current atom is within the right nucleotide
            if (rnaRes[i] == iAtoms->atom[cAtom].resind)
            {
                if (!taken)
                {
                    taken = true;
                    ri    = new resInfo(*iAtoms->resinfo[rnaRes[i]].name[0],
                                        iAtoms->resinfo[rnaRes[i]].nr,
                                        cAtom);
                }
                // Add an atom.
                ri->addAtom(*iAtoms->atomname[cAtom],
                            !getAtomType(*iAtoms->atomname[cAtom], &type));
                cAtom++;
            }
        }
        if (NULL != ri)
        {
            nucleotideInfo_.push_back(ri);
        }
    }

    fprintf(rnaLog, "There are %u RNA bases in input file.\n",
            static_cast<unsigned int>(rnaRes.size()));

    gmx_residuetype_destroy(rt);
}

void
RnaAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top)
{
    // Open log file
    rnaLog = gmx_fio_fopen(rnaLogFile.c_str(), "w");

    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);
    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data etc.
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(rnaAnalysis);
    plotm_->setTitle("RNA Analysis");
    plotm_->setXAxisIsTime();
    plotm_->setYLabel("RNA Analysis (%)");
    plotm_->appendLegend("RNA Analysis");
    data_.addModule(plotm_);

    readTopology(top);

    readTemplates();
}

static bool validateResInfo(resInfo *ri)
{
    return (((ri->residueType() == 'A') && (ri->nAtoms() == 30)) ||
            ((ri->residueType() == 'U') && (ri->nAtoms() == 28)) ||
            ((ri->residueType() == 'G') && (ri->nAtoms() == 31)) ||
            ((ri->residueType() == 'C') && (ri->nAtoms() == 29)));
}

void
RnaAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata)
{
    // When a base pair is found
    real                          RMSD           = 0;
    bool                          dimReset       = false;
    bool                          bRMSD          = true;
    int                           baseSize       = 0;

    std::list<pairInfo>           pairs;
    pairInfo                      tempPair   = {0, 0, 0};
    int                           adjustment = 0;
    rvec                          vec[200];
    rvec                          vec_swap[200];
    unsigned int                  base2 = 0;
    AnalysisDataHandle            dh    = pdata->dataHandle(data_);
    const Selection              &sel   = pdata->parallelSelection(sel_);

    // Saves all the coordinates (only done when making statistic, stats needs to be true!)
    if (statistic)
    {
        int   stop;
        frame frTemp;
        frTemp.vec = new rvec[iAtoms->nr];
        for (unsigned int r = 0; r < nucleotideInfo_.size(); r++)
        {
            stop = nucleotideInfo_[r]->atomStart() + nucleotideInfo_[r]->nAtoms();
            for (int i = nucleotideInfo_[r]->atomStart(); i < stop; i++)
            {
                copy_rvec(sel.coordinates()[i], frTemp.vec[i]);
            }
        }
        coords.push_back(frTemp);
    }

    GMX_RELEASE_ASSERT(NULL != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Set the default coordinates to search from
    rvec * repVec = new rvec[nucleotideInfo_.size()];
    int    i      = 0;
    for (std::vector<resInfo *>::iterator ri = nucleotideInfo_.begin();
         (ri < nucleotideInfo_.end()); ++ri, ++i)
    {
        copy_rvec(sel.coordinates()[(*ri)->atomStart()+offsetAtom], repVec[i]);
    }

    AnalysisNeighborhoodPositions pos(repVec, nucleotideInfo_.size());
    // Use neighborsearching tools
    AnalysisNeighborhoodSearch    nbsearch = nb_.initSearch(pbc, pos);
    // Find the first reference position within the cutoff.
    AnalysisNeighborhoodPair      pair;

    for (unsigned int base1 = 0; base1 < nucleotideInfo_.size(); base1++)
    {
        // Set the mode to be grid mode
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(repVec[base1]);

        // Start the search
        while (pairSearch.findNextPair(&pair))
        {

            base2 = pair.refIndex();
            if (base2 <= base1)
            {
                continue;
            }

            real        x, y, z;
            int         atomId1, atomId2;
            real        max        = 1.0;
            int         minIndex   = -1;
            adjustment = 0;

            unsigned int dimerSize = nucleotideInfo_[base1]->nAtoms() + nucleotideInfo_[base2]->nAtoms();
            // Concatenate the two arrays and set values
            for (unsigned int i = 0; i < dimerSize; i++)
            {
                // Exclude atoms with atom types which should not be used during analysis (e.g. hydrogen atoms)
                if (i < nucleotideInfo_[base1]->nAtoms())
                {
                    if (nucleotideInfo_[base1]->invalid(i))
                    {
                        adjustment++;
                        continue;
                    }

                    // Sets the first part of the arrays
                    copy_rvec(sel.coordinates()[i + nucleotideInfo_[base1]->atomStart()],
                              vec[(i - adjustment)]);
                }
                else
                {
                    if (nucleotideInfo_[base2]->invalid(i - nucleotideInfo_[base1]->nAtoms()))
                    {
                        adjustment++;
                        continue;
                    }
                    // Sets the second part of the arrays
                    copy_rvec(sel.coordinates()[i + nucleotideInfo_[base2]->atomStart() - nucleotideInfo_[base1]->nAtoms()],
                              vec[(i - adjustment)]);

                } // end else
            }     // end for

            adjustment = 0;
            for (unsigned int i = 0; i < dimerSize; i++)
            {
                if (nucleotideInfo_[base2]->nAtoms() > i)
                {
                    if (nucleotideInfo_[base2]->invalid(i))
                    {
                        adjustment++;
                        continue;
                    }
                    copy_rvec(sel.coordinates()[i +nucleotideInfo_[base2]->atomStart()], vec_swap[(i - adjustment)]);

                }
                else
                {
                    if (nucleotideInfo_[base1]->invalid(i - nucleotideInfo_[base2]->nAtoms()))
                    {
                        adjustment++;
                        continue;
                    }
                    copy_rvec(sel.coordinates()[i +nucleotideInfo_[base1]->atomStart() - nucleotideInfo_[base2]->nAtoms()],
                              vec_swap[(i - adjustment)]);
                }
            }
            dimerSize = dimerSize - adjustment;

            // Loops through the base pair template database
            dimReset  = false;
            bool bond                = false;
            bool needSwapForMoreInfo = false;

            int  base1_swap = 0;
            int  base2_swap = 0;

            for (std::vector<BasePair *>::iterator bpi = dataBase_.begin();
                 (bpi < dataBase_.end()); ++bpi)
            {

                // Check if there are missing atoms in nucleobase - if there are missing atoms, base won't be used
                if ( ((*bpi)->checkTemplate(nucleotideInfo_[base1]->residueType(),
                                            nucleotideInfo_[base2]->residueType()))/* &&
                                                                                      validateResInfo(nucleotideInfo_[base1]) &&
                                                                                      validateResInfo(nucleotideInfo_[base2]) */)
                {

                    base1_swap = base2;
                    base2_swap = base1;

                    // Get the current base size
                    if (dimerSize > (*bpi)->nrAtoms())
                    {
                        baseSize = (*bpi)->nrAtoms();
                    }
                    else
                    {
                        baseSize = dimerSize;
                    }

                    // reset the dimensions if needed
                    if (!dimReset && bRMSD)
                    {
                        reset_x_ndim(3, baseSize, NULL,
                                     baseSize, NULL, vec, (*bpi)->atomMass());
                        reset_x_ndim(3, baseSize, NULL,
                                     baseSize, NULL, vec_swap, (*bpi)->atomMass());
                        dimReset = true;
                    }

                    // check whether bonds exist for both cases (swap and non swap)
                    bool bondNonSwap = (*bpi)->isHBond(nucleotideInfo_[base1],
                                                       nucleotideInfo_[base2],
                                                       pbc, bondDist, sel);
                    bool bondSwap    = (*bpi)->isHBond(nucleotideInfo_[base1_swap],
                                                       nucleotideInfo_[base2_swap],
                                                       pbc, bondDist, sel);

                    //std::cout << "bases: " << (*bpi)->nucleotideType(0) << " " << (*bpi)->nucleotideType(1) << "; bondNonSwap = " << bondNonSwap << ", bondSwap = " << bondSwap << std::endl;

                    real RMSDNonSwap = 0;
                    real RMSDSwap    = 0;

                    // calculate the RMSD for both cases (swap and non swap)
                    if (bRMSD && (bondNonSwap || bondSwap))
                    {
                        RMSDNonSwap = (*bpi)->computeRootMeanSquareDeviation(vec, baseSize);
                        RMSDSwap    = (*bpi)->computeRootMeanSquareDeviation(vec_swap, baseSize);

                        if (RMSDNonSwap < RMSDSwap && bondNonSwap)
                        {
                            RMSD = RMSDNonSwap;
                            bond = bondNonSwap;
                        }
                        else if (bondSwap)
                        {
                            RMSD                = RMSDSwap;
                            bond                = bondSwap;
                            needSwapForMoreInfo = true;
                        }
                        else if (bondNonSwap)
                        {
                            RMSD = RMSDNonSwap;
                            bond = bondNonSwap;
                        }

                        //std::cout << "template: " << (*bpi)->bondName((*bpi)->bondType()) << "; bases: " << (*bpi)->nucleotideType(0) << nucleotideInfo_[base1]->residueNumber() << " " << (*bpi)->nucleotideType(1) << nucleotideInfo_[base2]->residueNumber() << "; bondNonSwap = " << bondNonSwap << ", bondSwap = " << bondSwap << "; RMSDNonSwap = " << RMSDNonSwap << "; RMSDSwap: " << RMSDSwap << std::endl;
                    }

                    if ((!bRMSD || (maxRMSD+(*bpi)->getTemplateRMSD() > RMSD)) && bond)
                    {
                        if (max > RMSD)
                        {
                            max            = RMSD;
                            minIndex       = bpi - dataBase_.begin();
                            if (needSwapForMoreInfo)
                            {
                                tempPair.base1 = base1_swap;
                                tempPair.base2 = base2_swap;
                            }
                            else
                            {
                                tempPair.base1 = base1;
                                tempPair.base2 = base2;
                            }

                            // Prints out list of base pairs every time the criteria is fullfilled
                            // Could be more than once for one single base pair
                            if (everyBPList)
                            {
                                if (nucleotideInfo_[base1_swap]->residueNumber() <= nucleotideInfo_[base2_swap]->residueNumber())
                                {
                                    fprintf(rnaLog, "base pair %c%i - %c%i",
                                            nucleotideInfo_[base1_swap]->residueType(),
                                            nucleotideInfo_[base1_swap]->residueNumber(),
                                            nucleotideInfo_[base2_swap]->residueType(),
                                            nucleotideInfo_[base2_swap]->residueNumber());

                                    if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Cis))
                                    {
                                        fprintf(rnaLog, "\t\t%s cis",
                                                dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Trans))
                                    {
                                        fprintf(rnaLog, "\t\t%s trans",
                                                dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                                    }
                                    if (Info)
                                    {
                                        fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                                    }
                                }
                                else
                                {
                                    fprintf(rnaLog, "base pair %c%i - %c%i", nucleotideInfo_[base2_swap]->residueType(), nucleotideInfo_[base2_swap]->residueNumber(), nucleotideInfo_[base1_swap]->residueType(), nucleotideInfo_[base1_swap]->residueNumber());

                                    if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Cis))
                                    {
                                        fprintf(rnaLog, "\t\t%s cis", dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Trans))
                                    {
                                        fprintf(rnaLog, "\t\t%s trans", dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                                    }
                                    if (Info)
                                    {
                                        fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                                    }
                                }

                                if (Info)
                                {
                                    for  (int con = 0; con < dataBase_[minIndex]->getNrBondAtoms(); con++)
                                    {
                                        atomId1 = findSameAtom(tempPair.base1, dataBase_[minIndex]->bondIndex(con, 0));
                                        atomId2 = findSameAtom(tempPair.base2, dataBase_[minIndex]->bondIndex(con, 1));

                                        atomId1 = nucleotideInfo_[tempPair.base1]->atomStart() + atomId1;
                                        atomId2 = nucleotideInfo_[tempPair.base2]->atomStart() + atomId2;

                                        // Get atom coordinates when Hbond exists and the calculate the distance

                                        x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                                        y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                                        z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                                        real currentDist = sqrt(x*x+y*y+z*z);

                                        if (currentDist <= bondDist)
                                        {
                                            fprintf(rnaLog, "bond distance %s-%s: %f\n", dataBase_[minIndex]->bondIndex(con, 0), dataBase_[minIndex]->bondIndex(con, 1), currentDist);
                                        }
                                    }
                                }
                                fprintf(rnaLog, "\n");

                            } // if (everyRMSD)
                        }     // if (max > RMSD)
                    }         // if ((!bRMSD || (maxRMSD+(*bpi)->getTemplateRMSD() > RMSD)) && bond)
                }             // if (check if there are missing atoms in base)
            }                 // for each template base pair


            // Only true if a templated base pair is found that
            // matches RMSD and hydrogen bond distance criterion
            if ((minIndex >= 0) && (max <= RMSD))
            {
                // Store the info for a new base pair hit
                tempPair.pairIndex = minIndex;

                // Prints out base pair list once (with the template that matches the base pair the best)
                if (oneBPList)
                {
                    if (nucleotideInfo_[base1_swap]->residueNumber() <= nucleotideInfo_[base2_swap]->residueNumber())
                    {
                        fprintf(rnaLog, "base pair %c%i - %c%i",
                                nucleotideInfo_[base1_swap]->residueType(),
                                nucleotideInfo_[base1_swap]->residueNumber(),
                                nucleotideInfo_[base2_swap]->residueType(),
                                nucleotideInfo_[base2_swap]->residueNumber());

                        if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Cis))
                        {
                            fprintf(rnaLog, "\t\t%s cis", dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Trans))
                        {
                            fprintf(rnaLog, "\t\t%s trans",
                                    dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                        }
                        if (Info)
                        {
                            fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                        }
                    }
                    else
                    {
                        fprintf(rnaLog, "base pair %c%i - %c%i",
                                nucleotideInfo_[base2_swap]->residueType(),
                                nucleotideInfo_[base2_swap]->residueNumber(),
                                nucleotideInfo_[base1_swap]->residueType(),
                                nucleotideInfo_[base1_swap]->residueNumber());

                        if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Cis))
                        {
                            fprintf(rnaLog, "\t\t%s cis",
                                    dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase_[minIndex]->isomerismType() == Trans))
                        {
                            fprintf(rnaLog, "\t\t%s trans",
                                    dataBase_[minIndex]->bondName(dataBase_[minIndex]->bondType()).c_str());
                        }
                        if (Info)
                        {
                            fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                        }
                    }


                    if (Info)
                    {
                        for  (int con = 0; con < dataBase_[minIndex]->getNrBondAtoms(); con++)
                        {
                            atomId1 = findSameAtom(tempPair.base1, dataBase_[minIndex]->bondIndex(con, 0));
                            atomId2 = findSameAtom(tempPair.base2, dataBase_[minIndex]->bondIndex(con, 1));

                            atomId1 = nucleotideInfo_[tempPair.base1]->atomStart() + atomId1;
                            atomId2 = nucleotideInfo_[tempPair.base2]->atomStart() + atomId2;

                            // Get atom coordinates when Hbond exists and calculate the distance
                            x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                            y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                            z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                            real currentDist = sqrt(x*x+y*y+z*z);
                            if ((currentDist <= bondDist))
                            {
                                fprintf(rnaLog, "bond distance %s-%s: %f\n", dataBase_[minIndex]->bondIndex(con, 0), dataBase_[minIndex]->bondIndex(con, 1), currentDist);
                            }
                        } // end for
                    }     // end if
                    fprintf(rnaLog, "\n");

                }           // end if
                tempPair.base1 = base1_swap;
                tempPair.base2 = base2_swap;
                // Saves the current pair
                pairs.push_back(tempPair);

            } // end if

        }     // while search for pairs

    }         // end for each residue in PDB structure

    fprintf(rnaLog, "\nFound %" GMX_PRId64 " RNA base pairs that match the given criterion.\n"
            "Double, triple, etc. count might be possible for one single base pair.\n",
            (gmx_int64_t) pairs.size());
    found.push_back(pairs);
    nbsearch.reset();
    // Magic
    dh.finishFrame();

    delete [] repVec;
}


void
RnaAnalysis::finishAnalysis(int /* nframes*/)
{
    if (NULL != rnaLog)
    {
        gmx_fio_fclose(rnaLog);
    }
}

void
RnaAnalysis::writeOutput()
{
    if ((Xpm && (DB == "bps")) || statistic)
    {
        // Set up the default colour
        int   adjustment;
        int   atomId[200];
        t_rgb rlo;
        t_rgb rhi;
        rhi.r = 1.0, rhi.g = 0.0, rhi.b = 1.0;
        rlo.r = 0.0, rlo.g = 0.0, rlo.b = 0.0;
        real * axis = new real[nucleotideInfo_.size()];

        // Set up axes for colour map, both axis are the same length
        for (unsigned int i = 0; i < nucleotideInfo_.size(); i++)
        {
            axis[i] = nucleotideInfo_[i]->residueNumber();
        }

        // Strings
        char    titel[]  = "RNA Analysis";
        char    legend[] = "legend";
        char    label[]  = "residues";
        int     nlevels  = 15;
        real  **matrix;
        matrix = new real*[nucleotideInfo_.size()];

        for (unsigned int x = 0; x < nucleotideInfo_.size(); x++)
        {
            matrix[x] = new real[nucleotideInfo_.size()];
            for (unsigned int y = 0; y < nucleotideInfo_.size(); y++)
            {
                matrix[x][y] = 0;
            }
        }

        /* Open files for storing all detected RNA base pairs */
        std::vector<FILE *>        pdbStatistic;
        std::vector<unsigned int>  model_nr;

        if (statistic)
        {
            for (std::vector<BasePair *>::iterator bpi = dataBase_.begin();
                 (bpi < dataBase_.end()); ++bpi)
            {
                char buf[256];

                // Final outpath depends on template
                if (outPath.size() > 0)
                {
                    char temp[256];
                    strcpy (temp, outPath.c_str());
                    strcat (temp, "/");

                    if (DB == "bps")
                    {
                        snprintf(buf, sizeof(buf), "%s%c%c_%s_%s.pdb",
                                 temp,
                                 (*bpi)->nucleotideType(0),
                                 (*bpi)->nucleotideType(1),
                                 (*bpi)->bondName((*bpi)->bondType()).c_str(),
                                 ((*bpi)->isomerismType() == Cis) ? "cis" : "trans");
                    }
                    else
                    {
                        snprintf(buf, sizeof(buf), "%s%c%c_temp.pdb",
                                 temp,
                                 (*bpi)->nucleotideType(0),
                                 (*bpi)->nucleotideType(1));
                    }
                }
                else
                {
                    if (DB == "bps")
                    {
                        snprintf(buf, sizeof(buf), "%c%c_%s_%s.pdb",
                                 (*bpi)->nucleotideType(0),
                                 (*bpi)->nucleotideType(1),
                                 (*bpi)->bondName((*bpi)->bondType()).c_str(),
                                 ((*bpi)->isomerismType() == Cis) ? "cis" : "trans");
                    }
                    else
                    {
                        snprintf(buf, sizeof(buf), "%c%c_temp.pdb",
                                 (*bpi)->nucleotideType(0),
                                 (*bpi)->nucleotideType(1));
                    }
                }
                pdbStatistic.push_back(gmx_fio_fopen(buf, "a+"));
                model_nr.push_back(1);
            }
        }


        // Statistic list
        std::list<frame>::const_iterator frVec = coords.begin();

        // Set up list pointers
        std::list<pairInfo>::const_iterator               posI;
        std::list<pairInfo>::const_iterator               posEnd;
        std::list <std::list <pairInfo> >::const_iterator frameI;
        std::list <std::list <pairInfo> >::const_iterator frameEnd = found.end();

        int ii = 0;
        // Start running through the whole trajectory
        for (frameI = found.begin(); frameI != frameEnd; ++frameI)
        {
            // Calculate the current marix
            posEnd = frameI->end();
            for (posI = frameI->begin(); posI != posEnd; ++posI)
            {
                // Saves the positions
                if (statistic)
                {
                    adjustment = 0;
                    unsigned int dimerSize = (nucleotideInfo_[posI->base1]->nAtoms() +
                                              nucleotideInfo_[posI->base2]->nAtoms());
                    // Concatenate the two arrays and set values to them
                    for (unsigned int i = 0; i < dimerSize; i++)
                    {
                        // Exclude invalid atom types ( hydrogen etc)
                        if (nucleotideInfo_[posI->base1]->nAtoms() > i)
                        {
                            if (nucleotideInfo_[posI->base1]->invalid(i))
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the first part of the array
                            atomId[(i - adjustment)] = i + nucleotideInfo_[posI->base1]->atomStart();
                        }
                        else
                        {
                            if (nucleotideInfo_[posI->base2]->invalid(i - nucleotideInfo_[posI->base1]->nAtoms()))
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the second part of the array
                            atomId[(i - adjustment)] = i +nucleotideInfo_[posI->base2]->atomStart() - nucleotideInfo_[posI->base1]->nAtoms();
                        }
                    }
                    dimerSize = dimerSize - adjustment;

                    // Write a pdb file
                    write_pdbfile_indexed(pdbStatistic[posI->pairIndex], " ", iAtoms,
                                          frVec->vec, dataBase_[posI->pairIndex]->getEpbc(),
                                          *dataBase_[posI->pairIndex]->box(),
                                          'A',
                                          model_nr[posI->pairIndex],
                                          dimerSize, atomId, NULL, false);
                    model_nr[posI->pairIndex]++;
                }
                if (Xpm && (DB == "bps"))
                {
                    // Get the index of the bond
                    matrix[posI->base2][posI->base1] = dataBase_[posI->pairIndex]->bondType() +1;
                    // Only got 16 colours
                    if (matrix[posI->base2][posI->base1] > 16)
                    {
                        matrix[posI->base2][posI->base1] = 17;
                    }
                    // Isomeri colour either 0 or 1
                    matrix[posI->base1][posI->base2] = dataBase_[posI->pairIndex]->isomerismType() * 10;
                }

            }
            // Print the current matrix
            if (Xpm && (DB == "bps"))
            {

                FILE             * file;
                char               temp[256] = "";

                std::ostringstream ss;
                ss << ii;
                std::string        s = ss.str();
                char const       * c = s.c_str();
                // Open a file to write to
                if (outPath != "")
                {
                    strcpy (temp, outPath.c_str());
                    strcat (temp, "/");
                    strcat (temp, c);
                    strcat (temp, "_Matrix.xpm");
                }
                else
                {
                    strcat (temp, c);
                    strcat (temp, "_Matrix.xpm");
                }
                file = gmx_fio_fopen(temp, "w");
                GMX_ASSERT((file != NULL), "Could not create an .xpm file");

                write_xpm_split(file, 0, titel, legend,
                                label, label, nucleotideInfo_.size(), nucleotideInfo_.size(),
                                axis, axis, matrix, 0.0, 16.0, &nlevels, rlo, rhi, 0.0, 16.0, &nlevels, true, rlo, rhi);

                // Reset the matix
                for (posI = frameI->begin(); posI != posEnd; ++posI)

                {
                    matrix[posI->base2][posI->base1] = 0.0;
                    matrix[posI->base1][posI->base2] = 0.0;
                }

                ii++;
                gmx_fio_fclose(file);

            }
            if (statistic)
            {
                // We dont need the current frame anymore
                delete [] frVec->vec;
                ++frVec;
            }


        }

        // Close files and release memory
        if (statistic)
        {
            for (unsigned int i = 0; i < dataBase_.size(); i++)
            {
                gmx_fio_fclose(pdbStatistic[i]);
            }
        }

        delete [] axis;
        for (unsigned int i = 0; i < nucleotideInfo_.size(); i++)
        {
            delete [] matrix[i];
        }
        delete [] matrix;

    }
    if (Xpm && (DB != "bps"))
    {
        printf("\nCalculation of the matrix is only available if the templateDB flag is set to \"bps\"");
    }
}

// find atom with the same name
int
RnaAnalysis::findSameAtom(int i, const char* name)
{
    for (unsigned int x = 0; x < nucleotideInfo_[i]->nAtoms(); x++)
    {
        if (nucleotideInfo_[i]->atomName(x).compare(name) == 0)
        {
            return x;
        }
    }
    return -1;
}

// Return true if its a valid atom, needs name and type of atom as input
bool RnaAnalysis::getAtomType(const char * name, char * type)
{
    bool set   = false;
    int  size  = strlen(name);
    for (int i = 0; i < size; i++)
    {
        // If there is an invalid atom type
        if ((!hydrogenRmsd_ && name[i] == 'H') ||
            (!phosphateRmsd_ && name[i] == 'P') ||
            (!sugarRmsd_ && name[i] == '\''))
        {
            return false;
        }
        // If its a char
        else if (isalpha(name[i]) && !set)
        {
            *type = name[i];
            set   = true;
        }
    }
    if (set)
    {
        return true;
    }
    fprintf(stderr, "no atom name found\n");
    return false;
}

const char RnaAnalysisInfo::name[]             = "rnaAnalysis";
const char RnaAnalysisInfo::shortDescription[] =
    "Analyze RNA PDB structures and trajectories";

TrajectoryAnalysisModulePointer RnaAnalysisInfo::create()
{
    return TrajectoryAnalysisModulePointer(new RnaAnalysis);
}

} // namespace analysismodules

} // namespace gmx
