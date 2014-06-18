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
 * \author Anders G√§rden√§s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "rnaAnalysis.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
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
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*! \brief
 * Class used to analyze RNA structures and trajectories.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 *
 * \ingroup module_trajectoryanalysis
 */

//! Information about residue
struct resInfo
{
    //! Residue type
    char   type;
    //! Start of something
    int    start;
    //! Whether or not this particle is invalid
    bool   invalid[50];
    //! Name of the particle
    char * name[50];
    //! Residue number
    int    resNr;
    //! Size of the residue
    int    size;
};
//! Our storage of residue information
resInfo * nucleotideInfo;

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

// Declaration of class BasePair before class definition
// (class definition see below class RnaAnalysis)
class BasePair;

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


        int                                 findSameAtom(int i, const char * name);
        Selection                           sel_;
        AnalysisData                        data_;
        AnalysisDataAverageModulePointer    adata_;
        AnalysisNeighborhood                nb_;

        gmx_atomprop_t                      atomprop;
        // Copy and assign disallowed by base.
        // A list list of a found pairs
        std::list < std::list < pairInfo> > found;
        std::list <frame>                   coords;

        std::vector<BasePair *>             dataBase_;
        std::string                         DB;
        std::string                         outPath;
        std::string                         rnaAnalysis;
        std::string                         rnaLogFile;
        FILE                               *rnaLog;
        int                                 nrTemplates;
        int                                 nrResidues;

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

//! Possible nucleotides including unknown
enum Nucleotide
{
    Cytosine, Guanine, Adenine, Uracil, Thymine, Nucnull
};

//! Designation of isomers
enum Isomerism
{
    Cis, Trans, Isonull
};

//! The different dimer types
enum BondType
{
    Ww, Ss, Hh, Bh, Bs, bondSize, Bondnull
};

class bondAtom
{
    public:
        std::string ai_, aj_;

        bondAtom(std::string ai, std::string aj) : ai_(ai), aj_(aj) {}

        ~bondAtom() {}

        std::string atomI() { return ai_; }

        std::string atomJ() { return aj_; }
};

//  RNA base pairs
class BasePair
{
    public:

        // Constructor
        BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd);

        // Destructor
        virtual ~BasePair();

        // Set the types of both nucleotide encoded in a two character string
        void setNucleotides(const char *  nucleotides);

        // Return number of atoms
        int nrAtoms()
        {
            return atomName_.size();
        }

        // Add an atom
        void addAtom(rvec x, int atomnum, std::string atomname, real m);

        // Return the atom positions
        rvec *atomPosition()
        {
            if (NULL == positionPtr_)
            {
                snew(positionPtr_, atomPos_.size());
                for (unsigned int i = 0; (i < atomPos_.size()); i++)
                {
                    for (int m = 0; (m < DIM); m++)
                    {
                        positionPtr_[i][m] = atomPos_[DIM*i+m];
                    }
                }
            }
            return positionPtr_;
        }

        // Add information about a pair of bonded atoms
        void addBondAtom(int ai, int aj)
        {
            std::string ai_name, aj_name;

            // Extract the atom names corresponding to the indices
            for (unsigned int i = 0; i < atomNumber_.size(); i++)
            {
                if (ai == atomNumber_[i])
                {
                    ai_name = atomName_[i];
                }
                if (aj == atomNumber_[i])
                {
                    aj_name = atomName_[i];
                }
            }

            bondAtoms.push_back(bondAtom(ai_name, aj_name));
        }

        // Return true if it is the right template
        bool checkTemplate(resInfo * mol0, resInfo * mol1)
        {
            return (((nucleotideType(0) == mol0->type) && (nucleotideType(1) == mol1->type)) ||
                    ((nucleotideType(0) == mol1->type) && (nucleotideType(1) == mol0->type)));
        }

        // Calculate the RMSD
        real computeRootMeanSquareDeviation(rvec * vec, int baseSize);

        // Return the maximun distance in nm TODO: between what?
        real maximumDistance()
        {
            return dist;
        }

        // Set the maximum distance to another atom
        void setAtomDist(int index);

        // Return phosphate
        bool phosphateRmsd()
        {
            return phosphateRmsd_;
        }

        // Check if a atom is valid
        bool checkAtom (const std::string &name);

        // Return the colour with the chosen bond type representative
        t_rgb getColurCode();

        // Return nucleotide type. Index should be 0 or 1 or the routine will throw.
        char nucleotideType (unsigned int index);

        // Set bond type
        void setBondtype (const char * name, int size);

        // Get bond type
        int getBondType();

        // Get bond name
        std::string getBondName(int type);

        // Check whether a bond exist by calculating the H bond distance
        bool getHBondDistance(int tempID1, int tempID2, real bondDist, Selection sel);

        // Check whether the C1' distance exceed the template threshold
        bool getC1distance(int tempID1, int tempID2, Selection sel);

        // Set the isomerism type
        void setIsomerism (std::string name);

        // Get the isomerism type
        int getIso();

        // Get the mass
        real * atomMass();

        // Returns the bond index nr bond from pair
        const char * getBondIndex(int bond, int pair);

        // Returns the number of bonded atoms
        int getNrBondAtoms() { return bondAtoms.size(); }

        // Returns a pointer to the matrix box
        matrix * getMatrix()
        {
            return &box;
        }

        //! Set the matrix diagonal
        void setMatrix(double x, double y, double z)
        {
            box[XX][XX] = x; box[YY][YY] = y; box[ZZ][ZZ] = z;
        }
        // Sets ePBC
        void setEpbc(int ep)
        {
            ePBC = ep;
        }

        // Get ePBC
        int  getEpbc()
        {
            return ePBC;
        }

        // Get the RMSD of the template
        real getTemplateRMSD() { return templateRMSD; }

        // set the RMSD of the template
        void setTemplateRMSD(real rmsd)
        {
            templateRMSD = rmsd;
            return;
        }

    private:
        real                      dist;
        Nucleotide                type[2];
        std::vector <std::string> atomName_;
        std::vector <int>         atomNumber_;
        std::vector<real>         atomPos_;
        std::vector<real>         atomMass_;
        rvec                     *positionPtr_;
        Isomerism                 iso;
        BondType                  bondtype[2];
        bool                      hydrogenRmsd_;
        bool                      sugarRmsd_;
        bool                      phosphateRmsd_;
        std::vector <bondAtom>    bondAtoms;
        matrix                    box;
        int                       ePBC;
        int                       resId[2];
        real                      templateRMSD;

        // Return type of atom
        char atomType(const char name[]);

        // Get the name of the bond
        static void getBondTypeName(BondType type, char * name)
        {
            switch (type)
            {
                case Ww: name[0] = 'W'; name[1] = 'w'; return;
                case Ss: name[0] = 'S'; name[1] = 's'; return;
                case Hh: name[0] = 'H'; name[1] = 'h'; return;
                case Bs: name[0] = 'B'; name[1] = 's'; return;
                case Bh: name[0] = 'B'; name[1] = 'h'; return;
                default: fprintf( stderr, " invalid bond type %s \n", name  );
            }
        }

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
    nrTemplates        = 0;
    nrResidues         = 0;
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
    delete [] nucleotideInfo;
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
                           .description("define the template DB you want to use (currentlzy just bps available)"));

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

void
RnaAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top)
{
    // Store all file names of RNA base pair templates within input directory
    // in vector pdbTemplates
    nrTemplates = 0;

    std::string str1   = "amber99sb-ildn.ff/RNA-" + DB + ".pdb";
    const char *inPath = gmxlibfn(str1.c_str());

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

    // Initialize variables
    iAtoms = &(top.topology()->atoms);
    int                * rnaRes = new int[iAtoms->nres];
    char                 type;
    const char          *grpnames;
    int                  currentNrAtoms = 0;

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
            rnaRes[nrResidues++] = i;
        }
    }

    // Set all the values of the nucleotide
    int       cAtom = 0;
    bool      taken = false;
    nucleotideInfo = new resInfo[nrResidues];

    for (int i = 0; i < nrResidues; i++)
    {
        for (; cAtom <= iAtoms->nr; cAtom++)
        {
            // If the last atom in the inputfile has been reached
            if (cAtom == iAtoms->nr)
            {
                // If the last atom is the end of a nucleotide
                if (taken)
                {
                    taken = false;
                    // Set the size of the nucleotide
                    nucleotideInfo[i].size      = cAtom - nucleotideInfo[i].start;
                    currentNrAtoms              = 0;
                    break;
                }

                break;
            }
            // If the last atom of the current nucleotide has been reached
            else if (taken && (rnaRes[i] != iAtoms->atom[cAtom].resind))
            {
                // New nucleotide
                taken = false;
                // Set the size of the nucleotide
                nucleotideInfo[i].size                    = cAtom - nucleotideInfo[i].start;
                currentNrAtoms = 0;
                break;
            }
            // Test if the current atom is within the right nucleotide
            else if (rnaRes[i] == iAtoms->atom[cAtom].resind)
            {
                if (!taken)
                {
                    taken                   = true;
                    nucleotideInfo[i].start = cAtom;
                    nucleotideInfo[i].type  = *iAtoms->resinfo[rnaRes[i]].name[0];
                    nucleotideInfo[i].resNr = iAtoms->resinfo[rnaRes[i]].nr;
                }
                // If it is an invalid atom
                nucleotideInfo[i].invalid[currentNrAtoms] = !getAtomType(*iAtoms->atomname[cAtom], &type);
                // Set atom name
                nucleotideInfo[i].name[currentNrAtoms] = *iAtoms->atomname[cAtom];
                currentNrAtoms++;
            }
        }
    }

    fprintf(rnaLog, "Number of RNA bases in input file: %i \n", nrResidues);

    gmx_residuetype_destroy(rt);

    // Set path and load file pointer
    real      maxDist = 0;
    char      line[100];
    FILE    * file;

    // Initialize RNA template data base to compare given RNA structure to
    file = gmx_fio_fopen(inPath, "r");
    BasePair *bp = NULL;
    while (NULL != fgets2(line, 100, file))
    {
        if (NULL == bp)
        {
            bp = new BasePair(hydrogenRmsd_, sugarRmsd_, phosphateRmsd_);
        }
        std::istringstream ISS(line);
        std::string        pdbType;
        ISS >> pdbType;

        if (pdbType.compare("MODEL") == 0)
        {
            continue;
        }
        else if (pdbType.compare("REMARK") == 0)
        {
            std::string aminoType;
            ISS >> aminoType;
            bp->setNucleotides(aminoType.c_str());
            if (DB == "bps")
            {
                std::string x, y, z;
                ISS >> x; ISS >> y; ISS >> z;
                bp->setBondtype(x.c_str(), x.size());
                bp->setIsomerism(y.c_str());
                bp->setTemplateRMSD(std::atof(z.c_str()));
            }
        }
        else if (pdbType.compare("CRYST1") == 0)
        {
            double x, y, z;
            ISS >> x; ISS >> y; ISS >> z;
            bp->setMatrix(x, y, z);
        }
        else if ( (pdbType.compare("ATOM") == 0) || (pdbType.compare("HETATM") == 0))
        {
            int         anum, rnum;
            std::string aname, rname;
            rvec        x;
            real        m;
            ISS >> anum; ISS >> aname; ISS >> rname; ISS >> rnum;
            ISS >> x[XX]; ISS >> x[YY]; ISS >> x[ZZ];
            if (gmx_atomprop_query(atomprop, epropMass, "???",
                                   aname.c_str(), &m))
            {
                bp->addAtom(x, anum, aname, m);
            }

        }
        else if (pdbType.compare("CONECT") == 0)
        {
            int ai, aj;
            ISS >> ai; ISS >> aj;
            bp->addBondAtom(ai, aj);
        }
        else if (pdbType.compare("ENDMDL") == 0)
        {
            bp->setEpbc(epbcNONE);

            // Assign the center of mass to origin of all atoms in three dimensions
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
            dataBase_.push_back(std::move(bp));
            bp           = NULL;
            nrTemplates += 1;
        }
    }
    gmx_fio_fclose(file);
    // TODO: maxDist and extraOffset have disappeared.
    // Initiate the neighborsearching code
    nb_.setCutoff(maxDist + extraOffset);
    // Did not work before, but we don't know why....
    //nb_.setMode(eSearchMode_Grid);

    delete [] rnaRes;
}

void
RnaAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata)
{
    // When a base pair is found
    real                 RMSD           = 0;
    //int                  minIndex       = -1;
    bool                 dimReset       = false;
    bool                 bRMSD          = true;
    int                  baseSize       = 0;

    std::list<pairInfo>  pairs;
    pairInfo             tempPair   = {0, 0, 0};
    int                  adjustment = 0;
    rvec                 vec[200];
    rvec                 vec_swap[200];
    int                  base2 = 0;
    AnalysisDataHandle   dh    = pdata->dataHandle(data_);
    const Selection     &sel   = pdata->parallelSelection(sel_);


    // Saves all the coordinates (only done when making statistic, stats needs to be true!)
    if (statistic)
    {
        int   stop;
        frame frTemp;
        frTemp.vec = new rvec[iAtoms->nr];
        for (int r = 0; r < nrResidues; r++)
        {
            stop = nucleotideInfo[r].start + nucleotideInfo[r].size;
            for (int i = nucleotideInfo[r].start; i < stop; i++)
            {
                frTemp.vec[i][0] = sel.coordinates()[i][0];
                frTemp.vec[i][1] = sel.coordinates()[i][1];
                frTemp.vec[i][2] = sel.coordinates()[i][2];
            }
        }
        coords.push_back(frTemp);
    }

    GMX_RELEASE_ASSERT(NULL != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);
    if (fr.time == 0)
    {
        fprintf(stderr, "Using an RNA structure file or just one frame of an MD trajectory.\n");
        fprintf(stderr, "Please, check your input trajectory if there are more than one time steps. \n\n");
    }
    else
    {
        printf("Using %f frames of input MD trajectory.\n\n", fr.time);
    }

    // Set the total number of frames
    // Set the default coordinates to search from
    rvec * repVec = new rvec[nrResidues];
    for (int i = 0; i < nrResidues; i++)
    {
        repVec[i][0] = sel.coordinates()[nucleotideInfo[i].start+offsetAtom][0];
        repVec[i][1] = sel.coordinates()[nucleotideInfo[i].start+offsetAtom][1];
        repVec[i][2] = sel.coordinates()[nucleotideInfo[i].start+offsetAtom][2];
    }

    AnalysisNeighborhoodPositions pos(repVec, nrResidues);
    // Use neighborsearching tools
    AnalysisNeighborhoodSearch    nbsearch = nb_.initSearch(pbc, pos);
    // Find the first reference position within the cutoff.
    AnalysisNeighborhoodPair      pair;

    for (int base1 = 0; base1 < nrResidues; base1++)
    {

        // Set the mode to be grind mode
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

            int totalSize = nucleotideInfo[base1].size + nucleotideInfo[base2].size;
            // Concatenate the two arrays and set values
            for (int i = 0; i < totalSize; i++)
            {
                // Exclude atoms with atom types which should not be used during analysis (e.g. hydrogen atoms)
                if (nucleotideInfo[base1].size > i)
                {
                    if (nucleotideInfo[base1].invalid[i])
                    {
                        adjustment++;
                        continue;
                    }

                    // Sets the first part of the arrays
                    vec[(i - adjustment)][0] = sel.coordinates()[i +nucleotideInfo[base1].start][0];
                    vec[(i - adjustment)][1] = sel.coordinates()[i +nucleotideInfo[base1].start][1];
                    vec[(i - adjustment)][2] = sel.coordinates()[i +nucleotideInfo[base1].start][2];

                }
                else
                {
                    if (nucleotideInfo[base2].invalid[i - nucleotideInfo[base1].size])
                    {
                        adjustment++;
                        continue;
                    }
                    // Sets the second part of the arrays
                    vec[(i - adjustment)][0] = sel.coordinates()[i +nucleotideInfo[base2].start - nucleotideInfo[base1].size][0];
                    vec[(i - adjustment)][1] = sel.coordinates()[i +nucleotideInfo[base2].start - nucleotideInfo[base1].size][1];
                    vec[(i - adjustment)][2] = sel.coordinates()[i +nucleotideInfo[base2].start - nucleotideInfo[base1].size][2];

                } // end else
            }     // end for

            adjustment = 0;
            for (int i = 0; i < totalSize; i++)
            {
                if (nucleotideInfo[base2].size > i)
                {
                    if (nucleotideInfo[base2].invalid[i])
                    {
                        adjustment++;
                        continue;
                    }
                    vec_swap[(i - adjustment)][0] = sel.coordinates()[i +nucleotideInfo[base2].start][0];
                    vec_swap[(i - adjustment)][1] = sel.coordinates()[i +nucleotideInfo[base2].start][1];
                    vec_swap[(i - adjustment)][2] = sel.coordinates()[i +nucleotideInfo[base2].start][2];

                }
                else
                {
                    if (nucleotideInfo[base1].invalid[i - nucleotideInfo[base2].size])
                    {
                        adjustment++;
                        continue;
                    }
                    vec_swap[(i - adjustment)][0] = sel.coordinates()[i +nucleotideInfo[base1].start - nucleotideInfo[base2].size][0];
                    vec_swap[(i - adjustment)][1] = sel.coordinates()[i +nucleotideInfo[base1].start - nucleotideInfo[base2].size][1];
                    vec_swap[(i - adjustment)][2] = sel.coordinates()[i +nucleotideInfo[base1].start - nucleotideInfo[base2].size][2];
                }
            }
            totalSize = totalSize - adjustment;

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
                if ( ((*bpi)->checkTemplate(&nucleotideInfo[base1], &nucleotideInfo[base2])) &&
                     (((nucleotideInfo[base1].type == 'A') && (nucleotideInfo[base1].size >= 30)) ||
                      ((nucleotideInfo[base1].type == 'U') && (nucleotideInfo[base1].size >= 28)) ||
                      ((nucleotideInfo[base1].type == 'G') && (nucleotideInfo[base1].size >= 31)) ||
                      ((nucleotideInfo[base1].type == 'C') && (nucleotideInfo[base1].size >= 29))) &&
                     (((nucleotideInfo[base2].type == 'A') && (nucleotideInfo[base2].size >= 30)) ||
                      ((nucleotideInfo[base2].type == 'U') && (nucleotideInfo[base2].size >= 28)) ||
                      ((nucleotideInfo[base2].type == 'G') && (nucleotideInfo[base2].size >= 32)) ||
                      ((nucleotideInfo[base2].type == 'C') && (nucleotideInfo[base2].size >= 29))))
                {

                    base1_swap = base2;
                    base2_swap = base1;

                    // Get the current base size
                    if (totalSize > (*bpi)->nrAtoms())
                    {
                        baseSize = (*bpi)->nrAtoms();
                    }
                    else
                    {
                        baseSize = totalSize;
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
                    bool bondNonSwap = (*bpi)->getHBondDistance(base1, base2, bondDist, sel);
                    bool bondSwap    = (*bpi)->getHBondDistance(base1_swap, base2_swap, bondDist, sel);

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

                        //std::cout << "template: " << (*bpi)->getBondName((*bpi)->getBondType()) << "; bases: " << (*bpi)->nucleotideType(0) << nucleotideInfo[base1].resNr << " " << (*bpi)->nucleotideType(1) << nucleotideInfo[base2].resNr << "; bondNonSwap = " << bondNonSwap << ", bondSwap = " << bondSwap << "; RMSDNonSwap = " << RMSDNonSwap << "; RMSDSwap: " << RMSDSwap << std::endl;
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
                                if (nucleotideInfo[base1_swap].resNr <= nucleotideInfo[base2_swap].resNr)
                                {
                                    fprintf(rnaLog, "base pair %c%i - %c%i",
                                            nucleotideInfo[base1_swap].type,
                                            nucleotideInfo[base1_swap].resNr,
                                            nucleotideInfo[base2_swap].type,
                                            nucleotideInfo[base2_swap].resNr);

                                    if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Cis))
                                    {
                                        fprintf(rnaLog, "\t\t%s cis",
                                                dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Trans))
                                    {
                                        fprintf(rnaLog, "\t\t%s trans",
                                                dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                                    }
                                    if (Info)
                                    {
                                        fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                                    }
                                }
                                else
                                {
                                    fprintf(rnaLog, "base pair %c%i - %c%i", nucleotideInfo[base2_swap].type, nucleotideInfo[base2_swap].resNr, nucleotideInfo[base1_swap].type, nucleotideInfo[base1_swap].resNr);

                                    if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Cis))
                                    {
                                        fprintf(rnaLog, "\t\t%s cis", dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Trans))
                                    {
                                        fprintf(rnaLog, "\t\t%s trans", dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
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
                                        atomId1 = findSameAtom(tempPair.base1, dataBase_[minIndex]->getBondIndex(con, 0));
                                        atomId2 = findSameAtom(tempPair.base2, dataBase_[minIndex]->getBondIndex(con, 1));

                                        atomId1 = nucleotideInfo[tempPair.base1].start + atomId1;
                                        atomId2 = nucleotideInfo[tempPair.base2].start + atomId2;

                                        // Get atom coordinates when Hbond exists and the calculate the distance

                                        x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                                        y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                                        z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                                        real currentDist = sqrt(x*x+y*y+z*z);

                                        if (currentDist <= bondDist)
                                        {
                                            fprintf(rnaLog, "bond distance %s-%s: %f\n", dataBase_[minIndex]->getBondIndex(con, 0), dataBase_[minIndex]->getBondIndex(con, 1), currentDist);
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
                    if (nucleotideInfo[base1_swap].resNr <= nucleotideInfo[base2_swap].resNr)
                    {
                        fprintf(rnaLog, "base pair %c%i - %c%i",
                                nucleotideInfo[base1_swap].type,
                                nucleotideInfo[base1_swap].resNr,
                                nucleotideInfo[base2_swap].type,
                                nucleotideInfo[base2_swap].resNr);

                        if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Cis))
                        {
                            fprintf(rnaLog, "\t\t%s cis", dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Trans))
                        {
                            fprintf(rnaLog, "\t\t%s trans",
                                    dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                        }
                        if (Info)
                        {
                            fprintf(rnaLog, "\nRMSD to template: %f\n", RMSD);
                        }
                    }
                    else
                    {
                        fprintf(rnaLog, "base pair %c%i - %c%i",
                                nucleotideInfo[base2_swap].type,
                                nucleotideInfo[base2_swap].resNr,
                                nucleotideInfo[base1_swap].type,
                                nucleotideInfo[base1_swap].resNr);

                        if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Cis))
                        {
                            fprintf(rnaLog, "\t\t%s cis",
                                    dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase_[minIndex]->getIso() == Trans))
                        {
                            fprintf(rnaLog, "\t\t%s trans",
                                    dataBase_[minIndex]->getBondName(dataBase_[minIndex]->getBondType()).c_str());
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
                            atomId1 = findSameAtom(tempPair.base1, dataBase_[minIndex]->getBondIndex(con, 0));
                            atomId2 = findSameAtom(tempPair.base2, dataBase_[minIndex]->getBondIndex(con, 1));

                            atomId1 = nucleotideInfo[tempPair.base1].start + atomId1;
                            atomId2 = nucleotideInfo[tempPair.base2].start + atomId2;

                            // Get atom coordinates when Hbond exists and calculate the distance
                            x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                            y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                            z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                            real currentDist = sqrt(x*x+y*y+z*z);
                            if ((currentDist <= bondDist))
                            {
                                fprintf(rnaLog, "bond distance %s-%s: %f\n", dataBase_[minIndex]->getBondIndex(con, 0), dataBase_[minIndex]->getBondIndex(con, 1), currentDist);
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
        real * axis = new real[nrResidues];

        // Set up axes for colour map, both axis are the same length
        for (int i = 0; i < nrResidues; i++)
        {
            axis[i] = nucleotideInfo[i].resNr;
        }

        // Strings
        char    titel[]  = "RNA Analysis";
        char    legend[] = "legend";
        char    label[]  = "residues";
        int     nlevels  = 15;
        real  **matrix;
        matrix = new real*[nrResidues];

        for (int x = 0; x < nrResidues; x++)
        {
            matrix[x] = new real[nrResidues];
            for (int y = 0; y < nrResidues; y++)
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
                char temp[256] = "";

                // Final outpath depends on template
                if (outPath != "")
                {
                    strcpy (temp, outPath.c_str());
                    strcat (temp, "/");

                    if (DB == "bps")
                    {
                        std::stringstream ss;
                        std::string       s;
                        ss << (*bpi)->nucleotideType(0) << (*bpi)->nucleotideType(1);
                        ss >> s;
                        strcat (temp, s.c_str());
                        strcat (temp, "_");

                        strcat (temp, (*bpi)->getBondName((*bpi)->getBondType()).c_str());
                        strcat (temp, "_");

                        if ((*bpi)->getIso() == Cis)
                        {
                            strcat (temp, "cis.pdb");
                        }
                        else
                        {
                            strcat (temp, "trans.pdb");
                        }
                    }
                    else
                    {
                        std::stringstream ss;
                        std::string       s;
                        ss << (*bpi)->nucleotideType(0) << (*bpi)->nucleotideType(1);
                        ss >> s;
                        strcat (temp, s.c_str());
                        strcat (temp, "_");

                        strcat(temp, "temp.pdb");
                    }
                }
                else
                {
                    if (DB == "bps")
                    {
                        std::stringstream ss;
                        std::string       s;
                        ss << (*bpi)->nucleotideType(0) << (*bpi)->nucleotideType(1);
                        ss >> s;
                        strcat (temp, s.c_str());
                        strcat (temp, "_");

                        strcat (temp, (*bpi)->getBondName((*bpi)->getBondType()).c_str()),
                        strcat (temp, "_");

                        if ((*bpi)->getIso() == Cis)
                        {
                            strcat (temp, "cis.pdb");
                        }
                        else
                        {
                            strcat (temp, "trans.pdb");
                        }
                    }
                    else
                    {
                        std::stringstream ss;
                        std::string       s;
                        ss << (*bpi)->nucleotideType(0) << (*bpi)->nucleotideType(1);
                        ss >> s;
                        strcat (temp, s.c_str());
                        strcat (temp, "_");

                        strcat(temp, "temp.pdb");
                    }
                }

                pdbStatistic.push_back(gmx_fio_fopen(temp, "a+"));
                model_nr.push_back(1);
            }
        }


        // Statistic list
        std::list<frame>::const_iterator frVec = coords.begin();

        // Set up list pointers
        std::list<pairInfo>::const_iterator                 posI;
        std::list<pairInfo>::const_iterator                 posEnd;
        std::list < std::list < pairInfo> >::const_iterator frameI;
        std::list < std::list < pairInfo> >::const_iterator frameEnd = found.end();

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
                    int totalSize = nucleotideInfo[posI->base1].size + nucleotideInfo[posI->base2].size;
                    // Concatenate the two arrays and set values to them
                    for (int i = 0; i < totalSize; i++)
                    {
                        // Exclude invalid atom types ( hydrogen etc)
                        if (nucleotideInfo[posI->base1].size > i)
                        {
                            if (nucleotideInfo[posI->base1].invalid[i])
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the first part of the array
                            atomId[(i - adjustment)] = i + nucleotideInfo[posI->base1].start;
                        }
                        else
                        {
                            if (nucleotideInfo[posI->base2].invalid[i - nucleotideInfo[posI->base1].size])
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the second part of the array
                            atomId[(i - adjustment)] = i +nucleotideInfo[posI->base2].start - nucleotideInfo[posI->base1].size;
                        }
                    }
                    totalSize = totalSize - adjustment;

                    // Write a pdb file
                    write_pdbfile_indexed(pdbStatistic[posI->pairIndex], " ", iAtoms,
                                          frVec->vec, dataBase_[posI->pairIndex]->getEpbc(),
                                          *dataBase_[posI->pairIndex]->getMatrix(), 'A',
                                          model_nr[posI->pairIndex], totalSize, atomId, NULL, false);
                    model_nr[posI->pairIndex]++;
                }
                if (Xpm && (DB == "bps"))
                {
                    // Get the index of the bond
                    matrix[posI->base2][posI->base1] = dataBase_[posI->pairIndex]->getBondType() +1;
                    // Only got 16 colours
                    if (matrix[posI->base2][posI->base1] > 16)
                    {
                        matrix[posI->base2][posI->base1] = 17;
                    }
                    // Isomeri colour either 0 or 1
                    matrix[posI->base1][posI->base2] = dataBase_[posI->pairIndex]->getIso() * 10;
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

                if (file == NULL)
                {
                    std::cerr << "Could not create an .xpm file" << std::endl;
                }

                write_xpm_split(file, 0, titel, legend,
                                label, label, nrResidues, nrResidues,
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
            for (int i = 0; i < nrTemplates; i++)
            {
                gmx_fio_fclose(pdbStatistic[i]);
            }
        }

        delete [] axis;
        for (int i = 0; i < nrResidues; i++)
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
    for (int x = 0; x < nucleotideInfo[i].size; x++)
    {
        if (strcmp(name, nucleotideInfo[i].name[x]) == 0)
        {
            return x;
        }
    }
    return -1;
}



// Constructor
BasePair::BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd)
{
    dist              = 0;
    type[0]           = Nucnull;
    type[1]           = Nucnull;
    bondtype[0]       = Bondnull;
    bondtype[1]       = Bondnull;
    iso               = Isonull;
    clear_mat(box);
    ePBC              = 0;
    resId[0]          = -1;
    resId[1]          = -1;
    templateRMSD      = 0;
    hydrogenRmsd_     = hydrogenRmsd;
    sugarRmsd_        = sugarRmsd;
    phosphateRmsd_    = phosphateRmsd;
    positionPtr_      = NULL;
}

// Destructor
BasePair::~BasePair()
{
    sfree(positionPtr_);
}

// Calculate the RMSD
real
BasePair::computeRootMeanSquareDeviation(rvec * vec, int baseSize)
{
    // TODO: add check for matching sizes of arrays
    // Set the matrix
    do_fit(baseSize, atomMass(), atomPosition(), vec);
    // Get the rmsd
    return rmsdev(baseSize, atomMass(), atomPosition(), vec);
}


// Add one atom.
void
BasePair::addAtom(rvec x, int atomnum, std::string atomname, real m)
{
    if (checkAtom(atomname))
    {
        atomNumber_.push_back(atomnum);
        atomName_.push_back(atomname);

        // convert coordinates from ≈ngstrˆm to nanometer
        // set the new atom position and mass
        for (unsigned int mm = 0; (mm < DIM); mm++)
        {
            atomPos_.push_back(x[mm]*0.1);
        }
        atomMass_.push_back(m);
    }
}


// Check whether the current atom is valid or not
bool
BasePair::checkAtom(const std::string &name)
{
    bool set  = false;
    int  size = name.size();
    for (int i = 0; i < size; ++i)
    {
        // If there is an invalid atom type
        if ((!hydrogenRmsd_ && name[i] == 'H') ||
            (!phosphateRmsd_ && name[i] == 'P') ||
            (!sugarRmsd_ && name[i] == '\''))
        {
            return false;
        }
        // If it's a char
        else if (isalpha(name[i]) && !set)
        {
            set = true;
        }
    }
    if (set)
    {
        return true;
    }
    fprintf(stderr, "%s is not a valid atom type!\n", name.c_str());
    return false;
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

// Return the mass array of an RNA base pair
real*
BasePair::atomMass()
{
    return &atomMass_[0];
}

// Set one letter code to the full name of the base
void
BasePair::setNucleotides(const char * bases)
{
    for (int i = 0; i < 2; i++)
    {
        if (bases[i] == 'A')
        {
            type[i] =  Adenine;
        }
        else if (bases[i] == 'U')
        {
            type[i] =  Uracil;
        }
        else if (bases[i] == 'T')
        {
            type[i] =  Thymine;
        }
        else if (bases[i] == 'G')
        {
            type[i] =  Guanine;
        }
        else if (bases[i] == 'C')
        {
            type[i] =  Cytosine;
        }
        else
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "invalid base type %c\n", bases[i]);
            GMX_THROW(InvalidInputError(buf));
        }
    }

}

// Return single letter code for each base type (either the 5' or 3')
char
BasePair::nucleotideType (unsigned int i)
{
    GMX_ASSERT(i < 2, "nucleotideType out of range");
    if (type[i] == Adenine)
    {
        return 'A';
    }
    else if (type[i] == Uracil)
    {
        return 'U';
    }
    else if (type[i] == Thymine)
    {
        return 'T';
    }
    else if (type[i] == Guanine)
    {
        return 'G';
    }
    else
    {
        return 'C';
    }
}

// Returns the bond index nr bond from pair
const char*
BasePair::getBondIndex(int bond, int pair)
{
    if (pair == 0)
    {
        return bondAtoms[bond].atomI().c_str();
    }
    else
    {
        return bondAtoms[bond].atomJ().c_str();
    }
}

// Sets the maximum distance to another atom
void
BasePair::setAtomDist(int offsetAtom)
{
    real  temp;
    rvec *x = atomPosition();
    dist = 0;
    for (unsigned int i = 0; i < atomName_.size(); i++)
    {
        temp = distance2(x[offsetAtom], x[i]);
        if (temp > dist)
        {
            dist = temp;
        }
    }
    dist = sqrt(dist);
}

// Get the type of bonds
int
BasePair::getBondType()
{
    int temp = 0;

    //sum formula start at 5
    temp = bondSize * bondtype[0];
    return (temp + bondtype[1]);
}

// Get the name of bonds
std::string
BasePair::getBondName(int type)
{
    switch (type)
    {
        case 0: return "WwWw";
        case 1: return "WwSs";
        case 2: return "WwHh";
        case 3: return "WwBh";
        case 4: return "WwBs";
        case 5: return "SsWw";
        case 6: return "SsSs";
        case 7: return "SsHh";
        case 8: return "SsBh";
        case 9: return "SsBs";
        case 10: return "HhWw";
        case 11: return "HhSs";
        case 12: return "HhHh";
        case 13: return "HhBh";
        case 14: return "HhBs";
        case 15: return "BhWw";
        case 16: return "BhSs";
        case 17: return "BhHh";
        case 18: return "BhBh";
        case 19: return "BhBs";
        case 20: return "BsWw";
        case 21: return "BsSs";
        case 22: return "BsHh";
        case 23: return "BsBh";
        case 24: return "BsBs";
        default: fprintf(stderr, " invalid bond type number!\n"); return "ERROR";
    }
}

// Set bond type
void
BasePair::setBondtype (const char * name, int size)
{
    int  found = 0;
    char bName[2];
    for (int i = 0; i < size - 1; i++)
    {
        for (int f = 0; f < bondSize; f++)
        {
            getBondTypeName(static_cast<BondType>(f), bName);
            if (name[i] == bName[0] && name[i+1] == bName[1])
            {
                this->bondtype[found++] = static_cast<BondType>(f);
                if (found == 2)
                {
                    return;
                }
            }
        }
    }
}


bool
BasePair::getC1distance(int tempId1, int tempId2, Selection s)
{
    int atomId1 = -1;
    int atomId2 = -1;

    for (int x = 0; x < nucleotideInfo[tempId1].size; x++)
    {
        if (strcmp(nucleotideInfo[tempId1].name[x], "C1'") == 0)
        {
            atomId1 = x;
        }
    }
    for (int x = 0; x < nucleotideInfo[tempId2].size; x++)
    {
        if (strcmp(nucleotideInfo[tempId2].name[x], "C1'") == 0)
        {
            atomId2 = x;
        }
    }

    if (atomId1 == -1 || atomId2 == -1)
    {
        std::cerr << "No C1' atom(s) found in base pair." << std::endl;
        return false;
    }
    else
    {
        real x, y, z;

        atomId1 = nucleotideInfo[tempId1].start + atomId1;
        atomId2 = nucleotideInfo[tempId2].start + atomId2;
        // Get atom coordinates when Hbond exists and the calculate the distance
        x =  s.coordinates()[atomId2][0] - s.coordinates()[atomId1][0];
        y =  s.coordinates()[atomId2][1] - s.coordinates()[atomId1][1];
        z =  s.coordinates()[atomId2][2] - s.coordinates()[atomId1][2];

        // Calculated the distance and store it in vector distances
        real currentDist = sqrt(x*x+y*y+z*z);
        if (currentDist < 1.2 /* templateC1Dist */)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool
BasePair::getHBondDistance(int tempId1, int tempId2, real dist, Selection s)
{
    std::vector<real>  distances;
    real               currentDist   = 0;
    int                maxBondNumber = this->getNrBondAtoms();

    for  (int con = 0; con < maxBondNumber; con++)
    {

        const char * name1 = this->getBondIndex(con, 0);
        const char * name2 = this->getBondIndex(con, 1);

        int          atomId1 = -1;
        int          atomId2 = -1;

        for (int x = 0; x < nucleotideInfo[tempId1].size; x++)
        {
            if (strcmp(name1, nucleotideInfo[tempId1].name[x]) == 0)
            {
                atomId1 = x;
                //std::cout << "Atom 1 " << name1 << " " << atomId1 << std::endl;
            }
        }
        for (int x = 0; x < nucleotideInfo[tempId2].size; x++)
        {
            if (strcmp(name2, nucleotideInfo[tempId2].name[x]) == 0)
            {
                atomId2 = x;
                //std::cout << "Atom 2 " << name2 << " " << atomId2 << std::endl;
            }

        }

        if (atomId1 != -1 && atomId2 != -1)
        {

            real x, y, z;
            atomId1 = nucleotideInfo[tempId1].start + atomId1;
            atomId2 = nucleotideInfo[tempId2].start + atomId2;
            // Get atom coordinates when Hbond exists and the calculate the distance
            x =  s.coordinates()[atomId2][0] - s.coordinates()[atomId1][0];
            y =  s.coordinates()[atomId2][1] - s.coordinates()[atomId1][1];
            z =  s.coordinates()[atomId2][2] - s.coordinates()[atomId1][2];

            // Calculated the distance and store it in vector distances
            currentDist = sqrt(x*x+y*y+z*z);
            distances.push_back(currentDist);

            //std::cout << "Atom 1 " << name1 << " " << atomId1 << " Atom 2 " << name2 << " " << atomId2 << std::endl;
        }
    }           // end for

    if (distances.size() > 0)
    {

        // sort distances in vector from smallest to largest distance
        real* first(&distances[0]);
        real* last(first + distances.size());
        std::sort(first, last);



        if (distances[0] <= dist)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

// Set the isomerism type
void
BasePair::setIsomerism (std::string name)
{
    // if "cis" is in name
    if (name.compare("cis") == 0)
    {
        iso = Cis;
    }
    else
    {
        iso = Trans;
    }
}

// Get the isomerism type
int
BasePair::getIso()
{
    return iso;
}




}       // end of namespace

const char RnaAnalysisInfo::name[]             = "rnaAnalysis";
const char RnaAnalysisInfo::shortDescription[] =
    "Analyze RNA PDB structures and trajectories";

TrajectoryAnalysisModulePointer RnaAnalysisInfo::create()
{
    return TrajectoryAnalysisModulePointer(new RnaAnalysis);
}

} // namespace analysismodules

} // namespace gmx
