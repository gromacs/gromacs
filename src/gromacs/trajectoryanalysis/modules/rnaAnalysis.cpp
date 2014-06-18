/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * \author Nina Fischer <nina.fischer@icm.uu.se>, Anders G√§rden√§s <anders.gardenas@gmail.co, Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"
#include <iostream>
#include "rnaAnalysis.h"

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <list>
#include <assert.h>

#include "config.h"

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#ifdef GMX_NATIVE_WINDOWS
#include <direct.h>
#include <io.h>
#endif

#include "gromacs/legacyheaders/copyrite.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/matio.h"

#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/uniqueptr.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/topology/atomprop.h"

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/random/random.h"

#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

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

// Declaration of class BasePair before class definition (class definition see below class RnaAnalysis)
class BasePair;

// Definition of class RnaAnalysis
class RnaAnalysis : public TrajectoryAnalysisModule
{
    public:

        //! Constructor
        RnaAnalysis();

        //! Destructor
        virtual ~RnaAnalysis();

        //! Set the options and setting
        virtual void initOptions(Options                    *options,
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

        int                                 findSameAtom(int i, const char * name);
        Selection                           sel_;
        AnalysisData                        data_;
        AnalysisDataAverageModulePointer    adata_;
        AnalysisNeighborhood                nb_;

        char                              * libfn;

        gmx_atomprop_t                      atomprop;
        // Copy and assign disallowed by base.
        // A list list of a found pairs
        std::list < std::list < pairInfo> > found;
        std::list <frame>                   coords;

        BasePair                          * dataBase;
        char			          * inPath;
        std::string                         DB;
        std::string                         outPath;
        std::string                         rnaAnalysis;
        int                                 nrTemplates;
        int                                 nrResidues;

        // Boolean options
        bool                               hydrogen;
        bool                               sugar;
        bool                               phosphate;
        bool                               everyBPList;
        bool                               statistic;
        bool                               oneBPList;
        bool                               Xpm;

        double                             extraOffset;
        double                             maxRMSD;
        double                             bondDist;
        bool   				   Info;
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

        std::string getAi() { return ai_; }

        std::string getAj() { return aj_; }
};

//  RNA base pairs
class BasePair
{
    public:

        // Constructor
        BasePair();

        // Destructor
        virtual ~BasePair();

        // Set the type of nucleotide
        void setNucleotides(const char *  nuc);

        // Get the number of atoms
        int getNrAtoms()
        {
            return atomName.size();
        }

        // Add an atom
        void addAtom(rvec x, int atomnum, std::string atomname, real m);

        // Return the atom positions
        rvec *getPos()
        {
            return atomPos;
        }

	// Add information about a pair of bonded atoms
        void addBondAtom(int ai, int aj)
        {
	    std::string ai_name, aj_name;

	    // Get the atom names corresponding to the indices
	    for (unsigned int i = 0; i < atomNumber.size(); i++)
	    {
	        if (ai == atomNumber[i])
	        {
		    ai_name = atomName[i];
	        }
	        if (aj == atomNumber[i])
	        {
		    aj_name = atomName[i];
	        }
	    }
	  
	    bondAtoms.push_back(bondAtom(ai_name, aj_name));
        }

	// Return true if it is the right template
	bool checkTemplate(resInfo * mol0, resInfo * mol1)
	{
	    return (((getNucleotides(0) == mol0->type) && (getNucleotides(1) == mol1->type)) || 
		    ((getNucleotides(0) == mol1->type) && (getNucleotides(1) == mol0->type)));
	}

        // Calculate the RMSD
        float getRMSD(rvec * vec, int baseSize);

        // Get the maximun distance in nm
        float getMaxDist()
        {
            return dist;
        }

        // Set the maximum distance to another atom
        void setAtomDist(int index);

        // Set if RMSD calculation should also hydrogen atoms
        static void setHydrogen(bool hydro)
        {
            hydrogen = hydro;
        }

        // Set if RMSD calculation should also consider sugar atoms
        static void setSugar(bool su)
        {
           sugar  = su;
        }

        // Set if RMSD calculation should also consider phoshpate atoms
        static void setPhosphate(bool phos)
        {
            phosphate = phos;
        }

        // Return phosphate
        static bool getPhosphate()
        {
            return phosphate;
        }

	// Check if a atom is valid
	bool checkAtom (const std::string & name);

        // Check if it is a valid atom to save and get the atom type
        static bool getAtomType(const char * name, char * type);

        // Return the colour with the chosen bond type representative
        t_rgb getColurCode();

        // Return nucleotide type
        char getNucleotides (int pair);

        // Set bond type
        void setBondtype (const char * name, int size);

        // Get bond type
        int getBondType();

	// Get bond name
	std::string getBondName(int type);

	// Check whether a bond exist by calculating the H bond distance
        bool getHBondDistance(int tempID1, int tempID2, float bondDist, Selection sel);
        
	// Check whether the C1' distance exceed the template threshold
	bool getC1distance(int tempID1, int tempID2, Selection sel);

        // Set the isomerism type
        void setIsomerism (std::string name);

        // Get the isomerism type
        int getIso();

        // Get the mass
        real * getMass();

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
	float getTemplateRMSD() { return templateRMSD; }

	// set the RMSD of the template
	void setTemplateRMSD(float rmsd)
	{
	    templateRMSD = rmsd;
	    return;
	}

    private:
        float                     dist;
        Nucleotide                type[2];
        std::vector <std::string> atomName;
        std::vector <int>         atomNumber;  
        rvec                    * atomPos;
        real                    * atomMass;
        Isomerism                 iso;
        BondType                  bondtype[2];
        static bool               hydrogen;
        static bool               sugar;
        static bool               phosphate;
        std::vector <bondAtom>    bondAtoms;
        matrix                    box;
        int                       ePBC;
	int                       resId[2];
	float                     templateRMSD;

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



// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.

bool BasePair::hydrogen;
bool BasePair::sugar;
bool BasePair::phosphate;

// Constructor
RnaAnalysis::RnaAnalysis()
    : TrajectoryAnalysisModule(RnaAnalysisInfo::name, RnaAnalysisInfo::shortDescription),
      adata_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(0, 1);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "rnaAnalysis");

    offsetAtom    = 0;
    libfn         = 0;
    nrTemplates   = 0;
    nrResidues    = 0;
    hydrogen      = false;
    sugar	  = false;
    phosphate     = false;
    everyBPList   = false;
    statistic     = false;
    oneBPList     = true;
    Xpm           = false;
    extraOffset   = 0.5;
    maxRMSD       = 0.25;
    Info          = false;
    iAtoms        = NULL;
    dataBase	  = NULL;
    inPath	  = NULL;
    bondDist      = 0.3;

    atomprop      = gmx_atomprop_init();
}

// Destructor
RnaAnalysis::~RnaAnalysis()
{
    // Destroy C structures where there is no automatic memory release
    // C++ takes care of memory in classes (hopefully)
    delete [] dataBase;
    delete [] nucleotideInfo;
    if (NULL != libfn)
    {
        sfree(libfn);
    }
    gmx_atomprop_destroy(atomprop);
}

void
RnaAnalysis::initOptions(Options                    *options,
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
    options->setDescription(concatenateStrings(desc));

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
                           .description("define the template DB you want to use (currentlz just bps available)"));

    // Print all results about RNA base pairs during analysis
    options->addOption(BooleanOption("bpLists").store(&everyBPList).defaultValue(false)
                           .description("print all possible base pairs during analysis"));

    // Print the result about RNA base pairs in the end
    options->addOption(BooleanOption("bpListFinal").store(&oneBPList).defaultValue(true)
                           .description("print all identified base pairs after analysis"));

    // Include hydrogen atoms in RMSD calculation
    options->addOption(BooleanOption("hydro").store(&hydrogen).defaultValue(true)
                           .description("include the hydrogen atoms of the nucleotide in RMSD calculations"));

    // Include  phosphate group atoms in RMSD calculation
    options->addOption(BooleanOption("phos").store(&phosphate).defaultValue(true)
                           .description("include the phosphate atoms of the nucleotide in RMSD calculation"));

    // Include ONLY base atoms of nucleotides in RMSD calculations
    options->addOption(BooleanOption("sugar").store(&sugar).defaultValue(true)
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
    // Store all file names of RNA base pair templates within input directory in vector pdbTemplates
    nrTemplates = 0;

    std::string str1   = "amber99sb-ildn.ff/RNA-" + DB + ".pdb";
    inPath = gmxlibfn(str1.c_str());
    
    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);
    // Set options which atoms to use for RMSD analysis
    BasePair::setHydrogen(hydrogen);
    BasePair::setSugar(sugar);
    BasePair::setPhosphate(phosphate);

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
		    nucleotideInfo[i].size		= cAtom - nucleotideInfo[i].start;
		    currentNrAtoms = 0;
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
                    taken                = true;
                    nucleotideInfo[i].start = cAtom;
                    nucleotideInfo[i].type  = *iAtoms->resinfo[rnaRes[i]].name[0];
                    nucleotideInfo[i].resNr = iAtoms->resinfo[rnaRes[i]].nr;
            	}
            	// If it is an invalid atom
            	nucleotideInfo[i].invalid[currentNrAtoms] = !BasePair::getAtomType(*iAtoms->atomname[cAtom], &type);
            	// Set atom name
            	nucleotideInfo[i].name[currentNrAtoms] = *iAtoms->atomname[cAtom];
            	currentNrAtoms++;
	    }
        }
    }

    printf("Number of RNA bases in input file: %i \n", nrResidues);

    gmx_residuetype_destroy(rt);

    // Set path and load file pointer
    float     maxDist = 0;
    char      line[100];
    FILE    * file;
    dataBase = new BasePair[100];

    // Initialize RNA template data base to compare given RNA structure to
    file = gmx_fio_fopen(inPath, "r");
    while (NULL != fgets2(line, 100, file))
    {
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
            dataBase[nrTemplates].setNucleotides(aminoType.c_str());
            if (DB == "bps")
            {
                std::string x, y, z;
                ISS >> x; ISS >> y; ISS >> z;
                dataBase[nrTemplates].setBondtype(x.c_str(), x.size());
                dataBase[nrTemplates].setIsomerism(y.c_str());
		dataBase[nrTemplates].setTemplateRMSD(std::atof(z.c_str()));
            }
        }
        else if (pdbType.compare("CRYST1") == 0)
        {
            double x, y, z;
            ISS >> x; ISS >> y; ISS >> z;
            dataBase[nrTemplates].setMatrix(x, y, z);
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
		dataBase[nrTemplates].addAtom(x, anum, aname, m);
	    }

	}
        else if (pdbType.compare("CONECT") == 0)
        {
            int ai, aj;
            ISS >> ai; ISS >> aj;
            dataBase[nrTemplates].addBondAtom(ai, aj);
        }
        else if (pdbType.compare("ENDMDL") == 0)
        {
            dataBase[nrTemplates].setEpbc(epbcNONE);
                
	    // Assign the center of mass to origin of all atoms in three dimensions
            reset_x_ndim(3, dataBase[nrTemplates].getNrAtoms(),  NULL,
                         dataBase[nrTemplates].getNrAtoms(), NULL,
                         dataBase[nrTemplates].getPos(), dataBase[nrTemplates].getMass());
		
	    // Set the distance between atoms
	    dataBase[nrTemplates].setAtomDist(offsetAtom);
	    if (dataBase[nrTemplates].getMaxDist() > maxDist)
	    {
		maxDist = dataBase[nrTemplates].getMaxDist();
	    }

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
    float                RMSD           = 0;
    int                  minIndex       = -1;
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
        printf("Using an RNA structure file or just one frame of an MD trajectory.\n");
	printf("Please, check your input trajectory if there are more than one time steps. \n\n");
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

            float       x, y, z;
            int         atomId1, atomId2;
            float       max       = 1.0;
            minIndex   = -1;
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
            totalSize = totalSize - adjustment;

            adjustment = 0;
            totalSize  = nucleotideInfo[base1].size + nucleotideInfo[base2].size;
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
            bool bond = false;
	    bool needSwapForMoreInfo = false;

            int  base1_new = 0;
            int  base2_new = 0;

            for (int i = 0; i < nrTemplates; i++)
            {

                // Check if there are missing atoms in nucleobase - if there are missing atoms, base won't be used
                if ( (dataBase[i].checkTemplate(&nucleotideInfo[base1], &nucleotideInfo[base2])) &&
                    (((nucleotideInfo[base1].type == 'A') && (nucleotideInfo[base1].size >= 30)) ||
                     ((nucleotideInfo[base1].type == 'U') && (nucleotideInfo[base1].size >= 28)) ||
                     ((nucleotideInfo[base1].type == 'G') && (nucleotideInfo[base1].size >= 31)) ||
                     ((nucleotideInfo[base1].type == 'C') && (nucleotideInfo[base1].size >= 29))) &&
                    (((nucleotideInfo[base2].type == 'A') && (nucleotideInfo[base2].size >= 30)) ||
                     ((nucleotideInfo[base2].type == 'U') && (nucleotideInfo[base2].size >= 28)) ||
                     ((nucleotideInfo[base2].type == 'G') && (nucleotideInfo[base2].size >= 32)) ||
                     ((nucleotideInfo[base2].type == 'C') && (nucleotideInfo[base2].size >= 29))))
                {

		    base1_new = base2;
		    base2_new = base1;

		    // Get the current base size
		    if (totalSize > dataBase[i].getNrAtoms())
		    {
			baseSize = dataBase[i].getNrAtoms();
		    }
		    else
		    {
			baseSize = totalSize;
		    }

		    // reset the dimensions if needed
		    if (!dimReset && bRMSD)
		    {
			reset_x_ndim(3, baseSize, NULL,
				     baseSize, NULL, vec, dataBase[i].getMass());
			reset_x_ndim(3, baseSize, NULL,
				     baseSize, NULL, vec_swap, dataBase[i].getMass());
			dimReset = true;
		    }

		    // check whether bonds exist for both cases (swap and non swap)
		    bool bondNonSwap = dataBase[i].getHBondDistance(base1, base2, bondDist, sel);
		    bool bondSwap = dataBase[i].getHBondDistance(base1_new, base2_new, bondDist, sel);

		    //std::cout << "bases: " << dataBase[i].getNucleotides(0) << " " << dataBase[i].getNucleotides(1) << "; bondNonSwap = " << bondNonSwap << ", bondSwap = " << bondSwap << std::endl;

		    float RMSDNonSwap = 0;
		    float RMSDSwap = 0;

		    // calculate the RMSD for both cases (swap and non swap)
		    if (bRMSD && (bondNonSwap || bondSwap))
		    {
			RMSDNonSwap = dataBase[i].getRMSD(vec, baseSize);
			RMSDSwap = dataBase[i].getRMSD(vec_swap, baseSize);

			if (RMSDNonSwap < RMSDSwap)
			{
			    RMSD = RMSDNonSwap;
			    bond = bondNonSwap;
			}
			else
			{
			    RMSD = RMSDSwap;
			    bond = bondSwap;
			    needSwapForMoreInfo = true;
			}

			//std::cout << "bases: " << dataBase[i].getNucleotides(0) << " " << dataBase[i].getNucleotides(1) << "; bondNonSwap = " << bondNonSwap << ", bondSwap = " << bondSwap << "; RMSD = " << RMSD << std::endl;
		    }

                    if ((!bRMSD || (maxRMSD+dataBase[i].getTemplateRMSD() > RMSD)) && bond)
                    {
			if (max > RMSD)
			{
                            max            = RMSD;
                            minIndex       = i;
			    if (needSwapForMoreInfo)
			    {
                                tempPair.base1 = base1_new;
                                tempPair.base2 = base2_new;
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
                                if (nucleotideInfo[base1_new].resNr <= nucleotideInfo[base2_new].resNr)
                                {
                                    printf("base pair %c%i - %c%i", nucleotideInfo[base1_new].type, nucleotideInfo[base1_new].resNr, nucleotideInfo[base2_new].type, nucleotideInfo[base2_new].resNr);

                                    if ((DB == "bps") && (dataBase[minIndex].getIso() == Cis))
                                    {
                                        printf("\t\t%s cis", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase[minIndex].getIso() == Trans))
                                    {
                                        printf("\t\t%s trans", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                                    }
                                    if (Info)
                                    {
                                        printf("\nRMSD to template: %f\n", RMSD);
                                    }
                                }
                                else
                                {
                                    printf("base pair %c%i - %c%i", nucleotideInfo[base2_new].type, nucleotideInfo[base2_new].resNr, nucleotideInfo[base1_new].type, nucleotideInfo[base1_new].resNr);

                                    if ((DB == "bps") && (dataBase[minIndex].getIso() == Cis))
                                    {
                                        printf("\t\t%s cis", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                                    }
                                    if ((DB == "bps") && (dataBase[minIndex].getIso() == Trans))
                                    {
                                        printf("\t\t%s trans", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                                    }
                                    if (Info)
                                    {
                                        printf("\nRMSD to template: %f\n", RMSD);
                                    }
                                }
			    
			        if (Info)
			        {
                                    for  (int con = 0; con < dataBase[minIndex].getNrBondAtoms(); con++)
                                    {
                                        atomId1 = findSameAtom(tempPair.base1, dataBase[minIndex].getBondIndex(con, 0));
                                        atomId2 = findSameAtom(tempPair.base2, dataBase[minIndex].getBondIndex(con, 1));

                                        atomId1 = nucleotideInfo[tempPair.base1].start + atomId1;
                                        atomId2 = nucleotideInfo[tempPair.base2].start + atomId2;

                                        // Get atom coordinates when Hbond exists and the calculate the distance

                                        x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                                        y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                                        z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                                        float currentDist = sqrt(x*x+y*y+z*z);

                                        if (currentDist <= bondDist)
                                        {
                                    	    printf("bond distance %s-%s: %f\n", dataBase[minIndex].getBondIndex(con, 0), dataBase[minIndex].getBondIndex(con, 1), currentDist);
                                        }
                            	    }
			        }
                                printf("\n");

                            } // if (everyRMSD)
			} // if (max > RMSD)
                    }     // if ((!bRMSD || (maxRMSD+dataBase[i].getTemplateRMSD() > RMSD)) && bond)
                }         // if (check if there are missing atoms in base)
            }             // for each template base pair


            // Only true if a templated base pair is found that
            // matches RMSD and hydrogen bond distance criterion
            if ((minIndex >= 0) && (max <= RMSD))
            {
                // Store the info for a new base pair hit
                tempPair.pairIndex = minIndex;

		// Prints out base pair list once (with the template that matches the base pair the best)
                if (oneBPList)
                {
                    if (nucleotideInfo[base1_new].resNr <= nucleotideInfo[base2_new].resNr)
                    {
                        printf("base pair %c%i - %c%i", nucleotideInfo[base1_new].type, nucleotideInfo[base1_new].resNr, nucleotideInfo[base2_new].type, nucleotideInfo[base2_new].resNr);

                        if ((DB == "bps") && (dataBase[minIndex].getIso() == Cis))
                        {
                            printf("\t\t%s cis", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase[minIndex].getIso() == Trans))
                        {
                            printf("\t\t%s trans", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                        }
                        if (Info)
                        {
                            printf("\nRMSD to template: %f\n", RMSD);
                        }
                    }
                    else
                    {
                        printf("base pair %c%i - %c%i", nucleotideInfo[base2_new].type, nucleotideInfo[base2_new].resNr, nucleotideInfo[base1_new].type, nucleotideInfo[base1_new].resNr);

                        if ((DB == "bps") && (dataBase[minIndex].getIso() == Cis))
                        {
                            printf("\t\t%s cis", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                        }
                        if ((DB == "bps") && (dataBase[minIndex].getIso() == Trans))
                        {
                            printf("\t\t%s trans", dataBase[minIndex].getBondName(dataBase[minIndex].getBondType()).c_str());
                        }
                        if (Info)
                        {
                            printf("\nRMSD to template: %f\n", RMSD);
                        }
                    }


		    if (Info)
		    {
                   	for  (int con = 0; con < dataBase[minIndex].getNrBondAtoms(); con++)
                    	{
                            atomId1 = findSameAtom(tempPair.base1, dataBase[minIndex].getBondIndex(con, 0));
                            atomId2 = findSameAtom(tempPair.base2, dataBase[minIndex].getBondIndex(con, 1));

                            atomId1 = nucleotideInfo[tempPair.base1].start + atomId1;
                            atomId2 = nucleotideInfo[tempPair.base2].start + atomId2;

                            // Get atom coordinates when Hbond exists and calculate the distance
                            x =  sel.coordinates()[atomId2][0] - sel.coordinates()[atomId1][0];
                            y =  sel.coordinates()[atomId2][1] - sel.coordinates()[atomId1][1];
                            z =  sel.coordinates()[atomId2][2] - sel.coordinates()[atomId1][2];
                            float currentDist = sqrt(x*x+y*y+z*z);
                            if ((currentDist <= bondDist))
                            {
                                printf("bond distance %s-%s: %f\n", dataBase[minIndex].getBondIndex(con, 0), dataBase[minIndex].getBondIndex(con, 1), currentDist);
                            }
                        }           // end for
		    }		    // end if
                    printf("\n");

                }           // end if
                tempPair.base1 = base1_new;
                tempPair.base2 = base2_new;
                // Saves the current pair
                pairs.push_back(tempPair);

            } // end if

        }     // while search for pairs

    }         // end for each residue in PDB structure

    printf("\nFound %lu RNA base pairs that match the given criterion.\nDouble, triple, etc. count might be possible for one single base pair.\n", pairs.size());
    found.push_back(pairs);
    nbsearch.reset();
    // Magic
    dh.finishFrame();

    delete [] repVec;
}


void
RnaAnalysis::finishAnalysis(int /* nframes*/)
{
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
        FILE ** pdbStatistic = NULL;
        int   * model_nr     = NULL;

        if (statistic)
        {
            pdbStatistic = new FILE * [nrTemplates];
            model_nr     = new int[nrTemplates];

            for (int i = 0; i < nrTemplates; i++)
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
		    	std::string s;
		    	ss << dataBase[i].getNucleotides(0) << dataBase[i].getNucleotides(1);
		    	ss >> s;
                    	strcat (temp, s.c_str());
		    	strcat (temp, "_");
	
		    	strcat (temp, dataBase[i].getBondName(dataBase[i].getBondType()).c_str());
		    	strcat (temp, "_");

		    	if (dataBase[i].getIso() == Cis)
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
		    	std::string s;
		    	ss << dataBase[i].getNucleotides(0) << dataBase[i].getNucleotides(1);
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
		    	std::string s;
		    	ss << dataBase[i].getNucleotides(0) << dataBase[i].getNucleotides(1);
		    	ss >> s;
                    	strcat (temp, s.c_str());
		    	strcat (temp, "_");

		    	strcat (temp, dataBase[i].getBondName(dataBase[i].getBondType()).c_str()),
		    	strcat (temp, "_");

		    	if (dataBase[i].getIso() == Cis)
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
		    	std::string s;
		    	ss << dataBase[i].getNucleotides(0) << dataBase[i].getNucleotides(1);
		    	ss >> s;
                    	strcat (temp, s.c_str());
		    	strcat (temp, "_");

			strcat(temp, "temp.pdb");
		    }
                }

                pdbStatistic[i] = gmx_fio_fopen(temp, "a+");
                model_nr[i]     = 1;
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
                                          frVec->vec, dataBase[posI->pairIndex].getEpbc(),
                                          *dataBase[posI->pairIndex].getMatrix(), 'A',
                                          model_nr[posI->pairIndex], totalSize, atomId, NULL, false);
                    model_nr[posI->pairIndex]++;
                }
                if (Xpm && (DB == "bps"))
                {
                    // Get the index of the bond
                    matrix[posI->base2][posI->base1] = dataBase[posI->pairIndex].getBondType() +1;
                    // Only got 16 colours
                    if (matrix[posI->base2][posI->base1] > 16)
                    {
                        matrix[posI->base2][posI->base1] = 17;
                    }
                    // Isomeri colour either 0 or 1
                    matrix[posI->base1][posI->base2] = dataBase[posI->pairIndex].getIso() * 10;
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
BasePair::BasePair()
{
    dist              = 0;
    atomPos           = NULL;
    atomMass          = NULL;
    type[0]           = Nucnull;
    type[1]           = Nucnull;
    bondtype[0]       = Bondnull;
    bondtype[1]       = Bondnull;
    iso               = Isonull;
    clear_mat(box);
    ePBC              = 0;
    resId[0]	      = -1;
    resId[1]          = -1;
}

// Destructor
BasePair::~BasePair()
{
    sfree(atomPos);
    sfree(atomMass);
}

// Calculate the RMSD
float
BasePair::getRMSD(rvec * vec, int baseSize)
{
    // TODO: add check for matching sizes of arrays
    // Set the matrix
    do_fit(baseSize, atomMass, atomPos, vec);
    // Get the rmsd
    return rmsdev(baseSize, atomMass, atomPos, vec);
}


// Add one atom.
void
BasePair::addAtom(rvec x, int atomnum, std::string atomname, real m)
{
    if (checkAtom(atomname))
    {
    	atomNumber.push_back(atomnum);
    	atomName.push_back(atomname);

	// convert coordinates from √ngstr√m to nanometer
    	x[XX] = x[XX] / 10;
    	x[YY] = x[YY] / 10;
    	x[ZZ] = x[ZZ] / 10;

	// allocate new space for the new atom position and mass
    	srenew(atomPos, atomName.size());
    	srenew(atomMass, atomName.size());

	// set the new atom position an mass
    	copy_rvec(x, atomPos[atomName.size()-1]);
    	atomMass[atomName.size()-1] = m;
    }
}


// Check whether the current atom is valid or not
bool
BasePair::checkAtom(const std::string & name)
{
    bool set = false;
    int size = name.size();
    for (int i = 0; i < size; ++i)
    {
	// If there is an invalid atom type
	if ((!BasePair::hydrogen && name[i] == 'H') || (!BasePair::phosphate && name[i] == 'P') || (!BasePair::sugar && name[i] == '\''))
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
    printf("%s is not a valid atom type!\n", name.c_str());
    return false;
}

// Return true if its a valid atom, needs name and type of atom as input
bool
BasePair::getAtomType(const char * name, char * type)
{
    bool set   = false;
    int  size  = strlen(name);
    for (int i = 0; i < size; i++)
    {
        // If there is an invalid atom type
        if ((!BasePair::hydrogen && name[i] == 'H') || (!BasePair::phosphate && name[i] == 'P') || (!BasePair::sugar && name[i] == '\''))
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
    printf("no atom name found\n");
    return false;
}

// Return the mass array of an RNA base pair
real*
BasePair::getMass()
{
    return atomMass;
}

// Set one letter code to the full name of the base
void
BasePair::setNucleotides(const char * base)
{
    for (int i = 0; i < 2; i++)
    {
        if (base[i] == 'A')
        {
            type[i] =  Adenine;
        }
        else if (base[i] == 'U')
        {
            type[i] =  Uracil;
        }
        else if (base[i] == 'T')
        {
            type[i] =  Thymine;
        }
        else if (base[i] == 'G')
        {
            type[i] =  Guanine;
        }
        else if (base[i] == 'C')
        {
            type[i] =  Cytosine;
        }
        else
        {
            printf("invalid base type %c\n", base[i]);
        }
    }

}

// Return single letter code for each base type (either the 5' or 3')
char
BasePair::getNucleotides (int i)
{
    if (i > 1)
    {
        std::cerr << "getNucleotides out of range " << std::endl;
    }
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
        return bondAtoms[bond].getAi().c_str();
    }
    else
    {
        return bondAtoms[bond].getAj().c_str();
    }
}

// Sets the maximum distance to another atom
void
BasePair::setAtomDist(int offsetAtom)
{
    float temp;
    float x, y, z;
    dist = 0;
    for (unsigned int i = 0; i < atomName.size(); i++)
    {
        x    = (atomPos[offsetAtom][0]) - (atomPos[i][0]);
        y    = (atomPos[offsetAtom][1]) - (atomPos[i][1]);
        z    = (atomPos[offsetAtom][2]) - (atomPos[i][2]);
        temp = x*x + y*y + z*z;
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
        float x, y, z;

        atomId1 = nucleotideInfo[tempId1].start + atomId1;
        atomId2 = nucleotideInfo[tempId2].start + atomId2;
        // Get atom coordinates when Hbond exists and the calculate the distance
        x =  s.coordinates()[atomId2][0] - s.coordinates()[atomId1][0];
        y =  s.coordinates()[atomId2][1] - s.coordinates()[atomId1][1];
        z =  s.coordinates()[atomId2][2] - s.coordinates()[atomId1][2];

        // Calculated the distance and store it in vector distances
        float currentDist = sqrt(x*x+y*y+z*z);
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
BasePair::getHBondDistance(int tempId1, int tempId2, float dist, Selection s)
{
    std::vector<float> distances;
    float              currentDist = 0;
    int              maxBondNumber = 0;
    if (this->getNucleotides(0) == this->getNucleotides(1))
    {
	maxBondNumber = round(this->getNrBondAtoms()/2);
    }
    else
    {
	maxBondNumber = this->getNrBondAtoms();
    }

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
            }
        }
        for (int x = 0; x < nucleotideInfo[tempId2].size; x++)
        {
            if (strcmp(name2, nucleotideInfo[tempId2].name[x]) == 0)
            {
                atomId2 = x;
            }

        }

        if (atomId1 != -1 && atomId2 != -1)
        {

            float x, y, z;
            atomId1 = nucleotideInfo[tempId1].start + atomId1;
            atomId2 = nucleotideInfo[tempId2].start + atomId2;
            // Get atom coordinates when Hbond exists and the calculate the distance
            x =  s.coordinates()[atomId2][0] - s.coordinates()[atomId1][0];
            y =  s.coordinates()[atomId2][1] - s.coordinates()[atomId1][1];
            z =  s.coordinates()[atomId2][2] - s.coordinates()[atomId1][2];

            // Calculated the distance and store it in vector distances
            currentDist = sqrt(x*x+y*y+z*z);
            distances.push_back(currentDist);
        }
    }           // end for

    if (distances.size() > 0)
    {

        // sort distances in vector from smallest to largest distance
        float* first(&distances[0]);
        float* last(first + distances.size());
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
