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
 * \author Nina Fischer <nina.fischer@icm.uu.se>, Anders Gärdenäs <anders.gardenas@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
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
#include "gromacs/utility/gmxassert.h"

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

// Information about residue
struct resInfo
{
    char   type;
    int    start;
    bool   invalid[50];
    char * name[50];
    int    resNr;
    int    size;
};
resInfo * residueInfo;

// Info about a base pair, base1 and base2 is the index of resInfo
struct pairInfo
{
    int base1;
    int base2;
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

        int                                findSameAtom(int i, char * name);
        Selection                          sel_;
        AnalysisData                       data_;
        AnalysisDataAverageModulePointer   adata_;
        AnalysisNeighborhood               nb_;

        char                             * libfn;

        // Copy and assign disallowed by base.
        // A list list of a found pairs
        std::list < std::list < pairInfo> > found;
        std::list <frame>                   coords;

        BasePair                          * dataBase;
        std::vector <std::string>           pdbTemplates;
        std::string                         inPath;
        std::string                         DB;
        std::string                         outPath;
        std::string                         rnaAnalysis;
        int                                 nrTemplates;
        int                                 nrResidues;

        // Boolean options
        bool                               hydrogen;
        bool                               onlyAmino;
        bool                               phosphate;
        bool                               printResult;
        bool                               statistic;
        bool                               everyRMSD;
        bool                               Xpm;
        bool                               bRMSD;

        double                             extraOffset;
        double                             maxRMSD;
        double                             bondDist;
        t_atoms                          * iAtoms;
        int                                offsetAtom;
};


enum Nucleotide
{
    Cytosine, Guanine, Adenine, Uracil, Thymine, Nucnull
};

enum Isomerism
{
    Cis, Trans, Isonull
};

enum BondType
{
    Ww, Ss, Hh, Bh, Bs, bondSize, Bondnull
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

        // Set the number of atoms and initialize their mass
        void setNrAtoms(int i)
        {
            nrAtoms     = i;
            this->atoms = new rvec[i];
            this->mass  = new real[i];
        }

        // Get the number of atoms
        int getNrAtoms()
        {
            return nrAtoms;
        }

        // Set the type and position of the atom also updates the current atom
        void setAtomPos(rvec vec);

        // Return the atom positions
        rvec *getPos()
        {
            return atoms;
        }

        // Return true if it is the same base pair
        bool sameNucleotides(resInfo * mol0, resInfo * mol1)
        {
            if (getNucleotides(0) ==  mol0->type && getNucleotides(1) ==  mol1->type)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        // If the order of the bases in a base pair should be flipped (it they are the same type, it should be true)
        bool flip(resInfo * mol0, resInfo * mol1)
        {
            return (getNucleotides(0) !=  mol0->type && getNucleotides(1) !=  mol1->type );
        }

        // Calculate the RMSD
        float getRMSD(rvec * vec, int baseSize);

        // Set the vector for the atoms
        void setMolecule(rvec * vec, t_atoms * atoms);

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

        // Set if RMSD calculation should only consider base atoms
        static void setOnlyAmino(bool amino)
        {
            onlyAmino = amino;
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

        // Return conect rows
        gmx_conect getConect()
        {
            return conections;
        }

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

        bool getHBondDistance(int tempID1, int tempID2, float bondDist, Selection sel);
        bool getC1distance(int tempID1, int tempID2, Selection sel);
        //float calculateRMSD(int tempId1, int tempId2, Selection s);

        // Set the isomerism type
        void setIsomerism (const char * name, int size);

        // Get the isomerism type
        int getIso();

        // Get the mass
        real * getMass();

        // Set the name of the atom bond (can be done multiple times)
        void setBondAtom(char * name1, char * name2);

        // Returns the bond index nr bond from pair
        char * getBondIndex(int bond, int pair);

        // Returns the number of bonded atoms
        int getNrBondAtoms();

        // Returns a pointer to the matrix box
        matrix * getMatrix()
        {
            return &box;
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

    private:
        int         nrAtoms;
        int         currentNrAtoms;
        float       dist;
        Nucleotide  type[2];
        rvec      * atoms;
        real      * mass;
        Isomerism   iso;
        BondType    bondtype[2];
        static bool hydrogen;
        static bool onlyAmino;
        static bool phosphate;
        char      * bondAtoms[500][2];
        int         nrBondAtoms;
        gmx_conect  conections;
        matrix      box;
        int         ePBC;

        // Return type of atom
        char atomType(const char name[]);

        // Set the mass at index i
        void setMass(char type);

        // Updates the nr of attoms (nrAtoms = currentNrAtoms)
        void UpdateNrAtoms();

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
bool BasePair::onlyAmino;
bool BasePair::phosphate;
//resInfo * residueInfo;

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
    dataBase      = NULL;
    nrTemplates   = 0;
    nrResidues    = 0;
    hydrogen      = false;
    onlyAmino     = false;
    phosphate     = false;
    printResult   = false;
    statistic     = false;
    everyRMSD     = false;
    Xpm           = false;
    bRMSD         = true;
    extraOffset   = 0.5;
    maxRMSD       = 0.25;
    iAtoms        = NULL;
    bondDist      = 0.3;
    //matrix        = NULL;

}

// Destructor
RnaAnalysis::~RnaAnalysis()
{
    // Destroy C structures where there is no automatic memory release
    // C++ takes care of memory in classes (hopefully)
    delete [] dataBase;
    delete [] residueInfo;
    if (NULL != libfn)
    {
        sfree(libfn);
    }

}

void
RnaAnalysis::initOptions(Options                    *options,
                         TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Reads and analyzes an RNA structure to detect all bases that are in contact (e.g., form base pairs),",
        "this is can be done in two ways:",
        "1. check if one hydrogen bond between two bases matches a certain distance criteria (use <hbonds> as templateDB).",
        "2. compare RNA structure bases to a certain base pair template database (use <bps> as templateDB).",
        "The data is displayed on the terminal as raw data and also as xtc file",
        "The xtc file can be coverted to an eps file displaying a matrix describing all base pairs."

    };

    // Add the descriptive text (program help text) to the options
    options->setDescription(concatenateStrings(desc));

    // Add option for optional output file
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&rnaAnalysis).defaultBasename("result")
                           .description("output file" ));

    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select").valueCount(1)
                           .store(&sel_).required().onlyAtoms()
                           .description("Use default (all); select which atoms should be used for analysis"));

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);

    options->addOption(StringOption("templateDB").store(&DB).defaultValue("")
                           .description("define the template DB you want to use (bps or hbonds)"));

    // Output path
    options->addOption(StringOption("outPath").store(&outPath).defaultValue("")
                           .description("output path (default current directory)"));

    // Print all results about RNA base pairs
    options->addOption(BooleanOption("allResult").store(&everyRMSD).defaultValue(false)
                           .description("print all possible base pairs and corresponding RMSD values after analysis"));

    // Include hydrogen atoms in RMSD calculation
    options->addOption(BooleanOption("hydro").store(&hydrogen).defaultValue(false)
                           .description("include hydrogen atoms in RMSD calculations"));

    // Include  phosphate group atoms in RMSD calculation
    options->addOption(BooleanOption("phos").store(&phosphate).defaultValue(false)
                           .description("include phosphate atoms in RMSD calculation"));

    // Include ONLY base atoms of nucleotides in RMSD calculations
    options->addOption(BooleanOption("base").store(&onlyAmino).defaultValue(false)
                           .description("include only the base atoms of the nucleotide in RMSD calculation"));

    // Print the result while running
    options->addOption(BooleanOption("pr").store(&printResult).defaultValue(false)
                           .description("print all identified base pairs and their characteristics during analysis"));

    // No output file
    options->addOption(BooleanOption("xpmFile").store(&Xpm).defaultValue(false)
                           .description("matrix xpm output file is generated by default"));

    // RMSD
    options->addOption(BooleanOption("RMSD").store(&bRMSD).defaultValue(true)
                           .description("calculate RMSD between template base pair and detected base pair"));

    // Saves coordinates of base pairs
    options->addOption(BooleanOption("outPDB").store(&statistic).defaultValue(false)
                           .description("save coordinates of all identified RNA base pairs in PDB files"));

    // Minium RMSD for identified base pars
    options->addOption(DoubleOption("maxRMSD").store(&maxRMSD).defaultValue(0.25)
                           .description("maximum RMSD (cut-off value) for identifying base pairs"));

    // Extra offset used when looking for base pairs
    options->addOption(DoubleOption("addOffset").store(&extraOffset).defaultValue(0.5)
                           .description("additional RMSD offset used when looking for base pairs"));

    // Distance of hydrogen bonds
    options->addOption(DoubleOption("bondDist").store(&bondDist).defaultValue(0.3)
                           .description("distance between hydrogen bond donor and acceptor atom"));
}

void
RnaAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top)
{
    // Store all file names of RNA base pair templates within input directory in vector pdbTemplates
    nrTemplates = 0;

#ifdef HAVE_DIRENT_H

    DIR           *dir = NULL;
    struct dirent *ent = NULL;

    // Input path option
    std::string str1 = "amber99sb-ildn.ff/RNAtemplates/" + DB + "/";
    inPath = gmxlibfn(str1.c_str());
    std::cout << inPath << std::endl;

    if ((dir = opendir (inPath.c_str())) != NULL)
    {
        while ((ent = readdir (dir)) != NULL)
        {
            if ((strcmp(ent->d_name, ".") != 0) && (strcmp(ent->d_name, "..") != 0))
            {
                nrTemplates++;
                pdbTemplates.push_back(ent->d_name);
            }
        }
        closedir (dir);
    }

#elif (defined GMX_NATIVE_WINDOWS)

    //    struct gmx_directory *  gmxdir;
    intptr_t             windows_handle = 0;
    struct _finddata_t   finddata       = 0;
    int                  first          = 0;
    long                 hFile          = 0;

    // Input path option
    char     * dirname  = "amber99sb-ildn.ff/RNAtemplates/" + DB + "/";
    int        len;
    char     * tmpname = 0;

    len = strlen(dirname);
    snew(tmpname, len+3);

    strncpy(tmpname, dirname, len+1);

    /* Remove possible trailing directory separator */
    if (tmpname[len] == '/' || tmpname[len] == '\\')
    {
        tmpname[len] = '\0';
    }

    /* Add wildcard */
    strcat(tmpname, "/*");

    hFile = _findfirst("*.*", &finddata);
    /*
    //windows_handle->first = 1;
    if ( (windows_handle = _findfirst(tmpname, &windows_handel->finddata)) != NULL)
    {
        nrTemplates++;
        pdbTemplates.push_back(finddata.name());
    }
    while (_findnext(hFile, &c_file) == 0)
    {
        nrTemplates++;
        pdbTemplates.push_back(finddata.name());
    }
    _findclose(hFile);
    */

#else
    {
        std::cerr << "error reading database" << std::endl;
    }
#endif

    // Allocate array with length of all possible template base pairs
    dataBase = new BasePair[nrTemplates];

    std::cout << nrTemplates << std::endl;

    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);
    // Set options which atoms to use for RMSD analysis
    BasePair::setHydrogen(hydrogen);
    BasePair::setOnlyAmino(onlyAmino);
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
    int                 rnaRes[iAtoms->nres];
    char                type;
    const char         *grpnames;
    int                 currentNrAtoms = 0;

    // Set up the residue type
    //gmx_residuetype_t rt = NULL;
    gmx_residuetype_t *rt;
    gmx_residuetype_init(&rt);
    assert(rt);

    // Get all the RNA residues
    nrResidues = 0;
    for (int i = 0; i < iAtoms->nres; i++)
    {
        gmx_residuetype_get_type(rt, *iAtoms->resinfo[i].name, &grpnames);
        if (strcmp("RNA", grpnames) == 0)
        {
            rnaRes[nrResidues++] = i;
        }
    }

    // Set all the values of the residues
    int       cAtom = 0;
    bool      taken = false;
    //resInfo * residueInfo;
    residueInfo = new resInfo[nrResidues];
    for (int i = 0; i < nrResidues; i++)
    {
        for (; cAtom < iAtoms->nr; cAtom++)
        {
            // Test if the current atom is within the right residue
            if (rnaRes[i] == iAtoms->atom[cAtom].resind && !taken)
            {
                taken                = true;
                residueInfo[i].start = cAtom;
                residueInfo[i].type  = *iAtoms->resinfo[rnaRes[i]].name[0];
                residueInfo[i].resNr = iAtoms->resinfo[rnaRes[i]].nr;
            }
            // If the last atom in the resiude has been found
            else if (taken && (rnaRes[i] != iAtoms->atom[cAtom].resind || cAtom == iAtoms->nr-1))
            {
                // New residue
                taken = false;
                // Also test the last one
                residueInfo[i].invalid[currentNrAtoms] = !BasePair::getAtomType(*iAtoms->atomname[cAtom], &type);
                residueInfo[i].name[currentNrAtoms]    = *iAtoms->atomname[cAtom];
                residueInfo[i].size                    = cAtom - residueInfo[i].start;
                currentNrAtoms = 0;
                break;
            }
            // If it is an invalid atom
            residueInfo[i].invalid[currentNrAtoms] = !BasePair::getAtomType(*iAtoms->atomname[cAtom], &type);
            // Set atom name
            residueInfo[i].name[currentNrAtoms] = *iAtoms->atomname[cAtom];
            currentNrAtoms++;
        }
    }

    std::cout << "Number of RNA residues in input file: " << nrResidues << std::endl;


    // Set path and load file pointer
    float   maxDist = 0;
    int     nrAtoms[2];
    char    aminoType[2];
    char    tempPath[2048];
    char    tempTitel[1000];
    FILE  * file;
    int     nameSize;
    int     pbc = 0;

    // Initialize RNA template data base to compare given RNA structure to
    for (int i = 0; i < nrTemplates; i++)
    {
        // Set up the path to RNA templates
        strcpy (tempPath, inPath.c_str());
        strcat (tempPath, pdbTemplates[i].c_str());
        file = gmx_fio_fopen(tempPath, "r");

        if (file == NULL)
        {
            std::cerr << "Could not open template database files " << i << std::endl;
        }

        // Count the number of atoms
        get_pdb_coordnum(file, nrAtoms);
        gmx_fio_fclose(file);

        // Initialize atoms and vector
        t_atoms * atoms = &(top.topology()->atoms);
        snew(atoms, 1);
        init_t_atoms(atoms, nrAtoms[0], TRUE);
        rvec vec[nrAtoms[0]];

        // Read the pdb file containing the RNA template structure
        read_pdb_conf(tempPath, tempTitel, atoms, vec, &pbc, *dataBase[i].getMatrix(), false, dataBase[i].getConect());
        dataBase[i].setEpbc(pbc);
        // Set the the molecule values
        aminoType[0] = *atoms->resinfo[0].name[0];
        aminoType[1] = *atoms->resinfo[1].name[0];
        dataBase[i].setNucleotides(aminoType);

        // Set the types of bonds and isomerism, wont be needed if there is no file to save
        nameSize = strlen(pdbTemplates[i].c_str());

        if (Xpm)
        {
            dataBase[i].setIsomerism(pdbTemplates[i].c_str(), nameSize);
            dataBase[i].setBondtype(pdbTemplates[i].c_str(), nameSize);
        }
        // Removes and set all the vectors
        dataBase[i].setNrAtoms(nrAtoms[0]);
        dataBase[i].setMolecule(vec, atoms);
        // Assign the center of mass to origin of all atoms in three dimensions
        reset_x_ndim(3, dataBase[i].getNrAtoms(),  NULL,
                     dataBase[i].getNrAtoms(), NULL,
                     dataBase[i].getPos(), dataBase[i].getMass());

        // Set the distance between atoms
        dataBase[i].setAtomDist(offsetAtom);
        if (dataBase[i].getMaxDist() > maxDist)
        {
            maxDist = dataBase[i].getMaxDist();
        }

        int in1 = 0;
        int in2 = 0;

        for (int j = 0; j < nrAtoms[0]; j++)
        {
            for (int k = j; k < nrAtoms[0]; k++)
            {
                in1 = atoms->pdbinfo[j].atomnr;
                in2 = atoms->pdbinfo[k].atomnr;

                if ((in1 < in2) && (gmx_conect_exist(dataBase[i].getConect(), in1, in2)))
                {
                    dataBase[i].setBondAtom(*atoms->atomname[j+1], *atoms->atomname[k+1]);
                }

                if ((in2 < in1) && (gmx_conect_exist(dataBase[i].getConect(), in2, in1)))
                {
                    dataBase[i].setBondAtom(*atoms->atomname[k+1], *atoms->atomname[j+1]);
                }

            }
        }
    }           //end for

    // Initiate the neighborsearching code
    nb_.setCutoff(maxDist + extraOffset);
    // Did not work before, but we don't know why....
    //nb_.setMode(eSearchMode_Grid);
}

void
RnaAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata)
{
    // When a base pair is found
    float                RMSD           = 0;
    //float                max            = maxRMSD;
    int                  minIndex       = -1;
    bool                 dimReset       = false;
    int                  baseSize       = 0;

    std::list<pairInfo>  pairs;
    pairInfo             tempPair   = {0, 0, 0};
    int                  adjustment = 0;
    rvec                 vec[200];
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
            stop = residueInfo[r].start + residueInfo[r].size;
            for (int i = residueInfo[r].start; i < stop; i++)
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
        std::cout << "Using a PDB structure file or just one frame of an MD trajectory." << std::endl << "Please, check your input trajectory." << std::endl << std::endl;
    }
    else
    {
        std::cout << "Using " << fr.time << " frames of input MD trajectory." << std::endl;
    }

    // Set the total number of frames
    // Set the default coordinates to search from
    rvec repVec[nrResidues];
    for (int i = 0; i < nrResidues; i++)
    {
        repVec[i][0] = sel.coordinates()[residueInfo[i].start+offsetAtom][0];
        repVec[i][1] = sel.coordinates()[residueInfo[i].start+offsetAtom][1];
        repVec[i][2] = sel.coordinates()[residueInfo[i].start+offsetAtom][2];
    }

    AnalysisNeighborhoodPositions pos(repVec, nrResidues);
    // Use neighborsearching tools
    AnalysisNeighborhoodSearch    nbsearch = nb_.initSearch(pbc, pos);
    // Find the first reference position within the cutoff.
    AnalysisNeighborhoodPair      pair;

    int tempBase1 = 0;
    int tempBase2 = 0;
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

            /*
               std::cout << "residue info " << residueInfo[base1].type << " " << residueInfo[base2].type << std::endl;
               std::cout << "residue info temp " << residueInfo[tempBase1].type << " " << residueInfo[tempBase2].type << std::endl;


                // Order base pairs
                if (residueInfo[base1].type == 'C' || (residueInfo[base1].type > residueInfo[base2].type && residueInfo[base2].type != 'C'))
                {
                    tempBase1 = base2;
                    tempBase2 = base1;
                }
                else
                {
                    tempBase1 = base1;
                    tempBase2 = base2;
                }

               std::cout << "residue info after " << residueInfo[base1].type << " " << residueInfo[base2].type << std::endl;
               std::cout << "residue info temp after " << residueInfo[tempBase1].type << " " << residueInfo[tempBase2].type << std::endl;
             */

            tempBase1 = base1;
            tempBase2 = base2;

            int         nrRuns = 1;
            float       x, y, z;
            int         atomId1, atomId2;
            float       max       = maxRMSD;
            minIndex  = -1;
            //std::vector<std::string> basePairNames;

            // If both bases are the same (e.g. A and A) the following has to be executed twice
            //if (residueInfo[base1].type == residueInfo[base2].type)
            if (residueInfo[tempBase1].type == residueInfo[tempBase2].type)
            {
                nrRuns = 2;
            }

            for (int run = 0; run < nrRuns; run++)
            {
                /*
                      if (run == 1)
                      {
                          // Change the order if both are the same base type
                          int temp = tempBase2;
                          tempBase2 = tempBase1;
                          tempBase1 = temp;
                   //std::cout << "change order of bases " << run << std::endl;
                      }
                 */
                adjustment = 0;

                int totalSize = residueInfo[tempBase1].size + residueInfo[tempBase2].size;
                // Concatenate the two arrays and set values
                for (int i = 0; i < totalSize; i++)
                {
                    // Exclude invalid atom types ( hydrogen etc)
                    if (residueInfo[tempBase1].size > i)
                    {
                        if (residueInfo[tempBase1].invalid[i])
                        {
                            adjustment++;
                            continue;
                        }

                        // Sets the first part of the arrays
                        vec[(i - adjustment)][0] = sel.coordinates()[i +residueInfo[tempBase1].start][0];
                        vec[(i - adjustment)][1] = sel.coordinates()[i +residueInfo[tempBase1].start][1];
                        vec[(i - adjustment)][2] = sel.coordinates()[i +residueInfo[tempBase1].start][2];

                    }           // end if
                    else
                    {
                        if (residueInfo[tempBase2].invalid[i - residueInfo[tempBase1].size])
                        {
                            adjustment++;
                            continue;
                        }
                        // Sets the second part of the arrays
                        vec[(i - adjustment)][0] = sel.coordinates()[i +residueInfo[tempBase2].start - residueInfo[tempBase1].size][0];
                        vec[(i - adjustment)][1] = sel.coordinates()[i +residueInfo[tempBase2].start - residueInfo[tempBase1].size][1];
                        vec[(i - adjustment)][2] = sel.coordinates()[i +residueInfo[tempBase2].start - residueInfo[tempBase1].size][2];

                    }           // end else
                }
                totalSize = totalSize - adjustment;


                // Loops through the base pair template database
                //float RMSD      = 0;
                dimReset  = false;
                bool bond = false;
                //bool bond = true;

                //bool C1dist = false;
                //C1dist = getC1distance(tempBase1, tempBase2, sel);

                for (int i = 0; i < nrTemplates; i++)
                {
                    //std::cout << "for loop templates " << i << " " << dataBase[i].getNrAtoms() << std::endl;

                    // If it is the same base pair (.. else .. continue)
                    //if (dataBase[i].sameNucleotides(&residueInfo[base1], &residueInfo[base2]))
                    if (dataBase[i].sameNucleotides(&residueInfo[tempBase1], &residueInfo[tempBase2]))
                    {

                        //std::cout << residueInfo[base1].type << residueInfo[base1].resNr << " " << residueInfo[base2].type << residueInfo[base2].resNr << " - " << residueInfo[tempBase1].type << residueInfo[tempBase1].resNr << " " << residueInfo[tempBase2].type << residueInfo[tempBase2].resNr << std::endl;
                        //std::cout << residueInfo[tempBase1].type << residueInfo[tempBase1].resNr << " " << residueInfo[tempBase2].type << residueInfo[tempBase2].resNr << std::endl;
                        //std::cout << "template " << dataBase[i].getNucleotides(0) << " " << dataBase[i].getNucleotides(1)<< std::endl;


                        /*
                                  // Flip current template base pair (e.g. UA to AU)
                                  if (dataBase[i].flip(&residueInfo[tempBase1], &residueInfo[tempBase2]))
                                  {
                                      int temp = tempBase2;
                                      tempBase2 = tempBase1;
                                      tempBase1 = temp;
                           //std::cout << "flip " << tempBase1 << " " << tempBase2 << std::endl;
                                  }
                         */
                        // Check if there is an atom missing
                        if (totalSize > dataBase[i].getNrAtoms())
                        {
                            baseSize = dataBase[i].getNrAtoms();
                        }
                        else
                        {
                            baseSize = totalSize;
                        }
                        // Reset the position to fit orgin to be able to calculate RMSD
                        if (!dimReset && bRMSD)
                        {
                            /*
                               int totalSize = residueInfo[tempBase1].size + residueInfo[tempBase2].size;
                               // Concatenate the two arrays and set values
                               for (int i = 0; i < totalSize; i++)
                               {
                                std::cout << *vec[i] << std::endl;
                               }
                             */
                            reset_x_ndim(3, baseSize,  NULL,
                                         baseSize, NULL, vec, dataBase[i].getMass());
                            dimReset = true;
                        }

                        //C1dist = dataBase[i].getC1distance(tempBase1, tempBase2, sel);


                        // bond is set to true if it matches the bond distance criterion
                        bond = dataBase[i].getHBondDistance(tempBase1, tempBase2, bondDist, sel);

                        //std::cout << bond << std::endl;

                        //float rms = dataBase[i].calculateRMSD(tempBase1, tempBase2, sel);

                        if (bRMSD && bond)
                        {
                            // Calculates RMSD
                            RMSD = dataBase[i].getRMSD(vec, baseSize);
                            //std::cout << "RMSD " << RMSD << std::endl;
                        }

                        if ((!bRMSD || (max > RMSD)) && bond)
                        {
                            max            = RMSD;
                            minIndex       = i;
                            tempPair.base1 = tempBase1;
                            tempPair.base2 = tempBase2;

                            if (everyRMSD)
                            {
                                for  (int con = 0; con < dataBase[minIndex].getNrBondAtoms(); con++)
                                {
                                    atomId1 = findSameAtom(tempPair.base1, dataBase[minIndex].getBondIndex(con, 0));
                                    atomId2 = findSameAtom(tempPair.base2, dataBase[minIndex].getBondIndex(con, 1));

                                    atomId1 = residueInfo[tempPair.base1].start + atomId1;
                                    atomId2 = residueInfo[tempPair.base2].start + atomId2;

                                    // Get atom coordinates when Hbond exists and the calculate the distance
                                    x =  sel.coordinates()[atomId1][0] - sel.coordinates()[atomId2][0];
                                    y =  sel.coordinates()[atomId1][1] - sel.coordinates()[atomId2][1];
                                    z =  sel.coordinates()[atomId1][2] - sel.coordinates()[atomId2][2];
                                    float currentDist = sqrt(x*x+y*y+z*z);
                                    //if ((currentDist <= bondDist) && (currentDist != 0))
                                    if (currentDist <= bondDist)
                                    {
                                        std::cout << "bond distance " << dataBase[minIndex].getBondIndex(con, 0) << "-" << dataBase[minIndex].getBondIndex(con, 1) << " " << currentDist << std::endl;
                                    }
                                }
                                std::cout << "base pair " << residueInfo[tempBase1].type << residueInfo[tempBase1].resNr << " - " << residueInfo[base2].type << residueInfo[base2].resNr << "     template " << pdbTemplates[minIndex].c_str() << " rmsd " << RMSD << std::endl << std::endl;
                            }


                        }
                    }
                }
            }

            // Only true if a templated base pair is found that
            // matches RMSD and hydrogen bond distance criterion
            if ((minIndex >= 0) && (max <= RMSD))
            {
                // Store the info for a new base pair hit
                tempPair.pairIndex = minIndex;

                if (printResult)
                {
                    for  (int con = 0; con < dataBase[minIndex].getNrBondAtoms(); con++)
                    {
                        atomId1 = findSameAtom(tempPair.base1, dataBase[minIndex].getBondIndex(con, 0));
                        atomId2 = findSameAtom(tempPair.base2, dataBase[minIndex].getBondIndex(con, 1));

                        atomId1 = residueInfo[tempPair.base1].start + atomId1;
                        atomId2 = residueInfo[tempPair.base2].start + atomId2;

                        // Get atom coordinates when Hbond exists and calculate the distance
                        x =  sel.coordinates()[atomId1][0] - sel.coordinates()[atomId2][0];
                        y =  sel.coordinates()[atomId1][1] - sel.coordinates()[atomId2][1];
                        z =  sel.coordinates()[atomId1][2] - sel.coordinates()[atomId2][2];
                        float currentDist = sqrt(x*x+y*y+z*z);
                        //if ((currentDist <= bondDist) && (currentDist != 0))
                        if (currentDist <= bondDist)
                        {
                            std::cout << "bond distance " << dataBase[minIndex].getBondIndex(con, 0) << "-" << dataBase[minIndex].getBondIndex(con, 1) << " " << currentDist << std::endl;
                        }
                    }
                    std::cout << "base pair " << residueInfo[tempBase1].type << residueInfo[tempBase1].resNr << " - " << residueInfo[tempBase2].type << residueInfo[tempBase2].resNr << "     template " << pdbTemplates[minIndex].c_str() << " rmsd " << max << std::endl << std::endl;
                }
                tempPair.base1 = base1;
                tempPair.base2 = base2;
                // Saves the current pair
                pairs.push_back(tempPair);
            }
        }
    }
    std::cout << std::endl << "Found " << pairs.size() << " RNA base pairs that match the given criterion." << std::endl << "Double, triple, etc. count might be possible for one single base pair." << std::endl;
    found.push_back(pairs);
    nbsearch.reset();
    // Magic
    dh.finishFrame();
}


void
RnaAnalysis::finishAnalysis(int /* nframes*/)
{
}

void
RnaAnalysis::writeOutput()
{
    if (Xpm || statistic)
    {
        // Set up the default colour
        int   adjustment;
        int   atomId[200];
        t_rgb rlo;
        t_rgb rhi;
        rhi.r = 1.0, rhi.g = 0.0, rhi.b = 1.0;
        rlo.r = 0.0, rlo.g = 0.0, rlo.b = 0.0;
        real axis[nrResidues];

        // Set up axes for colour map, both axis are the same length
        for (int i = 0; i < nrResidues; i++)
        {
            axis[i] = residueInfo[i].resNr;
        }

        // Strings
        char    titel[]  = "RNA Analysis";
        char    legend[] = "legend";
        char    label[]  = "residues";
        int     nlevels  = 15;

        // Make a x*x matrix
        real * matrix[nrResidues];
        for (int i = 0; i < nrResidues; i++)
        {
            for (int j = 0; j < nrResidues; j++)
            {
                matrix[i][j] = 0.0;
            }
        }
        //FILE  * file;
        //char    temp[256];
        if (Xpm)
        {
            for (int x = 0; x < nrResidues; x++)
            {
                matrix[x] = new real[nrResidues];
                for (int y = 0; y < nrResidues; y++)
                {
                    matrix[x][y] = 0;
                }
            }

            /*
                // Open a file to write to write matrix xpm file
               strcpy (temp, outPath.c_str());
                strcat (temp, "outMatrix.xpm");
               file = gmx_fio_fopen(temp, "w");
                if (file == NULL)
                {
                std::cerr << "could not create an .xpm file" << std::endl;
                }
             */
        }


        /* Open files for storing all detected RNA base pairs */
        FILE ** pdbStatistic = NULL;
        int   * model_nr     = NULL;
        //char  * chainid;

        if (statistic)
        {
            char temp[256];

            pdbStatistic = new FILE * [nrTemplates];
            model_nr     = new int[nrTemplates];
            //chainid      = new char[nrTemplates];

            //int nr = 1;

            for (int i = 0; i < nrTemplates; i++)
            {
                strcpy (temp, outPath.c_str());
                strcat (temp, "/");
                strcat (temp, pdbTemplates[i].c_str());
                pdbStatistic[i] = gmx_fio_fopen(temp, "a+");
                model_nr[i]     = 1;
                //nr++;
                //std::cout << "PDB File " << dataBase[i].getNucleotides(0) << std::endl;
                //chainid[i]      = 'A';
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
                    int totalSize = residueInfo[posI->base1].size + residueInfo[posI->base2].size;
                    // Concatenate the two arrays and set values to them
                    for (int i = 0; i < totalSize; i++)
                    {
                        // Exclude invalid atom types ( hydrogen etc)
                        if (residueInfo[posI->base1].size > i)
                        {
                            if (residueInfo[posI->base1].invalid[i])
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the first part of the array
                            atomId[(i - adjustment)] = i + residueInfo[posI->base1].start;
                        }
                        else
                        {
                            if (residueInfo[posI->base2].invalid[i - residueInfo[posI->base1].size])
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the second part of the array
                            atomId[(i - adjustment)] = i +residueInfo[posI->base2].start - residueInfo[posI->base1].size;
                        }
                    }
                    totalSize = totalSize - adjustment;

                    // Write a pdb file
                    write_pdbfile_indexed(pdbStatistic[posI->pairIndex], " ", iAtoms,
                                          frVec->vec, dataBase[posI->pairIndex].getEpbc(),
                                          *dataBase[posI->pairIndex].getMatrix(), 'A',
                                          model_nr[posI->pairIndex], totalSize, atomId, NULL, false);
                    model_nr[posI->pairIndex]++;
                    //model_nr[posI->pairIndex], totalSize, atomId, NULL, false);
                    //model_nr[posI->pairIndex]++;


                }
                if (Xpm)
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
            if (Xpm)
            {

                FILE             * file;
                char               temp[256];

                std::ostringstream ss;
                ss << ii;
                std::string        s = ss.str();
                char const       * c = s.c_str();
                // Open a file to write to
                strcpy (temp, outPath.c_str());
                strcat (temp, c);
                strcat (temp, "_outMatrix.xpm");
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
                for (int i = 0; i < nrResidues; i++)
                {
                    delete [] matrix[i];
                }
            }
            if (statistic)
            {
                // We dont need the current frame anymore
                delete [] frVec->vec;
                ++frVec;
            }


        }

        // Close files and release memory
        /*
           if (Xpm)
           {
            //gmx_fio_fclose(file);
            for (int i = 0; i < nrResidues; i++)
            {
                delete [] matrix[i];
            }
           }
         */
        if (statistic)
        {
            for (int i = 0; i < nrTemplates; i++)
            {
                gmx_fio_fclose(pdbStatistic[i]);
            }
            //delete [] model_nr;
        }
    }
}


// find atom with the same name
int
RnaAnalysis::findSameAtom(int i, char* name)
{
    for (int x = 0; x < residueInfo[i].size; x++)
    {
        if (strcmp(name, residueInfo[i].name[x]) == 0)
        {
            return x;
        }
    }
    return -1;
}



// Constructor
BasePair::BasePair()
{
    conections        = gmx_conect_init();
    currentNrAtoms    = 0;
    nrAtoms           = 0;
    dist              = 0;
    nrBondAtoms       = 0;
    mass              = 0;
    type[0]           = Nucnull;
    type[1]           = Nucnull;
    bondtype[0]       = Bondnull;
    bondtype[1]       = Bondnull;
    iso               = Isonull;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            box[i][j]         = 0.0;
        }
    }
    for (int i = 0; i < 500; i++)
    {
        bondAtoms[i][0] = NULL;
        bondAtoms[i][1] = NULL;
    }
    ePBC              = 0;
}

// Destructor
BasePair::~BasePair()
{
    gmx_conect_done(conections);
}

// Calculate the RMSD
float
BasePair::getRMSD(rvec * vec, int baseSize)
{
    // Set the matrix
    do_fit(baseSize, this->mass, this->getPos(), vec);
    // Get the rmsd
    return rmsdev(baseSize, this->mass, this->getPos(), vec);
}


// Set the vector of the atoms in "atoms"
// if onlyAmino is used, only the atoms belonging to a base will be saved
void
BasePair::setMolecule(rvec * vec, t_atoms * atoms)
{
    char name = 0;
    for (int i = 0; i < this->nrAtoms; i++)

    {
        // Get structure type and the set its posisions
        if (getAtomType(*atoms->atomname[i], &name))
        {
            this->setMass(name);
            this->setAtomPos(vec[i]);
        }
    }
    this->UpdateNrAtoms();
}

// Return true if its a valid atom, needs name and type of atom as input
bool
BasePair::getAtomType(const char * name, char * type)
{
    bool set   = false;
    bool amino = false;
    int  size  = strlen(name);
    for (int i = 0; i < size; i++)
    {
        // If there is an invalid atom type
        if ((!BasePair::hydrogen && name[i] == 'H') || (!BasePair::phosphate && name[i] == 'P'))
        {
            return false;
        }
        // Only if its a amino
        if (name[i] == '\'')
        {
            amino = true;
        }
        // If its a char
        else if (isalpha(name[i]) && !set)
        {
            *type = name[i];
            set   = true;
        }
    }
    // If the atom is not part of the base
    if (onlyAmino && !amino)
    {
        return false;
    }
    else if (set)
    {
        return true;
    }
    std::cout << "no atom name found" << std::endl;
    return false;
}




// Set the mass to index i
void
BasePair::setMass(char type)
{
    if (type == 'C')
    {
        this->mass[currentNrAtoms] =  12.0;
    }
    else if (type == 'O')
    {
        this->mass[currentNrAtoms] =  16.0;
    }
    else if (type == 'P')
    {
        this->mass[currentNrAtoms] =  30.0;
    }
    else if (type == 'H')
    {
        this->mass[currentNrAtoms] =  1.0;
    }
    else if (type == 'N')
    {
        this->mass[currentNrAtoms] =  14.0;
    }
    else
    {
        this->mass[currentNrAtoms] =  1.0;
        std::cerr << "invalid atom type" << std::endl;
    }
}

// Return the mass of an RNA base pair
real*
BasePair::getMass()
{
    return mass;
}

// Set the position of the atoms and uppdates also the current atom
void
BasePair::setAtomPos(rvec vec)
{
    atoms[currentNrAtoms][0]   = vec[0];
    atoms[currentNrAtoms][1]   = vec[1];
    atoms[currentNrAtoms++][2] = vec[2];
}


// Updates the number of atoms (nrAtoms = currentNrAtoms)
void
BasePair::UpdateNrAtoms()
{
    nrAtoms = currentNrAtoms;
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
            std::cout << "invalid base type " << base[i] << std::endl;
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

// Sets the corresponding Hbond from the BasePair and residueInfo
void
BasePair::setBondAtom(char * name1, char * name2)
{
    bondAtoms[nrBondAtoms][0]   = name1;
    bondAtoms[nrBondAtoms++][1] = name2;
}

// Returns the number of bonded atoms
int
BasePair::getNrBondAtoms()
{
    return nrBondAtoms;
}

// Returns the bond index nr bond from pair
char*
BasePair::getBondIndex(int bond, int pair)
{
    return bondAtoms[bond][pair];
}

// Sets the maximum distance to another atom
void
BasePair::setAtomDist(int offsetAtom)
{
    float temp;
    float x, y, z;
    dist = 0;
    for (int i = 0; i < nrAtoms; i++)
    {
        x    = (atoms[offsetAtom][0]) - (atoms[i][0]);
        y    = (atoms[offsetAtom][1]) - (atoms[i][1]);
        z    = (atoms[offsetAtom][2]) - (atoms[i][2]);
        temp = x*x+y*y+z*z;
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
    //always get the lowest one
    if (bondtype[0] < bondtype[1])
    {
        //sum formula start at 6 then -2
        temp = ((bondSize+1) - bondtype[0])*bondtype[0];
        return (temp + bondtype[1]);
    }
    else
    {
        //sum formula start at 6 then -2
        temp = ((bondSize+1) - bondtype[1])*bondtype[1];
        return (temp + bondtype[0]);
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

    for (int x = 0; x < residueInfo[tempId1].size; x++)
    {
        //      std::cout << residueInfo[tempId1].name[x] << std::endl;
        if (strcmp(residueInfo[tempId1].name[x], "C1'") == 0)
        {
            //std::cout << residueInfo[tempId1].name[x] << std::endl;

            atomId1 = x;
        }
    }
    for (int x = 0; x < residueInfo[tempId2].size; x++)
    {
        if (strcmp(residueInfo[tempId2].name[x], "C1'") == 0)
        {
            //std::cout << residueInfo[tempId2].name[x] << std::endl;
            atomId2 = x;
        }
    }
    //std::cout << atomId1 << " " << atomId2 << std::endl;
    if (atomId1 == -1 || atomId2 == -1)
    {
        std::cerr << "No C1' atom(s) found in base pair." << std::endl;
        return false;
    }
    else
    {
        float x, y, z;

        atomId1 = residueInfo[tempId1].start + atomId1;
        atomId2 = residueInfo[tempId2].start + atomId2;
        //std::cout << atomId1 << " " << atomId2 << " " << residueInfo[tempId1].start << " " << residueInfo[tempId2].start << std::endl;
        // Get atom coordinates when Hbond exists and the calculate the distance
        x =  s.coordinates()[atomId1][0] - s.coordinates()[atomId2][0];
        y =  s.coordinates()[atomId1][1] - s.coordinates()[atomId2][1];
        z =  s.coordinates()[atomId1][2] - s.coordinates()[atomId2][2];

        //std::cout << x << " " << y << " " << z << std::endl;
        // Calculated the distance and store it in vector distances
        float currentDist = sqrt(x*x+y*y+z*z);
        //      distances.push_back(currentDist);
        if (currentDist < 1.2)
        {
            std::cout << residueInfo[tempId1].type << " " << residueInfo[tempId1].resNr << " " << residueInfo[tempId2].type << " " << residueInfo[tempId2].resNr << " "<< currentDist << std::endl;
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

    //std::cout << "HBondD "  << " " << residueInfo[tempId1].type << " " << residueInfo[tempId2].type << std::endl;

    std::vector<float> distances;
    float              currentDist = 0;
    for  (int con = 0; con < this->getNrBondAtoms(); con++)
    {
        char * name1 = this->getBondIndex(con, 0);
        char * name2 = this->getBondIndex(con, 1);

        int    atomId1 = -1;
        int    atomId2 = -1;
        for (int x = 0; x < residueInfo[tempId1].size; x++)
        {

            //std::cout << residueInfo[tempId1].name[x] << std::endl;
            if (strcmp(name1, residueInfo[tempId1].name[x]) == 0)
            {
                atomId1 = x;
            }
        }
        for (int x = 0; x < residueInfo[tempId2].size; x++)
        {
            //std::cout << x << " " << name2 << " " << residueInfo[tempId2].name[x] << " ";

            if (strcmp(name2, residueInfo[tempId2].name[x]) == 0)
            {
                atomId2 = x;
            }
        }

        //std::cout << con << " " << this->getNrBondAtoms() << " "  << name1 << " " << name2 << " " << atomId1 << " " << atomId2 << " " << residueInfo[tempId1].size << " " << residueInfo[tempId2].size << std::endl;
        /*
           if (atomId1 == -1 || atomId2 == -1)
           {
           //	  std::cout << name1 << " " << name2 << std::endl;
            std::cerr << "No matching bond between " << name1 << " and " << name2 << std::endl;
            return false;
           }
         */

        //std::cout << name1 << " " << atomId1 << " "  << name2 << " " << atomId2 << std::endl;

        //if (atomId1 != -1 || atomId2 != -1)
        if (atomId1 != -1 && atomId2 != -1)
        {
            float x, y, z;

            atomId1 = residueInfo[tempId1].start + atomId1;
            atomId2 = residueInfo[tempId2].start + atomId2;
            // Get atom coordinates when Hbond exists and the calculate the distance
            x =  s.coordinates()[atomId1][0] - s.coordinates()[atomId2][0];
            y =  s.coordinates()[atomId1][1] - s.coordinates()[atomId2][1];
            z =  s.coordinates()[atomId1][2] - s.coordinates()[atomId2][2];

            // Calculated the distance and store it in vector distances
            currentDist = sqrt(x*x+y*y+z*z);
            distances.push_back(currentDist);
            //std::cout << currentDist << std::endl;

        }
    }           // end for

    // sort distances in vector from smallest to largest distance
    float* first(&distances[0]);
    float* last(first + distances.size());
    std::sort(first, last);
    /*
       if ((distances[0] <= dist) && (distances.size() < 3))
       {
        return true;
       }
       else if ((distances[0] <= dist) && (distances[1] <= dist) && (distances.size() == 3))
       {
        return true;
       }
     */
    if (distances[0] <= dist)
    {
        //std::cout << distances[0] << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}
/*
   float
   BasePair::calculateRMSD(int tempId1, int tempId2, Selection s)
   {
   float rmsd = 0.0;
   for (int i = 0; i < this->nrAtoms; i++)

    {

      std::cout << atoms->atomname[i] << std::endl;
    }

   //for  (int i = 0; i < this->getNrAtoms(); i++)
   //{
      //      std::cout << this->type << " " << this->resNr << " " << this->name[i] << std::endl;
      //std::cout << name[i] << " " << *this->type << " " << i << std::endl;
      //std::cout << i << " " << residueInfo[this].name[i] << std::endl;

   //}
   return rmsd;
   }
 */

// Set the isomerism type
void
BasePair::setIsomerism (const char * name, int size)
{
    // if "cis" is in name
    for (int i = 0; i < size - 2; i++)
    {
        if (name[i] == 'c' && name[i+1] == 'i' && name[i+2] == 's')
        {
            iso = Cis;
            return;
        }
    }
    iso = Trans;
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
