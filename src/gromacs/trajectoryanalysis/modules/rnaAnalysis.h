/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Declares trajectory analysis module for RNA analysis (PDB structures and trajectories).
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_RNAANALYSIS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_RNAANALYSIS_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <list>
#include <dirent.h>

#include "../analysismodule.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/matio.h"

#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/path.h"

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/legacyheaders/copyrite.h"
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
 * Class used to analyze RNA structures
 *
 * \ingroup module_trajectoryanalysis
 */

// Residue information
struct resInfo
{
    char   type;
    int    start;
    bool   invalid[50];
    char * name[50];
    int    resNr;
    int    size;

};
resInfo                          * residueInfo;

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
//virtual void finishAnalysis();

//! Routine to write output, that is additional over the built-in
        virtual void writeOutput();

    private:

        int                                findSameAtom(int i, char * name);
        Selection                          sel_;
        AnalysisData                       data_;
        AnalysisDataAverageModulePointer   adata_;

//resInfo                          * residueInfo;
        AnalysisNeighborhood               nb_;

        char                             * libfn;

// Copy and assign disallowed by base.
// A list list of a found pairs
        std::list < std::list < pairInfo>> found;
        std::list<frame>                    coords;

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
    Cytosine, Guanine, Adenine, Uracil, Thymine
};

enum Isomerism
{
    Cis, Trans
};

enum BondType
{
    Ww, Ss, Hh, Bh, Bs, bondSize
};


//  RNA base pairs
class BasePair
{
    public:

// Constructor
        BasePair()
        {
            conections     = gmx_conect_init();
            currentNrAtoms = 0;
            nrAtoms        = 0;
            dist           = 0;
            nrBondAtoms    = 0;
        }

// Destructor
        virtual ~BasePair();

// Set the type of nuclutide
        void setNucleotides(const char *  nuc);

// Set the number of atoms and initialize their mass
        void setNrAtoms(int i)
        {
            nrAtoms     = i;
            this->atoms = new rvec[i];
            this->mass  = new float[i];
        }

// Get the nummber of atoms
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
            //    return (getNucleotides(0) ==  mol0->type && getNucleotides(1) ==  mol1->type ) ||
            //         (getNucleotides(1) ==  mol0->type && getNucleotides(0) ==  mol1->type );
            //std::cout << getNucleotides(0) << " " << mol0->type << std::endl;
            if (getNucleotides(0) ==  mol0->type && getNucleotides(1) ==  mol1->type)
            {
                //std::cout << "true " << getNucleotides(0) << " " << mol0->type << " - " << getNucleotides(1) << " " << mol1->type << std::endl;
                return true;
            }
            else
            {
                //std::cout << "false " << getNucleotides(0) << " " << mol0->type << " - " << getNucleotides(1) << " " << mol1->type << std::endl;

                return false;
            }

            //return (getNucleotides(0) ==  mol0->type && getNucleotides(1) ==  mol1->type );


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

// Set if it should include hydrogen atoms
        static void setHydrogen(bool hydro)
        {
            hydrogen = hydro;
        }

// Set if it should include only base atoms
        static void setOnlyAmino(bool amino)
        {
            onlyAmino = amino;
        }

// Set if it should include phoshpate atoms
        static void setPhosphate(bool phos)
        {
            phosphate = phos;
        }

// Return phosphate
        static bool getPhosphate()
        {
            return phosphate;
        }

// Return conections
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

// Get the type of bonds
        int getBondType();

        bool getHBondDistance(int tempID1, int tempID2, float bondDist, Selection sel);
        bool getC1distance(int tempID1, int tempID2, Selection sel);
//float calculateRMSD(int tempId1, int tempId2, Selection s);


// Set the isomerism type
        void setIsomerism (const char * name, int size);

// Get the isomerism type
        int getIso();

// Get the mass
        float * getMass();

// Set the name of the atom bond (can be done multiple times)
        void setBondAtom(char * name1, char * name2);
//void setBondAtom(char name1, char name2);

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
        float     * mass;
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




}       // namespace

class RnaAnalysisInfo
{

    public:
        static const char name[];
        static const char shortDescription[];
        static TrajectoryAnalysisModulePointer create();


};

} // namespace analysismodules

} // namespace gmx

#endif
