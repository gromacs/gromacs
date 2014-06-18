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
 * Declares gmx::analysismodules::BasePair.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders Gärdenäs <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_BASEPAIR_H
#define GMX_TRAJECTORYANALYSIS_BASEPAIR_H

#include "gmxpre.h"

#include "gromacs/math/vec.h"
#include "gromacs/selection/selection.h"

namespace gmx
{

namespace analysismodules
{

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

//! Information about residues
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


/*! \libinternal \brief
 * Utility class to store pairs of atom names
 */
class bondAtom
{
    public:
        bondAtom(std::string ai, std::string aj) : ai_(ai), aj_(aj) {}

        ~bondAtom() {}

        std::string atomI() { return ai_; }

        std::string atomJ() { return aj_; }
    private:
        std::string ai_, aj_;
};

/*! \libinternal \brief
 * Utility class to store and operate on a nucleotide base pair
 */
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
        int nrAtoms() { return atomName_.size(); }

        // Add an atom
        void addAtom(rvec x, int atomnum, std::string atomname, real m);

        // Return the atom positions
        rvec *atomPosition();

        // Add information about a pair of bonded atoms
        void addBondAtom(int ai, int aj);

        // Return true if it is the right template
        bool checkTemplate(char mol0, char mol1)
        {
            return (((nucleotideType(0) == mol0) && (nucleotideType(1) == mol1)) ||
                    ((nucleotideType(0) == mol1) && (nucleotideType(1) == mol0)));
        }

        // Calculate the RMSD
        real computeRootMeanSquareDeviation(rvec * vec, int baseSize);

        // Return the maximun distance in nm TODO: between what?
        real maximumDistance() { return dist; }

        // Set the maximum distance to another atom
        void setAtomDist(int index);

        // Return phosphate
        bool phosphateRmsd() { return phosphateRmsd_; }

        // Check if an atom is valid
        bool checkAtom (const std::string &name);

        // Return the colour with the chosen bond type representative
        //t_rgb getColourCode();

        // Return nucleotide type. Index should be 0 or 1 or the routine will throw.
        char nucleotideType (unsigned int index);

        // Set bond type
        void setBondtype (const char * name, int size);

        // Get bond type
        int getBondType();

        // Get bond name
        std::string getBondName(int type);

        // Check whether a bond exist by calculating the H bond distance
        bool getHBondDistance(const resInfo &tempID1,
                              const resInfo &tempID2,
                              real           bondDist,
                              Selection      sel);

        // Check whether the C1' distance exceed the template threshold
        bool getC1distance(const resInfo &tempID1,
                           const resInfo &tempID2,
                           Selection      sel);

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

} // namespace analysismodules

} // namespace gmx

#endif
