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
 * Implements gmx::analysismodules::BasePair.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders G√§rden√§s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/math/do_fit.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "basepair.h"

namespace gmx
{

namespace analysismodules
{

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

rvec *BasePair::atomPosition()
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

void BasePair::addBondAtom(int ai, int aj)
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
BasePair::getC1distance(const resInfo &tempId1,
                        const resInfo &tempId2,
                        Selection      s)
{
    int atomId1 = -1;
    int atomId2 = -1;

    for (int x = 0; x < tempId1.size; x++)
    {
        if (strcmp(tempId1.name[x], "C1'") == 0)
        {
            atomId1 = x;
        }
    }
    for (int x = 0; x < tempId2.size; x++)
    {
        if (strcmp(tempId2.name[x], "C1'") == 0)
        {
            atomId2 = x;
        }
    }

    if (atomId1 == -1 || atomId2 == -1)
    {
        GMX_THROW(InvalidInputError("No C1' atom(s) found in base pair."));
    }
    else
    {
        real x, y, z;

        atomId1 = tempId1.start + atomId1;
        atomId2 = tempId2.start + atomId2;
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
BasePair::getHBondDistance(const resInfo &tempId1,
                           const resInfo &tempId2,
                           real           dist,
                           Selection      s)
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

        for (int x = 0; x < tempId1.size; x++)
        {
            if (strcmp(name1, tempId1.name[x]) == 0)
            {
                atomId1 = x;
                //std::cout << "Atom 1 " << name1 << " " << atomId1 << std::endl;
            }
        }
        for (int x = 0; x < tempId2.size; x++)
        {
            if (strcmp(name2, tempId2.name[x]) == 0)
            {
                atomId2 = x;
                //std::cout << "Atom 2 " << name2 << " " << atomId2 << std::endl;
            }

        }

        if (atomId1 != -1 && atomId2 != -1)
        {

            real x, y, z;
            atomId1 = tempId1.start + atomId1;
            atomId2 = tempId2.start + atomId2;
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

} // namespace analysismodules

} // namespace gmx
