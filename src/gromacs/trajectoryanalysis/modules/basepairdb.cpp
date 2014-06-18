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
 * \author Anders Gärdenäs <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "basepairdb.h"

namespace gmx
{

namespace analysismodules
{

int resInfo::searchAtom(const char *name) const
{
    for (unsigned int x = 0; x < nAtoms(); x++)
    {
        if (atomName(x).compare(name) == 0)
        {
            return static_cast<int>(x);
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
real BasePair::computeRootMeanSquareDeviation(rvec * vec, int baseSize)
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
void BasePair::addAtom(rvec x, int atomnum, std::string atomname, real m)
{
    if (checkAtom(atomname))
    {
        atomNumber_.push_back(atomnum);
        atomName_.push_back(atomname);

        // convert coordinates from Ångström to nanometer
        // set the new atom position and mass
        for (unsigned int mm = 0; (mm < DIM); mm++)
        {
            atomPos_.push_back(x[mm]*0.1);
        }
        atomMass_.push_back(m);
    }
}


// Check whether the current atom is valid or not
bool BasePair::checkAtom(const std::string &name)
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
real *BasePair::atomMass()
{
    return &atomMass_[0];
}

// Set one letter code to the full name of the base
void BasePair::setNucleotides(const std::string &bases)
{
    for (int i = 0; i < 2; i++)
    {
        switch (bases[i])
        {
            case 'A':
                type[i] =  Adenine;
                break;
            case 'U':
                type[i] =  Uracil;
                break;
            case 'T':
                type[i] =  Thymine;
                break;
            case 'G':
                type[i] =  Guanine;
                break;
            case 'C':
                type[i] =  Cytosine;
                break;
            default:
                char buf[256];
                snprintf(buf, sizeof(buf), "invalid base type %c\n", bases[i]);
                GMX_THROW(InvalidInputError(buf));
        }
    }
}

// Return single letter code for each base type (either the 5' or 3')
char BasePair::nucleotideType (unsigned int i)
{
    GMX_ASSERT(i < 2, "nucleotideType out of range");
    switch (type[i])
    {
        case Adenine:
            return 'A';
        case Uracil:
            return 'U';
        case Thymine:
            return 'T';
        case Guanine:
            return 'G';
        case Cytosine:
            return 'C';
        default:
            GMX_THROW(APIError("NucleotideType incorrect"));
    }
    return 'X';
}

// Returns the bond index nr bond from pair
const char *BasePair::bondIndex(int bond, int pair)
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
void BasePair::setAtomDist(int offsetAtom)
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
int BasePair::bondType()
{
    int temp = 0;

    //sum formula start at 5
    temp = bondSize * bondtype[0];
    return (temp + bondtype[1]);
}

// Get the name of bonds
std::string BasePair::bondName(unsigned int type)
{
    GMX_ASSERT((type < 25), "invalid bond type number");
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
    }
    return "";
}

// Get the type of the bond
BondType BasePair::searchBondType(std::string name)
{
    const char *bt[] = { "Ww", "Ss", "Hh", "Bh", "Bs" };

    for (unsigned int f = 0; (f < bondSize); f++)
    {
        if (name.compare(bt[f]) == 0)
        {
            return static_cast<BondType>(f);
        }
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "invalid bond type %s", name.c_str());
    GMX_THROW(InvalidInputError(buf));
}

// Set bond type
void BasePair::setBondType (const std::string &bondType)
{
    GMX_ASSERT((bondType.size() == 4),
               "bondType string should be length 4 characters. Check your RNA database");
    bondtype[0] = searchBondType(bondType.substr(0, 2));
    bondtype[1] = searchBondType(bondType.substr(2, 2));
}


bool BasePair::isC1Contact(const resInfo *tempId1,
                           const resInfo *tempId2,
                           Selection      s)
{
    int atomId1 = -1;
    int atomId2 = -1;

    for (unsigned int x = 0; x < tempId1->nAtoms(); x++)
    {
        if (tempId1->atomName(x).compare("C1'") == 0)
        {
            atomId1 = x;
        }
    }
    for (unsigned int x = 0; x < tempId2->nAtoms(); x++)
    {
        if (tempId2->atomName(x).compare("C1'") == 0)
        {
            atomId2 = x;
        }
    }

    GMX_ASSERT((atomId1 != -1 && atomId2 != -1),
               "No C1' atom(s) found in base pair.");

    atomId1 = tempId1->atomStart() + atomId1;
    atomId2 = tempId2->atomStart() + atomId2;
    // Get atom coordinates when Hbond exists and the calculate the distance
    real currentDist = sqrt(distance2(s.coordinates()[atomId1],
                                      s.coordinates()[atomId2]));

    //TODO: what should it be down here?
    return (currentDist < 1.2 /* templateC1Dist */);
}

bool BasePair::isHBond(const resInfo *tempId1,
                       const resInfo *tempId2,
                       const t_pbc   *pbc,
                       real           maxDist,
                       Selection      s)
{
    std::vector<real> distances;
    bool              isHB = false;

    for (int con = 0; con < getNrBondAtoms(); con++)
    {
        int atomId1 = tempId1->searchAtom(bondIndex(con, 0));
        int atomId2 = tempId2->searchAtom(bondIndex(con, 1));

        if (atomId1 != -1 && atomId2 != -1)
        {
            atomId1 = tempId1->atomStart() + atomId1;
            atomId2 = tempId2->atomStart() + atomId2;
            // Calculate the distance and store it in vector distances
            ConstArrayRef<rvec> x = s.coordinates();
            rvec                dx;
            pbc_dx(pbc, x[atomId1], x[atomId2], dx);
            real                currentDist = norm(dx);
            distances.push_back(currentDist);
        }
    }

    if (distances.size() > 0)
    {
        // sort distances in vector from smallest to largest distance
        std::vector<real>::iterator di = std::min_element(distances.begin(), distances.end());
        isHB = (*di <= maxDist);
    }
    return isHB;
}

// Set the isomerism type
void BasePair::setIsomerism (std::string name)
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

} // namespace analysismodules

} // namespace gmx
