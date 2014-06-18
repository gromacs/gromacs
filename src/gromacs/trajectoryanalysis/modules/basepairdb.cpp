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

#include <algorithm>

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

#include "basepairdb.h"

namespace gmx
{

namespace analysismodules
{

const char *isomerismName(Isomerism iso)
{
    static const char *isomers[1+IsomerismSize] = { "cis", "trans" };
    switch (iso)
    {
        case Cis:
            return isomers[0];
        case Trans:
            return isomers[1];
    }
    GMX_THROW(APIError("Unknown Isomerism"));
}

Nucleotide searchNucleotide(char c)
{
    switch (c)
    {
        case 'A':
            return Adenine;
        case 'C':
            return Cytosine;
        case 'G':
            return Guanine;
        case 'T':
            return Thymine;
        case 'U':
            return Uracil;
        default:
            break;
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "No such nucleotide type '%c'", c);
    GMX_THROW(APIError(buf));
}

static const char *btName[1+BondTypeSize] =  { "Ww", "Ss", "Hh", "Bh", "Bs" };
const char *bondTypeName(BondType bt)
{
    switch (bt)
    {
        case Ww:
            return btName[0];
        case Ss:
            return btName[1];
        case Hh:
            return btName[2];
        case Bh:
            return btName[3];
        case Bs:
            return btName[4];
    }
    GMX_THROW(APIError("Unknown BondType"));
}

BondType searchBondType(std::string name)
{
    for (unsigned int f = 0; (f <= BondTypeSize); f++)
    {
        if (name.compare(btName[f]) == 0)
        {
            return static_cast<BondType>(f);
        }
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "invalid bond type %s", name.c_str());
    GMX_THROW(InvalidInputError(buf));
}


unsigned int ResInfo::searchAtom(const char *name) const
{
    for (unsigned int x = 0; x < nAtoms(); x++)
    {
        if (atomName(x).compare(name) == 0)
        {
            return x;
        }
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "Can not find atom %s in ResInfo for %c%d", name,
             residueType(), residueNumber());
    GMX_THROW(APIError(buf));
}

// Constructor
BasePair::BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd)
{
    dist              = 0;
    //type[0]           = Nucnull;
    //type[1]           = Nucnull;
    //bondtype_[0]       = Bondnull;
    //bondtype_[1]       = Bondnull;
    //iso               = Isonull;
    resId[0]          = 0;
    resId[1]          = 0;
    templateRMSD      = 0;
    hydrogenRmsd_     = hydrogenRmsd;
    sugarRmsd_        = sugarRmsd;
    phosphateRmsd_    = phosphateRmsd;
    x_                = NULL;
    ePBC_             = 0;
    clear_mat(box_);
}

// Destructor
BasePair::~BasePair()
{
    sfree(x_);
}

// Calculate the RMSD
real BasePair::computeRootMeanSquareDeviation(rvec * vec, unsigned int baseSize)
{
    // TODO: add check for matching sizes of arrays
    char buf[256];
    snprintf(buf, sizeof(buf), "Size mismatch, got %u coordinates, expected %u for %s",
             baseSize, nrAtoms(), templateName().c_str());
    GMX_RELEASE_ASSERT((baseSize == nrAtoms()), buf);
    // Set the matrix
    do_fit(baseSize, atomMass(), x_, vec);
    // Get the rmsd
    return rmsdev(baseSize, atomMass(), x_, vec);
}

std::string BasePair::templateName()
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%c%c %s %s",
             nucleotideType(0), nucleotideType(1),
             bondName().c_str(), isomerismName(isomerismType()));
    return std::string(buf);
}

void BasePair::addBondAtoms(unsigned int ai, unsigned int aj)
{
    std::string  ai_name, aj_name;
    unsigned int ri = 0, rj = 0;
    // Extract the atom names corresponding to the indices
    for (unsigned int i = 0; i < atomNumber_.size(); i++)
    {
        if (ai == atomNumber_[i])
        {
            ai_name = atomName_[i];
            ri      = resNumber_[i];
        }
        if (aj == atomNumber_[i])
        {
            aj_name = atomName_[i];
            rj      = resNumber_[i];
        }
    }
    GMX_RELEASE_ASSERT(((ai_name.size() > 0) && (aj_name.size() > 0)),
                       "Could not find atoms in CONECT");
    GMX_RELEASE_ASSERT((ri != rj),
                       "The atoms in the CONECT record in the template database are in the same nucleotide");

    if (NULL != debug)
    {
        fprintf(debug, "Going to push a bond %s %s for %s\n",
                ai_name.c_str(), aj_name.c_str(),
                templateName().c_str());
    }
    BondAtom *ba = new BondAtom(ai_name, aj_name);
    bondAtoms.push_back(ba);
}

// Add one atom.
void BasePair::addBasePairAtom(rvec               x,
                               unsigned int       atomnum,
                               const std::string &atomname,
                               const std::string &resname,
                               unsigned int       resnum,
                               real               m)
{
    if (atomName_.size() == 0)
    {
        // Compare residue type
        std::string nt = "";
        nt += nucleotideType(0);
        char        buf[256];
        snprintf(buf, sizeof(buf), "Incorrect order of residues for template %s",
                 templateName().c_str());
        GMX_RELEASE_ASSERT((resname.compare(nt) == 0), buf);
    }
    if (checkAtom(atomname))
    {
        atomNumber_.push_back(atomnum);
        atomName_.push_back(atomname);
        resNumber_.push_back(resnum);

        // set the new atom position and mass
        srenew(x_, nrAtoms());
        copy_rvec(x, x_[nrAtoms()-1]);
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

BasePairComparison BasePair::checkTemplate(char nuc0, char nuc1)
{
    bool b1 = (nucleotideType(0) == nuc0) && (nucleotideType(1) == nuc1);
    bool b2 = (nucleotideType(0) == nuc1) && (nucleotideType(1) == nuc0);
    if (b1 && b2)
    {
        return BasePairComparison_Both;
    }
    else if (b1)
    {
        return BasePairComparison_Match;
    }
    else if (b2)
    {
        return BasePairComparison_Swap;
    }
    else
    {
        return BasePairComparison_MisMatch;
    }
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
    GMX_RELEASE_ASSERT(i < 2, "nucleotideType out of range");
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

const char *BasePair::bondIndex(unsigned int bond, int pair) const
{
    GMX_RELEASE_ASSERT((bond < nrBondAtoms()), "Bond index out of range");
    if (pair == 0)
    {
        return bondAtoms[bond]->atomI().c_str();
    }
    else
    {
        return bondAtoms[bond]->atomJ().c_str();
    }
}

// Sets the maximum distance to another atom
void BasePair::setAtomDist(int offsetAtom)
{
    real                temp;
    ConstArrayRef<rvec> x = atomPosition();
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

// Get the name of bonds
std::string BasePair::bondName()
{
    std::string bn(bondTypeName(bondtype_[0]));
    bn.append(bondTypeName(bondtype_[1]));

    return bn;
}

// Set bond type
void BasePair::setBondTypes(const std::string &bondTypes)
{
    GMX_RELEASE_ASSERT((bondTypes.size() == 4),
                       "bondTypes string should be length 4 characters. Check your RNA database");
    bondtype_[0] = searchBondType(bondTypes.substr(0, 2));
    bondtype_[1] = searchBondType(bondTypes.substr(2, 2));
}


bool BasePair::isC1Contact(const ResInfo *tempId1,
                           const ResInfo *tempId2,
                           Selection      s)
{
    unsigned int atomId1 = tempId1->searchAtom("C1'") + tempId1->atomStart();
    unsigned int atomId2 = tempId2->searchAtom("C1'") + tempId2->atomStart();

    // Get atom coordinates when Hbond exists and the calculate the distance
    real currentDist = sqrt(distance2(s.coordinates()[atomId1],
                                      s.coordinates()[atomId2]));

    //TODO: what should it be down here?
    return (currentDist < 1.2 /* templateC1Dist */);
}

bool BasePair::isHBond(const ResInfo *tempId1,
                       const ResInfo *tempId2,
                       const t_pbc   *pbc,
                       real           maxDist,
                       Selection      s)
{
    std::vector<real> distances;
    bool              isHB = false;

    for (std::vector<BondAtom *>::iterator bai = bondAtoms.begin();
         (bai < bondAtoms.end()); ++bai)
    {
        const char  *b_1     = (*bai)->atomI().c_str();
        const char  *b_2     = (*bai)->atomJ().c_str();
        unsigned int atomId1 = tempId1->searchAtom(b_1) + tempId1->atomStart();
        unsigned int atomId2 = tempId2->searchAtom(b_2) + tempId2->atomStart();

        // Calculate the distance and store it in vector distances
        ConstArrayRef<rvec> x = s.coordinates();
        rvec                dx;
        pbc_dx(pbc, x[atomId1], x[atomId2], dx);
        real                currentDist = norm(dx);
        if (NULL != debug)
        {
            fprintf(debug, "atomId1 %c%d %u %s atomId2 %c%d %u %s dist %g\n",
                    tempId1->residueType(), tempId1->residueNumber(), atomId1,
                    tempId1->atomName(atomId1-tempId1->atomStart()).c_str(),
                    tempId2->residueType(), tempId2->residueNumber(), atomId2,
                    tempId2->atomName(atomId2-tempId2->atomStart()).c_str(),
                    currentDist);
        }
        distances.push_back(currentDist);
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
