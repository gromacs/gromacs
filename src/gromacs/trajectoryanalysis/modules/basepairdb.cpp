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
 * Implements gmx::analysismodules::BasePair.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders G�rden�s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "basepairdb.h"

#include <cstdlib>

#include <algorithm>

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

namespace gmx
{

namespace analysismodules
{

const char *isomerismName(Isomerism iso)
{
    static const char *isomers[1 + IsomerismSize] = {"cis", "trans"};
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

char nucleotideToChar(Nucleotide nuc)
{
    switch (nuc)
    {
        case Adenine:
            return 'A';
        case Cytosine:
            return 'C';
        case Guanine:
            return 'G';
        case Thymine:
            return 'T';
        case Uracil:
            return 'U';
        default:
            break;
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "No charactor for nucleotide type %d", (int) nuc);
    GMX_THROW(APIError(buf));
}

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
    for (size_t f = 0; (f <= BondTypeSize); f++)
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


size_t ResInfo::searchAtom(const std::string &name) const
{
    for (size_t x = 0; x < nAtoms(); x++)
    {
        if (atomName(x).compare(name) == 0)
        {
            return x;
        }
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "Can not find atom %s in ResInfo for %c%zu", name.c_str(),
             residueType(), residueNumber());
    GMX_THROW(APIError(buf));
}

// Constructor
BasePair::BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd) : iso(IsomerismSize),
                                                                            hydrogenRmsd_(hydrogenRmsd),
                                                                            sugarRmsd_(sugarRmsd),
                                                                            phosphateRmsd_(phosphateRmsd)
{
    dist           = 0;
    type[0]        = NucleotideSize;
    type[1]        = NucleotideSize;
    bondtype_[0]   = BondTypeSize;
    bondtype_[1]   = BondTypeSize;
    resId[0]       = 0;
    resId[1]       = 0;
    templateRMSD   = 0;
    x_             = NULL;
    ePBC_          = 0;
    splitPoint_    = 0;
    clear_mat(box_);
}

// Destructor
BasePair::~BasePair()
{
    sfree(x_);
    for (auto &ba : bondAtoms)
    {
        delete ba;
    }
}

// Calculate the RMSD
real BasePair::computeRootMeanSquareDeviation(const ResInfo &base1,
                                              const ResInfo &base2,
                                              rvec *vec, size_t baseSize)
{
    // TODO: add check for matching sizes of arrays
//    char buf[256];
//    snprintf(buf, sizeof(buf), "Size mismatch, got %u coordinates, expected %u for %s",
//             baseSize, nrAtoms(), templateName().c_str());
//    GMX_RELEASE_ASSERT((baseSize == nrAtoms()), buf);
    bool   temp = false;
    real  *mass;
    rvec  *x;
    size_t nr = nrAtoms();
    size_t np = 0;
    if (baseSize == nr)
    {
        mass = &atomMass_[0];
        x    = x_;
    }
    else
    {
        temp = true;
        snew(mass, baseSize);
        snew(x, baseSize);
        if (phosphateRmsd_ && base1.endDirectionality() == FivePrimeEnd)
        {
            std::copy(&atomMass_[3], &atomMass_[splitPoint_], mass);
            copy_rvecn(x_ + 3, x, 0, splitPoint_ - 3);
            np += splitPoint_ - 3;
        }
        else
        {
            std::copy(&atomMass_[0], &atomMass_[splitPoint_], mass);
            copy_rvecn(x_, x, 0, splitPoint_);
            np += splitPoint_;
        }
        if (phosphateRmsd_ && base2.endDirectionality() == FivePrimeEnd)
        {
            std::copy(&atomMass_[splitPoint_ + 3], &atomMass_[nr], mass + np);
            copy_rvecn(x_ + splitPoint_ + 3, x + np, 0, nr - splitPoint_ - 3);
            np += nr - splitPoint_ - 3;
        }
        else
        {
            std::copy(&atomMass_[splitPoint_], &atomMass_[nr], mass + np);
            copy_rvecn(x_ + splitPoint_, x + np, 0, nr - splitPoint_);
            np += nr - splitPoint_;
        }
        GMX_RELEASE_ASSERT(np == baseSize, "Invalid structure used in calculating RMSD.");
    }

    // reset the coordinates to the origin if needed
    reset_x_ndim(3, baseSize, NULL, baseSize, NULL, vec, mass);
    // Set the matrix
    do_fit(baseSize, mass, x, vec);
    // Get the rmsd
    real rmsd = rmsdev(baseSize, mass, x, vec);
    if (temp)
    {
        sfree(mass);
        sfree(x);
    }
    return rmsd;
}

std::string BasePair::templateName() const
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%c%c %s %s",
             nucleotideTypeChar(0), nucleotideTypeChar(1),
             bondName().c_str(), isomerismName(isomerismType()));
    return std::string(buf);
}

void BasePair::addBondAtoms(size_t ai, size_t aj)
{
    std::string  ai_name, aj_name;
    size_t       ri = 0, rj = 0;
    // Extract the atom names corresponding to the indices
    for (size_t i = 0; i < atomNumber_.size(); i++)
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
    if (ai_name.size() == 0 || aj_name.size() == 0)
    {
        return;
    }
//    GMX_RELEASE_ASSERT(((ai_name.size() > 0) && (aj_name.size() > 0)),
//                       "Could not find atoms in CONECT");
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
                               size_t             atomnum,
                               const std::string &atomname,
                               const std::string &resname,
                               size_t             resnum,
                               real               m)
{
    if (atomName_.size() == 0)
    {
        // Compare residue type
        std::string nt = "";
        nt += nucleotideTypeChar(0);
        char        buf[256];
        snprintf(buf, sizeof(buf), "Incorrect order of residues for template %s",
                 templateName().c_str());
        GMX_RELEASE_ASSERT((resname.compare(nt) == 0), buf);
    }
    if (checkAtom(atomname))
    {
        atomNumber_.push_back(atomnum);
        atomName_.push_back(atomname);
        if (resNumber_.size() > 0 && resnum != resNumber_.back())
        {
            splitPoint_ = resNumber_.size();
        }
        resNumber_.push_back(resnum);

        // set the new atom position and mass
        srenew(x_, nrAtoms());
        copy_rvec(x, x_[nrAtoms() - 1]);
        atomMass_.push_back(m);
    }
}


// Check whether the current atom is valid or not
bool BasePair::checkAtom(const std::string &name)
{
    bool set = false;
    for (size_t i = 0; i < name.size(); ++i)
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

BasePairComparison BasePair::checkTemplate(Nucleotide nuc0, Nucleotide nuc1)
{
    bool b1 = (type[0] == nuc0) && (type[1] == nuc1);
    bool b2 = (type[0] == nuc1) && (type[1] == nuc0);
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
        type[i] = searchNucleotide(bases[i]);
    }
}

// Return base type (either the 5' or 3')
Nucleotide BasePair::nucleotideType(size_t i) const
{
    GMX_RELEASE_ASSERT(i < 2, "nucleotideType out of range");
    return type[i];
}

char BasePair::nucleotideTypeChar(size_t index) const
{
    return nucleotideToChar(nucleotideType(index));
}

const char *BasePair::bondIndex(size_t bond, int pair) const
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
void BasePair::setAtomDist(size_t offsetAtom)
{
    real                temp;
    ConstArrayRef<rvec> x = atomPosition();
    dist = 0;
    for (size_t i = 0; i < atomName_.size(); i++)
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
std::string BasePair::bondName(bool swap) const
{
    std::string bn;
    if (swap)
    {
        bn.append(bondTypeName(bondtype_[1]));
        bn.append(bondTypeName(bondtype_[0]));
    }
    else
    {
        bn.append(bondTypeName(bondtype_[0]));
        bn.append(bondTypeName(bondtype_[1]));
    }
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
    size_t atomId1 = tempId1->searchAtom("C1'") + tempId1->atomStart();
    size_t atomId2 = tempId2->searchAtom("C1'") + tempId2->atomStart();

    // Get atom coordinates when Hbond exists and the calculate the distance
    real currentDist = sqrt(distance2(s.coordinates()[atomId1],
                                      s.coordinates()[atomId2]));

    //TODO: what should it be down here?
    return (currentDist < 1.2 /* templateC1Dist */);
}

bool BasePair::checkHydrogen(const ResInfo       &base1,
                             const ResInfo       &base2,
                             const std::string   &b_1,
                             const std::string   &b_2,
                             std::vector<size_t> &hs)
{
    hs.clear();
    std::string p1("H");
    if (b_1.compare("O2'") == 0)
    {
        p1 += "O'2";
    }
    else
    {
        p1 += (b_1.c_str() + 1);
    }
    for (size_t i = 0; i < base1.nAtoms(); ++i)
    {
        if (base1.atomName(i).compare(p1) == 0 ||
            base1.atomName(i).compare(p1 + "1") == 0 ||
            base1.atomName(i).compare(p1 + "2") == 0)
        {
            hs.push_back(i + base1.atomStart());
        }
    }
    std::string p2("H");
    if (b_2.compare("O2'") == 0)
    {
        p2 += "O'2";
    }
    else
    {
        p2 += (b_2.c_str() + 1);
    }
    for (size_t i = 0; i < base2.nAtoms(); ++i)
    {
        if (base2.atomName(i).compare(p2) == 0 ||
            base2.atomName(i).compare(p2 + "1") == 0 ||
            base2.atomName(i).compare(p2 + "2") == 0)
        {
            hs.push_back(i + base2.atomStart());
        }
    }
    return hs.size() > 0;
}

bool BasePair::isHBond(const ResInfo   &base1,
                       const ResInfo   &base2,
                       const t_pbc     *pbc,
                       double           maxDist,
                       const Selection &s)
{
    std::vector<real> distances;
    size_t            va = 0, ia = 0;

    for (std::vector<BondAtom *>::iterator bai = bondAtoms.begin();
         (bai < bondAtoms.end()); ++bai)
    {
        const std::string &b_1     = (*bai)->atomI();
        const std::string &b_2     = (*bai)->atomJ();
        size_t             atomId1 = base1.searchAtom(b_1) + base1.atomStart();
        size_t             atomId2 = base2.searchAtom(b_2) + base2.atomStart();

        // Calculate the distance and store it in vector distances
        ConstArrayRef<rvec> x = s.coordinates();

        // TODO: check hydrogen atom.
        std::vector<size_t> hs;
        if (!checkHydrogen(base1, base2, b_1, b_2, hs))
        {
            return false;
        }

        bool   has_c = false;
        size_t nc, no;
        nc = no = x.size();
        if (b_1.find("C") < b_1.npos)
        {
            nc    = atomId1;
            no    = atomId2;
            has_c = true;
        }
        else if (b_2.find("C") < b_2.npos)
        {
            nc    = atomId2;
            no    = atomId1;
            has_c = true;
        }
        bool valid_angle = false;
        for (const size_t &h : hs)
        {
            rvec dx1, dx2;
            if (has_c)
            {
                pbc_dx(pbc, x[no], x[h], dx1);
                pbc_dx(pbc, x[no], x[nc], dx2);
                if (norm(dx1) <= 0.27 && norm(dx2) <= 0.37)
                {
                    valid_angle = true;
                    break;
                }
            }
            pbc_dx(pbc, x[atomId1], x[h], dx1);
            pbc_dx(pbc, x[atomId2], x[h], dx2);
            real angle = gmx_angle(dx1, dx2);
            if (angle <= -min_angle || angle >= min_angle)
            {
                valid_angle = true;
                break;
            }
        }
        // if (!valid_angle)
        // TODO: Don't ignore C!
        if (!valid_angle && !(has_c && bondAtoms.size() > 1))
        {
            ia++;
            continue;
//            return false;
        }
        else
        {
            va++;
        }
        rvec dx;
        pbc_dx(pbc, x[atomId1], x[atomId2], dx);


        real currentDist = norm(dx);
        if (NULL != debug)
        {
            fprintf(debug, "atomId1 %c%zu %zu %s atomId2 %c%zu %zu %s dist %g\n",
                    base1.residueType(), base1.residueNumber(), atomId1,
                    base1.atomName(atomId1 - base1.atomStart()).c_str(),
                    base2.residueType(), base2.residueNumber(), atomId2,
                    base2.atomName(atomId2 - base2.atomStart()).c_str(),
                    currentDist);
        }
        distances.push_back(currentDist);
    }

    if (distances.size() > 0)
    {
        std::vector<real>::iterator di = std::min_element(distances.begin(), distances.end());
        return *di <= maxDist && (ia == 0 || va >= 2);
    }
    return false;
}

// Set the isomerism type
void BasePair::setIsomerism(std::string name)
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
