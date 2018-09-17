/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// structure related stuff

#include "structure.h"

// #include "primitives-3d.h"

#include "gromacs/math/utilities.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <cmath>
#include <deque>
#include <memory>
#include <set>
#include <vector>
#include <numeric>
#include <functional>
#include <sstream>
// --------------------------------------------------------------------

const double
    kSSBridgeDistance     = 3.0,
    kMinimalDistance      = 0.5,
    kMinimalCADistance    = 9.0,
    kMinHBondEnergy       = -9.9,
    kMaxHBondEnergy       = -0.5,
    kCouplingConstant     = -27.888, //  = -332 * 0.42 * 0.2
    kMaxPeptideBondLength = 2.5;

const double
    kRadiusN        = 1.65,
    kRadiusCA       = 1.87,
    kRadiusC        = 1.76,
    kRadiusO        = 1.4,
    kRadiusSideAtom = 1.8,
    kRadiusWater    = 1.4;

// --------------------------------------------------------------------
enum class SecondaryStructure : char
{
    loop       = ' ',
    alphahelix = 'H',
    betabridge = 'B',
    strand     = 'E',
    helix_3    = 'G',
    helix_5    = 'I',
    turn       = 'T',
    bend       = 'S'
};

// 22 real letters and 1 dummy
const char    kResidues[]       = "ACDEFGHIKLMNPQRSTVWYBZX";
const uint8_t kResidueNrTable[] = {
//  A   B   C   D   E   F   G   H   I       K   L   M   N       P   Q   R   S   T  U=X  V   W   X   Y   Z
//  0,  1,  2,  3,  4,  5,  6,  7,  8, 23,  9, 10, 11, 12, 23, 13, 14, 15, 16, 17, 22, 18, 19, 22, 20, 21
    0, 20,  1,  2,  3,  4,  5,  6,  7, 23,  8,  9, 10, 11, 23, 12, 13, 14, 15, 16, 22, 17, 18, 22, 19, 21
};

sequence encode(const std::string &s)
{
    sequence result(s.length(), 0);
    for (unsigned int i = 0; i < s.length(); ++i)
    {
        result[i] = is_gap(s[i]) ? '-' : ResidueNr(s[i]);
    }
    return result;
}

std::string decode(const sequence &s)
{
    std::string result;
    for (const auto &ss : s)
    {
        result += ss >= 23 ? '.' : kResidues[ss];
    }
    return result;
}


namespace
{

double Distance(const gmx::RVec &p1, const gmx::RVec &p2)
{
    return (p1-p2).norm();
}
double DistanceSquared(const gmx::RVec &p1, const gmx::RVec &p2)
{
    return (p1-p2).norm2();
}


double DihedralAngle(const gmx::RVec &p1, const gmx::RVec &p2, const gmx::RVec &p3, const gmx::RVec &p4)
{
    gmx::RVec v12 = p1 - p2; // std::vector from p2 to p1
    gmx::RVec v43 = p4 - p3; // std::vector from p3 to p4

    gmx::RVec z = p2 - p3;   // std::vector from p3 to p2

    gmx::RVec p = z.cross(v12);
    gmx::RVec x = z.cross(v43);
    gmx::RVec y = z.cross(x);

    float     u = x.dot(x);
    float     v = y.dot(y);

    double    result = 360;
    if (u > 0 and v > 0)
    {
        u = p.dot(x) / sqrt(u);
        v = p.dot(y) / sqrt(v);
        if (u != 0 or v != 0)
        {
            result = atan2(v, u) * 180 / M_PI;
        }
    }

    return result;
}
double CosinusAngle(const gmx::RVec &p1, const gmx::RVec &p2, const gmx::RVec &p3, const gmx::RVec &p4)
{
    const gmx::RVec v12 = p1 - p2;
    const gmx::RVec v34 = p3 - p4;

    double          result = 0;

    float           x = v12.dot(v12) * v34.dot(v34);
    if (x > 0)
    {
        result = v12.dot(v34) / sqrt(x);
    }

    return result;
}

// we use a fibonacci spheres to calculate the even distribution of the dots
class MSurfaceDots
{
    public:
        static MSurfaceDots &Instance();

        size_t size() const { return mPoints.size(); }
        const gmx::RVec &operator[](uint32_t inIx) const  { return mPoints[inIx]; }
        double weight() const { return mWeight; }

    private:
        MSurfaceDots(int32_t inN);

        std::vector<gmx::RVec> mPoints;
        double                 mWeight;
};

MSurfaceDots &MSurfaceDots::Instance()
{
    const u_int32_t      kN = 200;

    static MSurfaceDots  sInstance(kN);
    return sInstance;
}

MSurfaceDots::MSurfaceDots(int32_t N)
{
    int32_t     P = 2 * N + 1;

    const float kGoldenRatio = (1 + std::sqrt(5.0f)) / 2;

    mWeight = (4 * M_PI) / P;

    for (int32_t i = -N; i <= N; ++i)
    {
        float     lat = std::asin((2.0f * i) / P);
        float     lon = fmod(i, kGoldenRatio) * 2 * M_PI / kGoldenRatio;

        gmx::RVec p;
        p[XX] = sin(lon) * cos(lat);
        p[YY] = cos(lon) * cos(lat);
        p[ZZ] = sin(lat);

        mPoints.push_back(p);
    }
}

}

// --------------------------------------------------------------------

namespace
{

void trim_left(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
                                        return !std::isspace(ch);
                                    }));
}

void trim_right(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
                             return !std::isspace(ch);
                         }).base(), s.end());
}

void trim(std::string &s)
{
    trim_left(s);
    trim_right(s);
}

std::string trim_copy(const std::string &s)
{
    std::string copyString {
        s
    };
    trim(copyString);
    return copyString;
}

}

MAtomType MapElement(std::string inElement)
{
    trim(inElement);
    std::transform(inElement.begin(), inElement.end(), inElement.begin(), ::toupper);

    MAtomType result = kUnknownAtom;
    if (inElement == "H")
    {
        result = kHydrogen;
    }
    else if (inElement == "C")
    {
        result = kCarbon;
    }
    else if (inElement == "N")
    {
        result = kNitrogen;
    }
    else if (inElement == "O")
    {
        result = kOxygen;
    }
    else if (inElement == "F")
    {
        result = kFluorine;
    }
    else if (inElement == "P")
    {
        result = kPhosphorus;
    }
    else if (inElement == "S")
    {
        result = kSulfur;
    }
    else if (inElement == "CL")
    {
        result = kChlorine;
    }
    else if (inElement == "K")
    {
        result = kPotassium;
    }
    else if (inElement == "MG")
    {
        result = kMagnesium;
    }
    else if (inElement == "CA")
    {
        result = kCalcium;
    }
    else if (inElement == "ZN")
    {
        result = kZinc;
    }
    else if (inElement == "SE")
    {
        result = kSelenium;
    }
    return result;
}

const std::array<MResidueInfo, 21> kResidueInfo = {
    {    { ResidueType::UnknownResidue, 'X', "UNK" },
         { ResidueType::Alanine,        'A', "ALA" },
         { ResidueType::Arginine,       'R', "ARG" },
         { ResidueType::Asparagine,     'N', "ASN" },
         { ResidueType::AsparticAcid,   'D', "ASP" },
         { ResidueType::Cysteine,       'C', "CYS" },
         { ResidueType::GlutamicAcid,   'E', "GLU" },
         { ResidueType::Glutamine,      'Q', "GLN" },
         { ResidueType::Glycine,        'G', "GLY" },
         { ResidueType::Histidine,      'H', "HIS" },
         { ResidueType::Isoleucine,     'I', "ILE" },
         { ResidueType::Leucine,        'L', "LEU" },
         { ResidueType::Lysine,         'K', "LYS" },
         { ResidueType::Methionine,     'M', "MET" },
         { ResidueType::Phenylalanine,  'F', "PHE" },
         { ResidueType::Proline,        'P', "PRO" },
         { ResidueType::Serine,         'S', "SER" },
         { ResidueType::Threonine,      'T', "THR" },
         { ResidueType::Tryptophan,     'W', "TRP" },
         { ResidueType::Tyrosine,       'Y', "TYR" },
         { ResidueType::Valine,         'V', "VAL" }}
};


ResidueType MapResidue(std::string inName)
{
    trim(inName);

    for (auto residueInfo : kResidueInfo)
    {
        if (inName == residueInfo.name)
        {
            return residueInfo.type;
        }
    }

    return ResidueType::UnknownResidue;
}

ResidueType MapResidue(char inCode)
{
    for (auto residueInfo : kResidueInfo)
    {
        if (inCode == residueInfo.code)
        {
            return residueInfo.type;
        }
    }

    return ResidueType::UnknownResidue;
}

// --------------------------------------------------------------------
// a custom float parser, optimised for speed (and the way floats are represented in a PDB file)

double ParseFloat(const std::string &s)
{
    double result = 0;
    bool   negate = false;
    double div    = 10;

    enum State {
        pStart, pSign, pFirst, pSecond
    } state = pStart;

    for (std::string::const_iterator ch = s.begin(); ch != s.end(); ++ch)
    {
        switch (state)
        {
            case pStart:
                if (isspace(*ch))
                {
                    continue;
                }
                if (*ch == '-')
                {
                    negate = true;
                    state  = pSign;
                }
                else if (*ch == '+')
                {
                    state = pSign;
                }
                else if (*ch == '.')
                {
                    state = pSecond;
                }
                else if (isdigit(*ch))
                {
                    result = *ch - '0';
                    state  = pFirst;
                }
                break;

            case pSign:
                if (*ch == '.')
                {
                    state = pSecond;
                }
                else if (isdigit(*ch))
                {
                    state  = pFirst;
                    result = *ch - '0';
                }
                break;

            case pFirst:
                if (*ch == '.')
                {
                    state = pSecond;
                }
                else if (isdigit(*ch))
                {
                    result = 10 * result + (*ch - '0');
                }
                break;

            case pSecond:
                if (isdigit(*ch))
                {
                    result += (*ch - '0') / div;
                    div    *= 10;
                }
                break;
        }
    }

    if (negate)
    {
        result = -result;
    }

    return result;
}


struct MBridge
{
    MBridgeType           type;
    u_int32_t             sheet, ladder;
    std::set<MBridge*>    link;
    std::deque<u_int32_t> i, j;
    std::string           chainI, chainJ;

    bool operator<(const MBridge &b) const
    {
        return chainI < b.chainI or (chainI == b.chainI and i.front() < b.i.front());
    }
};

std::ostream &operator<<(std::ostream &os, const MBridge &b)
{
    os << '[' << (b.type == btParallel ? "p" : "a") << ':' << b.i.front()
    << '-' << b.i.back() << '/' << b.j.front() << '-' << b.j.back() << ']';
    return os;
}

// return true if any of the residues in bridge a is identical to any of the
// residues in bridge b
bool Linked(const MBridge &a, const MBridge &b)
{
    return
        find_first_of(a.i.begin(), a.i.end(), b.i.begin(), b.i.end()) != a.i.end() or
            find_first_of(a.i.begin(), a.i.end(), b.j.begin(), b.j.end()) != a.i.end() or
            find_first_of(a.j.begin(), a.j.end(), b.i.begin(), b.i.end()) != a.j.end() or
            find_first_of(a.j.begin(), a.j.end(), b.j.begin(), b.j.end()) != a.j.end();
}

// --------------------------------------------------------------------

MResidue::MResidue(int32_t inNumber, MResidue* inPrevious, const std::vector<MAtom> &inAtoms)
    : mPrev(inPrevious),
      mNext(nullptr),
      mSeqNumber(inAtoms.front().mResSeq),
      mNumber(inNumber),
      mInsertionCode(inAtoms.front().mICode),
      mType(MapResidue(inAtoms.front().mResName)),
      mSSBridgeNr(0),
      mAccessibility(0),
      mSecondaryStructure(SecondaryStructure::loop),
      mSheet(0),
      mBend(false)
{
    if (mPrev != nullptr)
    {
        mPrev->mNext = this;
    }

    std::fill(mHelixFlags, mHelixFlags + 3, helixNone);

    mBetaPartner[0].residue = mBetaPartner[1].residue = nullptr;

    mHBondDonor[0].energy  = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
    mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nullptr;

    static const MAtom kNullAtom = {};
    mN = mCA = mC = mO = kNullAtom;

    for (const auto &atom : inAtoms)
    {
        if (mChainID.empty())
        {
            mChainID = atom.mChainID;
        }

        if (atom.GetName() == "N")
        {
            mN = atom;
        }
        else if (atom.GetName() == "CA")
        {
            mCA = atom;
        }
        else if (atom.GetName() == "C")
        {
            mC = atom;
        }
        else if (atom.GetName() == "O")
        {
            mO = atom;
        }
        else
        {
            mSideChain.push_back(atom);
        }
    }

    // assign the Hydrogen
    mH = GetN();

    if (mType != ResidueType::Proline and mPrev != nullptr)
    {
        const MAtom &pc = mPrev->GetC();
        const MAtom &po = mPrev->GetO();

        double       CODistance = Distance(pc, po);

        mH.mLoc += (pc.mLoc - po.mLoc).scale(1./CODistance);
    }

    // update the box containing all atoms
    mBox[0][XX] = mBox[0][YY] = mBox[0][ZZ] =  std::numeric_limits<double>::max();
    mBox[1][XX] = mBox[1][YY] = mBox[1][ZZ] = -std::numeric_limits<double>::max();

    ExtendBox(mN, kRadiusN + 2 * kRadiusWater);
    ExtendBox(mCA, kRadiusCA + 2 * kRadiusWater);
    ExtendBox(mC, kRadiusC + 2 * kRadiusWater);
    ExtendBox(mO, kRadiusO + 2 * kRadiusWater);
    for (const auto &atom : mSideChain)
    {
        ExtendBox(atom, kRadiusSideAtom + 2 * kRadiusWater);
    }

    mRadius = mBox[1][XX] - mBox[0][XX];
    if (mRadius < mBox[1][YY] - mBox[0][YY])
    {
        mRadius = mBox[1][YY] - mBox[0][YY];
    }
    if (mRadius < mBox[1][ZZ] - mBox[0][ZZ])
    {
        mRadius = mBox[1][ZZ] - mBox[0][ZZ];
    }

    mCenter[XX] = (mBox[0][XX] + mBox[1][XX]) / 2;
    mCenter[YY] = (mBox[0][YY] + mBox[1][YY]) / 2;
    mCenter[ZZ] = (mBox[0][ZZ] + mBox[1][ZZ]) / 2;

}

MResidue::MResidue(int32_t inNumber, char inTypeCode, MResidue* inPrevious)
    : mPrev(nullptr),
      mNext(nullptr),
      mSeqNumber(inNumber),
      mNumber(inNumber),
      mType(MapResidue(inTypeCode)),
      mSSBridgeNr(0),
      mAccessibility(0),
      mSecondaryStructure(SecondaryStructure::loop),
      mSheet(0),
      mBend(false),
      mRadius(0),
      mH(MAtom())
{
    std::fill(mHelixFlags, mHelixFlags + 3, helixNone);

    mBetaPartner[0].residue = mBetaPartner[1].residue = nullptr;

    mHBondDonor[0].energy  = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
    mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nullptr;

    static const MAtom kNullAtom = {};
    mN = mCA = mC = mO = kNullAtom;

    mCA.mResSeq  = inTypeCode;
    mCA.mChainID = "A";
}

MResidue::MResidue(const MResidue &residue)
    : mChainID(residue.mChainID),
      mPrev(nullptr),
      mNext(nullptr),
      mSeqNumber(residue.mSeqNumber),
      mNumber(residue.mNumber),
      mType(residue.mType),
      mSSBridgeNr(residue.mSSBridgeNr),
      mAccessibility(residue.mAccessibility),
      mSecondaryStructure(residue.mSecondaryStructure),
      mC(residue.mC),
      mN(residue.mN),
      mCA(residue.mCA),
      mO(residue.mO),
      mH(residue.mH),
      mSideChain(residue.mSideChain),
      mSheet(residue.mSheet),
      mBend(residue.mBend),
      mCenter(residue.mCenter),
      mRadius(residue.mRadius)
{
    std::copy(residue.mHBondDonor, residue.mHBondDonor + 2, mHBondDonor);
    std::copy(residue.mHBondAcceptor, residue.mHBondAcceptor + 2, mHBondAcceptor);
    std::copy(residue.mBetaPartner, residue.mBetaPartner + 2, mBetaPartner);
    std::copy(residue.mHelixFlags, residue.mHelixFlags + 3, mHelixFlags);
    std::copy(residue.mBox, residue.mBox + 2, mBox);
}

void MResidue::SetPrev(MResidue* inResidue)
{
    mPrev        = inResidue;
    mPrev->mNext = this;
}

bool MResidue::NoChainBreak(const MResidue* from, const MResidue* to)
{
    bool result = true;
    for (const MResidue* r = from; result and r != to; r = r->mNext)
    {
        MResidue* next = r->mNext;
        if (next == nullptr)
        {
            result = false;
        }
        else
        {
            result = next->mNumber == r->mNumber + 1;
        }
    }
    return result;
}

void MResidue::SetChainID(const std::string &inChainID)
{
    mChainID = inChainID;

    mC.SetChainID(inChainID);
    mCA.SetChainID(inChainID);
    mO.SetChainID(inChainID);
    mN.SetChainID(inChainID);
    mH.SetChainID(inChainID);
    std::for_each(mSideChain.begin(), mSideChain.end(),
                  [inChainID](MAtom &atom ){atom.SetChainID(inChainID); });
}

bool MResidue::ValidDistance(const MResidue &inNext) const
{
    return Distance(GetC(), inNext.GetN()) <= kMaxPeptideBondLength;
}

bool MResidue::TestBond(const MResidue* other) const
{
    return
        (mHBondAcceptor[0].residue == other and mHBondAcceptor[0].energy < kMaxHBondEnergy)or
            (mHBondAcceptor[1].residue == other and mHBondAcceptor[1].energy < kMaxHBondEnergy);
}

double MResidue::Phi() const
{
    double result = 360;
    if (mPrev != nullptr and NoChainBreak(mPrev, this))
    {
        result = DihedralAngle(mPrev->GetC(), GetN(), GetCAlpha(), GetC());
    }
    return result;
}

double MResidue::Psi() const
{
    double result = 360;
    if (mNext != nullptr and NoChainBreak(this, mNext))
    {
        result = DihedralAngle(GetN(), GetCAlpha(), GetC(), mNext->GetN());
    }
    return result;
}

std::tuple<double, char> MResidue::Alpha() const
{
    double          alhpa     = 360;
    char            chirality = ' ';

    const MResidue* nextNext = mNext ? mNext->Next() : nullptr;
    if (mPrev != nullptr and
        nextNext != nullptr and
        NoChainBreak(mPrev, nextNext))
    {
        alhpa = DihedralAngle(mPrev->GetCAlpha(), GetCAlpha(), mNext->GetCAlpha(),
                              nextNext->GetCAlpha());
        if (alhpa < 0)
        {
            chirality = '-';
        }
        else
        {
            chirality = '+';
        }
    }
    return std::make_tuple(alhpa, chirality);
}

double MResidue::Kappa() const
{
    double          result   = 360;
    const MResidue* prevPrev = mPrev ? mPrev->Prev() : nullptr;
    const MResidue* nextNext = mNext ? mNext->Next() : nullptr;
    if (prevPrev != nullptr and
        nextNext != nullptr and
        NoChainBreak(prevPrev, nextNext))
    {
        float  ckap = CosinusAngle(GetCAlpha(), prevPrev->GetCAlpha(),
                                   nextNext->GetCAlpha(), GetCAlpha());
        double skap = sqrt(1 - ckap * ckap);
        result = atan2(skap, ckap) * 180 / M_PI;
    }
    return result;
}

double MResidue::TCO() const
{
    double result = 0;
    if (mPrev != nullptr and NoChainBreak(mPrev, this))
    {
        result = CosinusAngle(GetC(), GetO(), mPrev->GetC(), mPrev->GetO());
    }
    return result;
}

void MResidue::SetBetaPartner(u_int32_t n,
                              MResidue* inResidue, u_int32_t inLadder, bool inParallel)
{
    // assert(n == 0 or n == 1);

    mBetaPartner[n].residue  = inResidue;
    mBetaPartner[n].ladder   = inLadder;
    mBetaPartner[n].parallel = inParallel;
}

MBridgeParner MResidue::GetBetaPartner(u_int32_t n) const
{
    // assert(n == 0 or n == 1);
    return mBetaPartner[n];
}

MHelixFlag MResidue::GetHelixFlag(u_int32_t inHelixStride) const
{
    // assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
    return mHelixFlags[inHelixStride - 3];
}

bool MResidue::IsHelixStart(u_int32_t inHelixStride) const
{
    // assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
    return mHelixFlags[inHelixStride - 3] == helixStart or
           mHelixFlags[inHelixStride - 3] == helixStartAndEnd;
}

void MResidue::SetHelixFlag(u_int32_t inHelixStride, MHelixFlag inHelixFlag)
{
    // assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
    mHelixFlags[inHelixStride - 3] = inHelixFlag;
}

void MResidue::SetSSBridgeNr(u_int8_t inBridgeNr)
{
    mSSBridgeNr = inBridgeNr;
}

u_int8_t MResidue::GetSSBridgeNr() const
{
    return mSSBridgeNr;
}

// TODO: use the angle to improve bond energy calculation.
double MResidue::CalculateHBondEnergy(MResidue &inDonor, MResidue &inAcceptor)
{
    float result = 0;

    if (inDonor.mType != ResidueType::Proline)
    {
        double distanceHO = Distance(inDonor.GetH(), inAcceptor.GetO());
        double distanceHC = Distance(inDonor.GetH(), inAcceptor.GetC());
        double distanceNC = Distance(inDonor.GetN(), inAcceptor.GetC());
        double distanceNO = Distance(inDonor.GetN(), inAcceptor.GetO());

        if (distanceHO < kMinimalDistance or distanceHC < kMinimalDistance or
            distanceNC < kMinimalDistance or distanceNO < kMinimalDistance)
        {
            result = kMinHBondEnergy;
        }
        else
        {
            result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;
        }

        // DSSP compatibility mode:
        result = std::round(result * 1000) / 1000;

        if (result < kMinHBondEnergy)
        {
            result = kMinHBondEnergy;
        }
    }

    // update donor
    if (result < inDonor.mHBondAcceptor[0].energy)
    {
        inDonor.mHBondAcceptor[1]         = inDonor.mHBondAcceptor[0];
        inDonor.mHBondAcceptor[0].residue = &inAcceptor;
        inDonor.mHBondAcceptor[0].energy  = result;
    }
    else if (result < inDonor.mHBondAcceptor[1].energy)
    {
        inDonor.mHBondAcceptor[1].residue = &inAcceptor;
        inDonor.mHBondAcceptor[1].energy  = result;
    }

    // and acceptor
    if (result < inAcceptor.mHBondDonor[0].energy)
    {
        inAcceptor.mHBondDonor[1]         = inAcceptor.mHBondDonor[0];
        inAcceptor.mHBondDonor[0].residue = &inDonor;
        inAcceptor.mHBondDonor[0].energy  = result;
    }
    else if (result < inAcceptor.mHBondDonor[1].energy)
    {
        inAcceptor.mHBondDonor[1].residue = &inDonor;
        inAcceptor.mHBondDonor[1].energy  = result;
    }

    return result;
}

MBridgeType MResidue::TestBridge(MResidue* test) const
{                                        // I.  a  d  II.  a  d    parallel
    const MResidue* a = mPrev;           //      \        /
    const MResidue* b = this;            //    b  e    b  e
    const MResidue* c = mNext;           //       /        \                      ..
    const MResidue* d = test->mPrev;     //    c  f    c  f
    const MResidue* e = test;            //
    const MResidue* f = test->mNext;     // III.  a <- f  IV. a    f    antiparallel
    //
    MBridgeType     result = btNoBridge; //    b   e      b <-> e
    //
    //    c -> d    c     d

    if (a and c and NoChainBreak(a, c) and d and f and NoChainBreak(d, f))
    {
        if ((TestBond(c, e) and TestBond(e, a))or
                (TestBond(f, b) and TestBond(b, d)))
        {
            result = btParallel;
        }
        else if ((TestBond(c, d) and TestBond(f, a))or
                     (TestBond(e, b) and TestBond(b, e)))
        {
            result = btAntiParallel;
        }
    }

    return result;
}

void MResidue::ExtendBox(const MAtom &atom, double inRadius)
{
    if (mBox[0][XX] > atom.mLoc[XX] - inRadius)
    {
        mBox[0][XX] = atom.mLoc[XX] - inRadius;
    }
    if (mBox[0][YY] > atom.mLoc[YY] - inRadius)
    {
        mBox[0][YY] = atom.mLoc[YY] - inRadius;
    }
    if (mBox[0][ZZ] > atom.mLoc[ZZ] - inRadius)
    {
        mBox[0][ZZ] = atom.mLoc[ZZ] - inRadius;
    }
    if (mBox[1][XX] < atom.mLoc[XX] + inRadius)
    {
        mBox[1][XX] = atom.mLoc[XX] + inRadius;
    }
    if (mBox[1][YY] < atom.mLoc[YY] + inRadius)
    {
        mBox[1][YY] = atom.mLoc[YY] + inRadius;
    }
    if (mBox[1][ZZ] < atom.mLoc[ZZ] + inRadius)
    {
        mBox[1][ZZ] = atom.mLoc[ZZ] + inRadius;
    }
}

inline
bool MResidue::AtomIntersectsBox(const MAtom &atom, double inRadius) const
{
    return
        atom.mLoc[XX] + inRadius >= mBox[0][XX] and
        atom.mLoc[XX] - inRadius <= mBox[1][XX] and
        atom.mLoc[YY] + inRadius >= mBox[0][YY] and
        atom.mLoc[YY] - inRadius <= mBox[1][YY] and
        atom.mLoc[ZZ] + inRadius >= mBox[0][ZZ] and
        atom.mLoc[ZZ] - inRadius <= mBox[1][ZZ];
}

void MResidue::CalculateSurface(const std::vector<MResidue*> &inResidues)
{
    std::vector<MResidue*> neighbours;

    for (const auto &r : inResidues)
    {
        gmx::RVec center;
        double    radius;
        r->GetCenterAndRadius(center, radius);

        if (Distance(mCenter, center) < mRadius + radius)
        {
            neighbours.push_back(r);
        }
    }

    mAccessibility = CalculateSurface(mN, kRadiusN, neighbours) +
        CalculateSurface(mCA, kRadiusCA, neighbours) +
        CalculateSurface(mC, kRadiusC, neighbours) +
        CalculateSurface(mO, kRadiusO, neighbours);

    for (const auto &atom : mSideChain)
    {
        mAccessibility += CalculateSurface(atom, kRadiusSideAtom, neighbours);
    }
}

class MAccumulator
{
    public:

        struct candidate
        {
            gmx::RVec  location;
            double     radius;
            double     distance;

            bool operator<(const candidate &rhs) const
            { return distance < rhs.distance; }
        };

        void operator()(const gmx::RVec &a, const gmx::RVec &b, double d, double r)
        {
            double distance = DistanceSquared(a, b);

            d += kRadiusWater;
            r += kRadiusWater;

            double test = d + r;
            test *= test;

            if (distance < test and distance > 0.0001)
            {
                candidate c = { b - a, r * r, distance };

                m_x.push_back(c);
                push_heap(m_x.begin(), m_x.end());
            }
        }

        void sort()
        {
            sort_heap(m_x.begin(), m_x.end());
        }

        std::vector<candidate>  m_x;
};

double MResidue::CalculateSurface(const MAtom &inAtom, double inRadius,
                                  const std::vector<MResidue*> &inResidues)
{
    MAccumulator accumulate;

    for (const auto &r :  inResidues)
    {
        if (r->AtomIntersectsBox(inAtom, inRadius))
        {
            accumulate(inAtom, r->mN, inRadius, kRadiusN);
            accumulate(inAtom, r->mCA, inRadius, kRadiusCA);
            accumulate(inAtom, r->mC, inRadius, kRadiusC);
            accumulate(inAtom, r->mO, inRadius, kRadiusO);

            for (const auto &atom : r->mSideChain)
            {
                accumulate(inAtom, atom, inRadius, kRadiusSideAtom);
            }
        }
    }

    accumulate.sort();

    double        radius  = inRadius + kRadiusWater;
    double        surface = 0;

    MSurfaceDots &surfaceDots = MSurfaceDots::Instance();

    for (u_int32_t i = 0; i < surfaceDots.size(); ++i)
    {
        gmx::RVec xx = surfaceDots[i].scale(radius);

        bool      free = true;
        for (u_int32_t k = 0; free and k < accumulate.m_x.size(); ++k)
        {
            free = accumulate.m_x[k].radius < DistanceSquared(xx, accumulate.m_x[k].location);
        }

        if (free)
        {
            surface += surfaceDots.weight();
        }
    }

    return surface * radius * radius;
}

// --------------------------------------------------------------------

MChain::MChain(const MChain &chain)
    : mChainID(chain.mChainID)
{
    MResidue* previous = nullptr;

    for (const auto &residue : chain.mResidues)
    {
        MResidue* newResidue = new MResidue(*residue);
        newResidue->SetPrev(previous);
        mResidues.push_back(newResidue);
        previous = newResidue;
    }
}

MChain::~MChain()
{
}

MChain &MChain::operator=(const MChain &chain)
{
    mResidues.clear();

    for (const auto &residue : chain.mResidues)
    {
        mResidues.push_back(new MResidue(*residue));
    }

    mChainID = chain.mChainID;

    return *this;
}

void MChain::SetChainID(const std::string &inChainID)
{
    mChainID = inChainID;
    std::for_each(mResidues.begin(), mResidues.end(), [inChainID](MResidue * residue){residue->SetChainID(inChainID); });
}

void MChain::SetAuthChainID(const std::string &inAuthChainID)
{
    mAuthChainID = inAuthChainID;
}

std::string MChain::GetAuthChainID(void) const
{
    return mAuthChainID;
}

const MResidue* MChain::GetResidueBySeqNumber(u_int16_t           inSeqNumber,
                                              const std::string  &inInsertionCode) const
{
    const auto r = std::find_if(mResidues.begin(), mResidues.end(),
                                [inSeqNumber, inInsertionCode](const MResidue * residue){
                                    return (residue->GetSeqNumber() == inSeqNumber) &&
                                    (residue->GetInsertionCode() == inInsertionCode);
                                }
                                );
    return *r;
}

// --------------------------------------------------------------------

struct MResidueID
{
    std::string  chain;
    u_int16_t    seqNumber;
    std::string  insertionCode;

    bool operator<(const MResidueID &o) const
    {
        return
            chain < o.chain or
                (chain == o.chain and seqNumber < o.seqNumber) or
                (chain == o.chain and seqNumber == o.seqNumber and
                insertionCode < o.insertionCode);
    }

    bool operator!=(const MResidueID &o) const
    {
        return chain != o.chain or
               seqNumber != o.seqNumber or
               insertionCode != o.insertionCode;
    }
};

MProtein::MProtein()
    : mResidueCount(0),
      mID("unknown"),
      mChainBreaks(0)
{
}

MProtein::MProtein(const std::string &inID, MChain* inChain)
    : mID(inID),
      mChainBreaks(0),
      mResidueCount(0)
{
    mChains.push_back(inChain);
}

MProtein::~MProtein()
{
}

void MProtein::ReadPDB(const std::string &pdbString, bool cAlphaOnly)
{
    std::istringstream is(pdbString);
    mResidueCount                    = 0;
    mChainBreaks                     = 0;

    std::vector < std::pair < MResidueID, MResidueID>> ssbonds;
    std::set<char>         terminatedChains;

    bool                   model = false;
    std::vector<MAtom>     atoms;
    char                   firstAltLoc = 0;
    std::unique_ptr<MAtom> prevAtom;

    while (not is.eof())
    {
        std::string line;
        getline(is, line);
        if (line.empty() and is.eof())
        {
            break;
        }

        if (line.find( "HEADER") == 0)
        {
            mHeader = line;
            trim(mHeader);
            if (line.length() >= 66)
            {
                mID = line.substr(62, 4);
            }
            else
            {
                mID = "UNDF";
            }
            continue;
        }

        if (line.find( "COMPND") == 0)
        {
            trim_right(line);
            if (line.length() >= 10)
            {
                mCompound = mCompound + line.substr(10);
            }
            continue;
        }

        if (line.find( "SOURCE") == 0)
        {
            trim_right(line);
            if (line.length() >= 10)
            {
                mSource = mSource + line.substr(10);
            }
            continue;
        }

        if (line.find( "AUTHOR") == 0)
        {
            trim_right(line);
            if (line.length() >= 10)
            {
                mAuthor = mAuthor + line.substr(10);
            }
            continue;
        }

        if (line.find( "DBREF") == 0)
        {
            trim(line);
            mDbRef.push_back(line);
            continue;
        }

        // brain dead support for only the first model in the file (NMR)
        if (line.find( "MODEL") == 0)
        {
            model = true;
            continue;
        }

        if (line.find( "ENDMDL") == 0 and model == true)
        {
            break;
        }

        if (line.find( "SSBOND") == 0)
        {
            //SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.03
            std::pair<MResidueID, MResidueID> ssbond;
            ssbond.first.chain     = line[15];
            ssbond.first.seqNumber = std::stoi(
                        trim_copy(line.substr(16, 5)));
            ssbond.first.insertionCode = line[21];
            ssbond.second.chain        = line[29];
            ssbond.second.seqNumber    = std::stoi(
                        trim_copy(line.substr(30, 5)));
            ssbond.second.insertionCode = line[35];

            ssbonds.push_back(ssbond);
            continue;
        }

        if (line.find( "TER   ") == 0)
        {
            if (atoms.empty())
            {
                std::cerr << "no atoms read before TER record " << std::endl
                << line << std::endl;
                continue;
            }

            AddResidue(atoms);
            atoms.clear();
            firstAltLoc = 0;
            prevAtom    = nullptr;

            terminatedChains.insert(line[21]);

            continue;
        }

        if (line.find( "ATOM  ") == 0 or line.find( "HETATM") == 0)
        //  1 - 6  Record name "ATOM "
        {
            if (cAlphaOnly and line.substr(12, 4) != " CA ")
            {
                continue;
            }

            MAtom atom = {};

            //  7 - 11  Integer serial Atom serial number.
            atom.mSerial = std::stoi(
                        trim_copy(line.substr(6, 5)));
            //  13 - 16  Atom name Atom name.
            atom.mName = trim_copy(line.substr(12, 4));
            //  17    Character altLoc Alternate location indicator.
            atom.mAltLoc = line[16];
            //  18 - 20  Residue name resName Residue name.
            atom.mResName = trim_copy(line.substr(17, 4));
            //  22    Character chainID Chain identifier.
            atom.mChainID     = line[21];
            atom.mAuthChainID = atom.mChainID;
            //  23 - 26  Integer resSeq Residue sequence number.
            atom.mResSeq = std::stoi(
                        trim_copy(line.substr(22, 4)));
            //  27    AChar iCode Code for insertion of residues.
            atom.mICode = line.substr(26, 1);

            //  31 - 38  Real(8.3) x Orthogonal coordinates for X in Angstroms.
            atom.mLoc[XX] = ParseFloat(line.substr(30, 8));
            //  39 - 46  Real(8.3) y Orthogonal coordinates for Y in Angstroms.
            atom.mLoc[YY] = ParseFloat(line.substr(38, 8));
            //  47 - 54  Real(8.3) z Orthogonal coordinates for Z in Angstroms.
            atom.mLoc[ZZ] = ParseFloat(line.substr(46, 8));

            //  55 - 60  Real(6.2) occupancy Occupancy.
            if (line.length() > 54)
            {
                atom.mOccupancy = ParseFloat(line.substr(54, 6));
            }

            //  61 - 66  Real(6.2) tempFactor Temperature factor.
            if (line.length() > 60)
            {
                atom.mTempFactor = ParseFloat(line.substr(60, 6));
            }

            //  77 - 78  LString(2) element Element symbol, right-justified.
            if (line.length() > 76)
            {
                atom.mElement = trim_copy(line.substr(76, 3));
            }
            //  79 - 80  LString(2) charge Charge on the atom.
            atom.mCharge = 0;

//      alternative test, check chain ID as well.
            if (prevAtom
                &&
                (
                    atom.mChainID != prevAtom->mChainID
                    ||
                    atom.mResSeq  != prevAtom->mResSeq
                    ||
                    atom.mICode   != prevAtom->mICode
                )
                )
//      if (not atoms.empty() and
//        (atom.mResSeq != atoms.back().mResSeq or (atom.mResSeq == atoms.back().mResSeq and atom.mICode != atoms.back().mICode)))
            {
                if (!atoms.empty() )
                {
                    AddResidue(atoms);
                    atoms.clear();
                }
                firstAltLoc = 0;
                prevAtom    = nullptr;
            }

            try
            {
                atom.mType = MapElement(line.substr(76, 2));
            }
            catch (const std::exception &e)
            {
                atom.mType = kUnknownAtom;
            }

            if (atom.mType == kHydrogen)
            {
                continue;
            }
            prevAtom = std::unique_ptr<MAtom>(new MAtom(atom));

            if (atom.mAltLoc != ' ')
            {
                if (firstAltLoc == 0)
                {
                    firstAltLoc = atom.mAltLoc;
                }
                if (atom.mAltLoc == firstAltLoc)
                {
                    atom.mAltLoc = 'A';
                }
            }

            if (firstAltLoc != 0 and
                atom.mAltLoc != ' ' and
                atom.mAltLoc != firstAltLoc)
            {
                continue;
            }

            atoms.push_back(atom);
        }
    }

    if (not atoms.empty()) // we have read atoms without a TER
    {
        AddResidue(atoms);
    }
    mChains.erase(
            std::remove_if(mChains.begin(), mChains.end(), [](MChain *chain){return chain->Empty(); }),
            mChains.end());

}

std::string MProtein::GetCompound() const
{
    std::string result("COMPND    ");
    result += mCompound;
    return result.substr(0, 80);
}

std::string MProtein::GetSource() const
{
    std::string result("SOURCE    ");
    result += mSource;
    return result.substr(0, 80);
}

std::string MProtein::GetAuthor() const
{
    std::string result("AUTHOR    ");
    result += mAuthor;
    return result.substr(0, 80);
}

void MProtein::GetStatistics(u_int32_t &outNrOfResidues, u_int32_t &outNrOfChains,
                             u_int32_t &outNrOfIntraChainSSBridges,
                             u_int32_t &outNrOfHBonds,
                             u_int32_t outNrOfHBondsPerDistance[11]) const
{
    outNrOfResidues  = mResidueCount;
    outNrOfChains    = mChains.size() + mChainBreaks;

    outNrOfIntraChainSSBridges = 0;
    for (std::vector<std::pair<MResidue*, MResidue*> >::const_iterator ri = mSSBonds.begin(); ri != mSSBonds.end(); ++ri)
    {
        if (ri->first->GetChainID() == ri->second->GetChainID() and
                (MResidue::NoChainBreak(ri->first, ri->second) or
                MResidue::NoChainBreak(ri->first, ri->second)))
        {
            ++outNrOfIntraChainSSBridges;
        }
    }

    outNrOfHBonds = 0;
    for (const MChain* chain : mChains)
    {
        for (const MResidue* r : chain->GetResidues())
        {
            const HBond* donor = r->Donor();

            for (u_int32_t i = 0; i < 2; ++i)
            {
                if (donor[i].residue != nullptr and donor[i].energy < kMaxHBondEnergy)
                {
                    ++outNrOfHBonds;
                    int32_t k = donor[i].residue->GetNumber() - r->GetNumber();
                    if (k >= -5 and k <= 5)
                    {
                        outNrOfHBondsPerDistance[k + 5] += 1;
                    }
                }
            }
        }
    }
}

void MProtein::AddResidue(const std::vector<MAtom> &inAtoms)
{
    bool hasN = false, hasCA = false, hasC = false, hasO = false;
    for (const MAtom &atom : inAtoms)
    {
        if (not hasN and atom.GetName() == "N")
        {
            hasN = true;
        }
        if (not hasCA and atom.GetName() == "CA")
        {
            hasCA = true;
        }
        if (not hasC and atom.GetName() == "C")
        {
            hasC = true;
        }
        if (not hasO and atom.GetName() == "O")
        {
            hasO = true;
        }
    }

    if (hasN and hasCA and hasC and hasO)
    {
        MChain &chain = GetChain(inAtoms.front().mChainID);
        chain.SetAuthChainID(inAtoms.front().mAuthChainID);

        std::vector<MResidue*> &residues(chain.GetResidues());

        MResidue              * prev = nullptr;
        if (not residues.empty())
        {
            prev = residues.back();
        }

        int32_t   resNumber = mResidueCount + mChains.size() + mChainBreaks;
        MResidue* r         = new MResidue(resNumber, prev, inAtoms);
        // check for chain breaks
        if (prev != nullptr and not prev->ValidDistance(*r))
        {
            ++mChainBreaks;
            r->SetNumber(resNumber + 1);
        }

        residues.push_back(r);
        ++mResidueCount;
    }

}

const MChain &MProtein::GetChain(const std::string &inChainID) const
{
    for (u_int32_t i = 0; i < mChains.size(); ++i)
    {
        if (mChains[i]->GetChainID() == inChainID)
        {
            return *mChains[i];
        }
    }

    return *mChains.front();
}

MChain &MProtein::GetChain(const std::string &inChainID)
{
    for (u_int32_t i = 0; i < mChains.size(); ++i)
    {
        if (mChains[i]->GetChainID() == inChainID)
        {
            return *mChains[i];
        }
    }

    mChains.push_back(new MChain(inChainID));
    return *mChains.back();
}

void MProtein::CalculateSecondaryStructure(bool inPreferPiHelices)
{
    std::vector<MResidue*> residues;
    residues.reserve(mResidueCount);
    for (const MChain* chain : mChains)
    {
        residues.insert(residues.end(), chain->GetResidues().begin(),
                        chain->GetResidues().end());
    }

    CalculateAccessibilities(residues);
    CalculateHBondEnergies(residues);
    CalculateBetaSheets(residues);
    CalculateAlphaHelices(residues, inPreferPiHelices);
}

void MProtein::CalculateHBondEnergies(const std::vector<MResidue*> &inResidues)
{
    // Calculate the HBond energies
    for (u_int32_t i = 0; i + 1 < inResidues.size(); ++i)
    {
        MResidue* ri = inResidues[i];

        for (u_int32_t j = i + 1; j < inResidues.size(); ++j)
        {
            MResidue* rj = inResidues[j];

            if (Distance(ri->GetCAlpha(), rj->GetCAlpha()) < kMinimalCADistance)
            {
                MResidue::CalculateHBondEnergy(*ri, *rj);
                if (j != i + 1)
                {
                    MResidue::CalculateHBondEnergy(*rj, *ri);
                }
            }
        }
    }
}

// TODO: improve alpha helix calculation by better recognizing pi-helices
void MProtein::CalculateAlphaHelices(const std::vector<MResidue*> &inResidues,
                                     bool                          inPreferPiHelices)
{
    // Helix and Turn
    for (const MChain* chain : mChains)
    {
        for (u_int32_t stride = 3; stride <= 5; ++stride)
        {
            std::vector<MResidue*> res(chain->GetResidues());
            if (res.size() < stride)
            {
                continue;
            }

            for (u_int32_t i = 0; i + stride < res.size(); ++i)
            {
                if (MResidue::TestBond(res[i + stride], res[i]) and
                    MResidue::NoChainBreak(res[i], res[i + stride]))
                {
                    res[i + stride]->SetHelixFlag(stride, helixEnd);
                    for (u_int32_t j = i + 1; j < i + stride; ++j)
                    {
                        if (res[j]->GetHelixFlag(stride) == helixNone)
                        {
                            res[j]->SetHelixFlag(stride, helixMiddle);
                        }
                    }

                    if (res[i]->GetHelixFlag(stride) == helixEnd)
                    {
                        res[i]->SetHelixFlag(stride, helixStartAndEnd);
                    }
                    else
                    {
                        res[i]->SetHelixFlag(stride, helixStart);
                    }
                }
            }
        }
    }

    for (MResidue* r : inResidues)
    {
        double kappa = r->Kappa();
        r->SetBend(kappa != 360 and kappa > 70);
    }

    for (u_int32_t i = 1; i + 4 < inResidues.size(); ++i)
    {
        if (inResidues[i]->IsHelixStart(4) and inResidues[i - 1]->IsHelixStart(4))
        {
            for (u_int32_t j = i; j <= i + 3; ++j)
            {
                inResidues[j]->SetSecondaryStructure(SecondaryStructure::alphahelix);
            }
        }
    }

    for (u_int32_t i = 1; i + 3 < inResidues.size(); ++i)
    {
        if (inResidues[i]->IsHelixStart(3) and inResidues[i - 1]->IsHelixStart(3))
        {
            bool empty = true;
            for (u_int32_t j = i; empty and j <= i + 2; ++j)
            {
                empty = inResidues[j]->GetSecondaryStructure() == SecondaryStructure::loop or
                    inResidues[j]->GetSecondaryStructure() == SecondaryStructure::helix_3;
            }
            if (empty)
            {
                for (u_int32_t j = i; j <= i + 2; ++j)
                {
                    inResidues[j]->SetSecondaryStructure(SecondaryStructure::helix_3);
                }
            }
        }
    }

    for (u_int32_t i = 1; i + 5 < inResidues.size(); ++i)
    {
        if (inResidues[i]->IsHelixStart(5) and inResidues[i - 1]->IsHelixStart(5))
        {
            bool empty = true;
            for (u_int32_t j = i; empty and j <= i + 4; ++j)
            {
                empty = inResidues[j]->GetSecondaryStructure() == SecondaryStructure::loop or
                    inResidues[j]->GetSecondaryStructure() == SecondaryStructure::helix_5 or
                        (inPreferPiHelices and
                        inResidues[j]->GetSecondaryStructure() == SecondaryStructure::alphahelix);
            }
            if (empty)
            {
                for (u_int32_t j = i; j <= i + 4; ++j)
                {
                    inResidues[j]->SetSecondaryStructure(SecondaryStructure::helix_5);
                }
            }
        }
    }

    for (u_int32_t i = 1; i + 1 < inResidues.size(); ++i)
    {
        if (inResidues[i]->GetSecondaryStructure() == SecondaryStructure::loop)
        {
            bool isTurn = false;
            for (u_int32_t stride = 3; stride <= 5 and not isTurn; ++stride)
            {
                for (u_int32_t k = 1; k < stride and not isTurn; ++k)
                {
                    isTurn = (i >= k) and inResidues[i - k]->IsHelixStart(stride);
                }
            }

            if (isTurn)
            {
                inResidues[i]->SetSecondaryStructure(SecondaryStructure::turn);
            }
            else if (inResidues[i]->IsBend())
            {
                inResidues[i]->SetSecondaryStructure(SecondaryStructure::bend);
            }
        }
    }
}

void MProtein::CalculateBetaSheets(const std::vector<MResidue*> &inResidues)
{
    // Calculate Bridges
    std::vector<MBridge> bridges;
    if (inResidues.size() > 4)
    {
        for (u_int32_t i = 1; i + 4 < inResidues.size(); ++i)
        {
            MResidue* ri = inResidues[i];

            for (u_int32_t j = i + 3; j + 1 < inResidues.size(); ++j)
            {
                MResidue  * rj = inResidues[j];

                MBridgeType type = ri->TestBridge(rj);
                if (type == btNoBridge)
                {
                    continue;
                }

                bool found = false;
                for (MBridge &bridge : bridges)
                {
                    if (type != bridge.type or i != bridge.i.back() + 1)
                    {
                        continue;
                    }

                    if (type == btParallel and bridge.j.back() + 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_back(j);
                        found = true;
                        break;
                    }

                    if (type == btAntiParallel and bridge.j.front() - 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_front(j);
                        found = true;
                        break;
                    }
                }

                if (not found)
                {
                    MBridge bridge = {};

                    bridge.type = type;
                    bridge.i.push_back(i);
                    bridge.chainI = ri->GetChainID();
                    bridge.j.push_back(j);
                    bridge.chainJ = rj->GetChainID();

                    bridges.push_back(bridge);
                }
            }
        }
    }

    // extend ladders
    sort(bridges.begin(), bridges.end());

    for (u_int32_t i = 0; i < bridges.size(); ++i)
    {
        for (u_int32_t j = i + 1; j < bridges.size(); ++j)
        {
            u_int32_t ibi = bridges[i].i.front();
            u_int32_t iei = bridges[i].i.back();
            u_int32_t jbi = bridges[i].j.front();
            u_int32_t jei = bridges[i].j.back();
            u_int32_t ibj = bridges[j].i.front();
            u_int32_t iej = bridges[j].i.back();
            u_int32_t jbj = bridges[j].j.front();
            u_int32_t jej = bridges[j].j.back();

            if (bridges[i].type != bridges[j].type or
                MResidue::NoChainBreak(inResidues[std::min(ibi, ibj)],
                                       inResidues[std::max(iei, iej)]) == false or
                MResidue::NoChainBreak(inResidues[std::min(jbi, jbj)],
                                       inResidues[std::max(jei, jej)]) == false or
                ibj - iei >= 6 or
                    (iei >= ibj and ibi <= iej))
            {
                continue;
            }

            bool bulge;
            if (bridges[i].type == btParallel)
            {
                bulge = ((jbj - jei < 6 and ibj - iei < 3)or (jbj - jei < 3));
            }
            else
            {
                bulge = ((jbi - jej < 6 and ibj - iei < 3)or (jbi - jej < 3));
            }

            if (bulge)
            {
                bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(),
                                    bridges[j].i.end());
                if (bridges[i].type == btParallel)
                {
                    bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(),
                                        bridges[j].j.end());
                }
                else
                {
                    bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(),
                                        bridges[j].j.end());
                }
                bridges.erase(bridges.begin() + j);
                --j;
            }
        }
    }

    // Sheet
    std::set<MBridge*> ladderset;
    for (MBridge &bridge : bridges)
    {
        ladderset.insert(&bridge);
    }

    u_int32_t sheet = 1, ladder = 0;
    while (not ladderset.empty())
    {
        std::set<MBridge*> sheetset;
        sheetset.insert(*ladderset.begin());
        ladderset.erase(ladderset.begin());

        bool done = false;
        while (not done)
        {
            done = true;
            for (MBridge* a : sheetset)
            {
                for (MBridge* b : ladderset)
                {
                    if (Linked(*a, *b))
                    {
                        sheetset.insert(b);
                        ladderset.erase(b);
                        done = false;
                        break;
                    }
                }
                if (not done)
                {
                    break;
                }
            }
        }

        for (MBridge* bridge : sheetset)
        {
            bridge->ladder = ladder;
            bridge->sheet  = sheet;
            bridge->link   = sheetset;

            ++ladder;
        }

        ++sheet;
    }

    for (MBridge &bridge : bridges)
    {
        // find out if any of the i and j set members already have
        // a bridge assigned, if so, we're assigning bridge 2

        u_int32_t betai = 0, betaj = 0;

        for (u_int32_t l : bridge.i)
        {
            if (inResidues[l]->GetBetaPartner(0).residue != nullptr)
            {
                betai = 1;
                break;
            }
        }

        for (u_int32_t l : bridge.j)
        {
            if (inResidues[l]->GetBetaPartner(0).residue != nullptr)
            {
                betaj = 1;
                break;
            }
        }

        SecondaryStructure ss = SecondaryStructure::betabridge;
        if (bridge.i.size() > 1)
        {
            ss = SecondaryStructure::strand;
        }

        if (bridge.type == btParallel)
        {
            std::deque<u_int32_t>::iterator j = bridge.j.begin();
            for (u_int32_t i : bridge.i)
            {
                inResidues[i]->SetBetaPartner(betai, inResidues[*j++], bridge.ladder,
                                              true);
            }

            j = bridge.i.begin();
            for (u_int32_t i : bridge.j)
            {
                inResidues[i]->SetBetaPartner(betaj, inResidues[*j++], bridge.ladder,
                                              true);
            }
        }
        else
        {
            std::deque<u_int32_t>::reverse_iterator j = bridge.j.rbegin();
            for (u_int32_t i : bridge.i)
            {
                inResidues[i]->SetBetaPartner(betai, inResidues[*j++], bridge.ladder,
                                              false);
            }

            j = bridge.i.rbegin();
            for (u_int32_t i : bridge.j)
            {
                inResidues[i]->SetBetaPartner(betaj, inResidues[*j++], bridge.ladder,
                                              false);
            }
        }

        for (u_int32_t i = bridge.i.front(); i <= bridge.i.back(); ++i)
        {
            if (inResidues[i]->GetSecondaryStructure() != SecondaryStructure::strand)
            {
                inResidues[i]->SetSecondaryStructure(ss);
            }
            inResidues[i]->SetSheet(bridge.sheet);
        }

        for (u_int32_t i = bridge.j.front(); i <= bridge.j.back(); ++i)
        {
            if (inResidues[i]->GetSecondaryStructure() != SecondaryStructure::strand)
            {
                inResidues[i]->SetSecondaryStructure(ss);
            }
            inResidues[i]->SetSheet(bridge.sheet);
        }
    }
}

void MProtein::CalculateAccessibilities(
        const std::vector<MResidue*> &inResidues)
{
    for (MResidue* residue : inResidues)
    {
        residue->CalculateSurface(inResidues);
    }
}

void MProtein::SetChain(const std::string &inChainID, const MChain &inChain)
{
    MChain &chain(GetChain(inChainID));
    chain = inChain;
    chain.SetChainID(inChainID);
}

// Non-const overload, implemented in terms of the const overload
MResidue* MProtein::GetResidue(const std::string  &inChainID,
                               u_int16_t           inSeqNumber,
                               const std::string  &inInsertionCode)
{
    return const_cast<MResidue *>( static_cast<const MProtein &>( *this ).GetResidue(
                                           inChainID,
                                           inSeqNumber,
                                           inInsertionCode
                                           ) );
}

// Const overload
const MResidue* MProtein::GetResidue(const std::string  &inChainID,
                                     u_int16_t           inSeqNumber,
                                     const std::string  &inInsertionCode) const
{
    const MChain &chain = GetChain(inChainID);
    return chain.GetResidueBySeqNumber(inSeqNumber, inInsertionCode);
}

void MProtein::GetCAlphaLocations(const std::string      &inChainID,
                                  std::vector<gmx::RVec> &outPoints) const
{
    std::string chainID = inChainID;
    if (chainID.empty())
    {
        chainID = mChains.front()->GetChainID();
    }

    for (const MResidue* r : GetChain(chainID).GetResidues())
    {
        outPoints.push_back(gmx::RVec(r->GetCAlpha()));
    }
}

gmx::RVec MProtein::GetCAlphaPosition(const std::string &inChainID, int16_t inPDBResSeq) const
{
    std::string chainID = inChainID;
    if (chainID.empty())
    {
        chainID = mChains.front()->GetChainID();
    }

    gmx::RVec result;
    for (const MResidue* r : GetChain(chainID).GetResidues())
    {
        if (r->GetSeqNumber() != inPDBResSeq)
        {
            continue;
        }

        result = r->GetCAlpha();
    }

    return result;
}


std::string ResidueToDSSPLine(const MResidue &residue)
{
/*
   This is the header line for the residue lines in a DSSP file:

 #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA           CHAIN AUTHCHAIN
 */

    const MAtom &ca = residue.GetCAlpha();

    char         ss = static_cast<char>(residue.GetSecondaryStructure());

    char         helix[3];
    for (u_int32_t stride = 3; stride <= 5; ++stride)
    {
        switch (residue.GetHelixFlag(stride))
        {
            case helixNone:      helix[stride - 3]    = ' '; break;
            case helixStart:    helix[stride - 3]     = '>'; break;
            case helixEnd:      helix[stride - 3]     = '<'; break;
            case helixStartAndEnd:  helix[stride - 3] = 'X'; break;
            case helixMiddle:    helix[stride - 3]    = '0' + stride; break;
        }
    }

    char bend = ' ';
    if (residue.IsBend())
    {
        bend = 'S';
    }

    double alpha;
    char   chirality;
    std::tie(alpha, chirality) = residue.Alpha();

    u_int32_t bp[2] = {};
    char      bridgelabel[2] = { ' ', ' ' };
    for (u_int32_t i = 0; i < 2; ++i)
    {
        MBridgeParner p = residue.GetBetaPartner(i);
        if (p.residue != nullptr)
        {
            bp[i]          = p.residue->GetNumber();
            bp[i]         %= 10000; // won't fit otherwise...
            bridgelabel[i] = 'A' + p.ladder % 26;
            if (p.parallel)
            {
                bridgelabel[i] = tolower(bridgelabel[i]);
            }
        }
    }

    char sheet = ' ';
    if (residue.GetSheet() != 0)
    {
        sheet = 'A' + (residue.GetSheet() - 1) % 26;
    }

    std::string  NHO[2], ONH[2];
    const HBond* acceptors = residue.Acceptor();
    const HBond* donors    = residue.Donor();
    for (u_int32_t i = 0; i < 2; ++i)
    {
        NHO[i] = ONH[i] = "0, 0.0";

        if (acceptors[i].residue != nullptr)
        {
            int32_t d = acceptors[i].residue->GetNumber() - residue.GetNumber();
            char    buf[1000];
            sprintf(buf, "%d,%3.1f", d, acceptors[i].energy);
            NHO[i] = buf;
        }

        if (donors[i].residue != nullptr)
        {
            int32_t d = donors[i].residue->GetNumber() - residue.GetNumber();
            char    buf[1000];
            sprintf(buf, "%d,%3.1f", d, donors[i].energy);
            ONH[i] = buf;
        }
    }

    std::string chainChar     = ca.mChainID,
                long_ChainID1 = ca.mChainID,
                long_ChainID2 = ca.mAuthChainID;
    if (ca.mChainID.length () > 1)
    {
        // For mmCIF compatibility

        chainChar = ">";
    }

    char buf[1000];
    sprintf(buf, "%5.5d%5.5d%1.1s%1.1s %c  %c ", residue.GetNumber(), ca.mResSeq, ca.mICode.c_str(), chainChar.c_str(), ' ', ss);

    char buf2[1000];
    sprintf(buf2, "%s%c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d ", buf, helix[0], helix[1], helix[2], bend, chirality, bridgelabel[0], bridgelabel[1],
            bp[0], bp[1], sheet, static_cast<int>(residue.Accessibility() + 0.5) );

    char buf3[1000];
    sprintf(buf3, "%s%11s%11s%11s%11s  %6.3f%6.1f%6.1f", buf2,   NHO[0].c_str(), ONH[0].c_str(), NHO[1].c_str(), ONH[1].c_str(),
            residue.TCO(), residue.Kappa(), alpha );

    char buf4[1000];
    sprintf(buf4, "%s%6.1f%6.1f %6.1f %6.1f %6.1f             %4.4s      %4.4s", buf3, residue.Phi(), residue.Psi(),
            ca.mLoc[XX], ca.mLoc[YY], ca.mLoc[ZZ], long_ChainID1.c_str(), long_ChainID2.c_str());

    return buf4;
}

std::vector<std::string> MProtein::writeSecondaryStructure()
{

    std::vector<const MResidue*> residues;
    for (const auto &chain : GetChains())
    {
        for (const auto &residue : chain->GetResidues())
        {
            residues.push_back(residue);
        }
    }
    // keep residues sorted by residue number as assigned during reading the PDB file
    sort(residues.begin(), residues.end(), [](const MResidue * a, const MResidue * b){return a->GetNumber() < b->GetNumber(); });

    const MResidue         * last                     = nullptr;
    std::string              secondaryStructureString = {};
    std::vector<std::string> result;
    for (const auto &residue : residues)
    {
        // insert a break line whenever we detect missing residues
        // can be the transition to a different chain, or missing residues in the current chain
        if (last != nullptr and last->GetNumber() + 1 != residue->GetNumber())
        {
            if (last->GetChainID() != residue->GetChainID())
            {
                result.push_back(secondaryStructureString);
                secondaryStructureString = {};
            }
            else
            {
                secondaryStructureString += ' ';
            }
        }
        secondaryStructureString += static_cast<char>(residue->GetSecondaryStructure());
        last                      = residue;
    }
    return result;
}

std::string WriteDSSP(MProtein &protein)
{
    std::string  dsspString;
    u_int32_t    nrOfResidues, nrOfChains, nrOfIntraChainSSBridges, nrOfHBonds;
    u_int32_t    nrOfHBondsPerDistance[11] = {};

    protein.GetStatistics(nrOfResidues, nrOfChains, nrOfIntraChainSSBridges, nrOfHBonds, nrOfHBondsPerDistance);

    std::string                  kDSSPResidueLine("%5.5d        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

    std::vector<const MResidue*> residues;

    for (const auto &chain : protein.GetChains())
    {
        for (const auto &residue : chain->GetResidues())
        {
            residues.push_back(residue);
        }
    }
    // keep residues sorted by residue number as assigned during reading the PDB file
    sort(residues.begin(), residues.end(), [](const MResidue * a, const MResidue * b){return a->GetNumber() < b->GetNumber(); });

    const MResidue* last = nullptr;
    for (const auto &residue : residues)
    {
        // insert a break line whenever we detect missing residues
        // can be the transition to a different chain, or missing residues in the current chain
        if (last != nullptr and last->GetNumber() + 1 != residue->GetNumber())
        {
            char breaktype = ' ';
            if (last->GetChainID() != residue->GetChainID())
            {
                breaktype = '*';
            }
            char buf[1000];
            sprintf(buf, kDSSPResidueLine.c_str(), last->GetNumber() + 1, breaktype);
            dsspString += std::string(buf) + "\n";
        }
        dsspString += ResidueToDSSPLine(*residue) + "\n";
        last        = residue;
    }
    return dsspString;
}
