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

#ifndef GMX_DSSP_STRUCTURE
#define GMX_DSSP_STRUCTURE

// #include "primitives-3d.h"
#include <cmath>
#include <vector>
#include "gromacs/math/vectypes.h"

struct MAtom;
class MResidue;
class MChain;
class MProtein;

const u_int32_t kHistogramSize = 30;

#include <string>

// Code for amino acid sequences

typedef std::basic_string<u_int8_t> sequence;

// 22 real letters and 1 dummy (X is the dummy, B and Z are pseudo letters)
extern const char     kResidues[]; // = "ACDEFGHIKLMNPQRSTVWYBZX";
extern const u_int8_t kResidueNrTable[];

inline u_int8_t ResidueNr(char inAA)
{
    int result = 23;

    inAA |= 040;
    if (inAA >= 'a' and inAA <= 'z')
    {
        result = kResidueNrTable[inAA - 'a'];
    }

    return result;
}

inline bool is_gap(char aa)
{
    return aa == ' ' or aa == '.' or aa == '-';
}

sequence encode(const std::string &s);
std::string decode(const sequence &s);


// a limited set of known atoms. This is an obvious candidate for improvement
// of DSSP.
enum MAtomType
{
    kUnknownAtom,
    kHydrogen,
    // ...
    kCarbon,
    kNitrogen,
    kOxygen,
    kFluorine,
    // ...
    kPhosphorus,
    kSulfur,
    kChlorine,
    kMagnesium,
    kPotassium,
    kCalcium,
    kZinc,
    kSelenium,

    kAtomTypeCount
};

MAtomType MapElement(std::string inElement);

// for now, MAtom contains exactly what the ATOM line contains in a PDB file
struct MAtom
{
    u_int32_t        mSerial;
    std::string      mName;
    char             mAltLoc;
    std::string      mResName;
    std::string      mChainID, mAuthChainID;
    int32_t          mResSeq;
    std::string      mICode;
    MAtomType        mType;
    gmx::RVec        mLoc;
    double           mOccupancy;
    double           mTempFactor;
    std::string      mElement;
    int              mCharge;

    void             SetChainID(const std::string &inChainID){ mChainID = inChainID; }
    std::string      GetName() const              { return mName; }
    operator const gmx::RVec &() const      { return mLoc; }
    operator gmx::RVec &()            { return mLoc; }
};

enum class ResidueType : char
{
    UnknownResidue = 'U',
    Alanine        = 'A', //   ala
    Arginine       = 'R', //   arg
    Asparagine     = 'N', //   asn
    AsparticAcid   = 'D', //   asp
    Cysteine       = 'C', //   cys
    GlutamicAcid   = 'E', //   glu
    Glutamine      = 'Q', //   gln
    Glycine        = 'G', //   gly
    Histidine      = 'H', //   his
    Isoleucine     = 'I', //   ile
    Leucine        = 'L', //   leu
    Lysine         = 'K', //   lys
    Methionine     = 'M', //   met
    Phenylalanine  = 'F', //   phe
    Proline        = 'P', //   pro
    Serine         = 'S', //   ser
    Threonine      = 'T', //   thr
    Tryptophan     = 'W', //   trp
    Tyrosine       = 'Y', //   tyr
    Valine         = 'V', //   val
};

struct MResidueInfo
{
    ResidueType type;
    char        code;
    char        name[4];
};


ResidueType MapResidue(std::string inName);

struct HBond
{
    MResidue*    residue;
    double       energy;
};

enum MBridgeType
{
    btNoBridge, btParallel, btAntiParallel
};

struct MBridgeParner
{
    MResidue  *    residue;
    u_int32_t      ladder;
    bool           parallel;
};

enum MHelixFlag
{
    helixNone, helixStart, helixEnd, helixStartAndEnd, helixMiddle
};
enum class SecondaryStructure : char;

class MResidue
{
    public:
        MResidue(const MResidue &residue);
        MResidue(int32_t inNumber, char inTypeCode, MResidue* inPrevious);
        MResidue(int32_t inNumber, MResidue* inPrevious, const std::vector<MAtom> &inAtoms);

        void        SetChainID(const std::string &inChainID);
        std::string GetChainID() const { return mChainID; }

        ResidueType GetType() const { return mType; }

        const MAtom &GetCAlpha() const { return mCA; }
        const MAtom &GetC() const   { return mC; }
        const MAtom &GetN() const   { return mN; }
        const MAtom &GetO() const   { return mO; }
        const MAtom &GetH() const   { return mH; }

        double Phi() const;
        double Psi() const;
        std::tuple<double, char> Alpha() const;
        double Kappa() const;
        double TCO() const;

        double Accessibility() const      { return mAccessibility; }

        void SetSecondaryStructure(SecondaryStructure inSS) { mSecondaryStructure = inSS; }
        SecondaryStructure  GetSecondaryStructure() const
        {
            return mSecondaryStructure;
        }

        const MResidue* Next() const { return mNext; }
        const MResidue* Prev() const { return mPrev; }

        void SetPrev(MResidue* inResidue);

        void SetBetaPartner(u_int32_t n, MResidue* inResidue, u_int32_t inLadder, bool inParallel);
        MBridgeParner GetBetaPartner(u_int32_t n) const;

        void SetSheet(u_int32_t inSheet)  { mSheet = inSheet; }
        u_int32_t GetSheet() const { return mSheet; }

        bool        IsBend() const { return mBend; }
        void        SetBend(bool inBend) { mBend = inBend; }

        MHelixFlag      GetHelixFlag(u_int32_t inHelixStride) const;
        bool        IsHelixStart(u_int32_t inHelixStride) const;
        void        SetHelixFlag(u_int32_t inHelixStride, MHelixFlag inHelixFlag);

        void        SetSSBridgeNr(u_int8_t inBridgeNr);
        u_int8_t        GetSSBridgeNr() const;

        void        AddAtom(MAtom &inAtom);

        HBond*        Donor()         { return mHBondDonor; }
        HBond*        Acceptor()       { return mHBondAcceptor; }

        const HBond*    Donor() const     { return mHBondDonor; }
        const HBond*    Acceptor() const   { return mHBondAcceptor; }

        bool        ValidDistance(const MResidue &inNext) const;

        static bool      TestBond(const MResidue* a, const MResidue* b)
        {
            return a->TestBond(b);
        }

        // bridge functions
        MBridgeType      TestBridge(MResidue* inResidue) const;

        int16_t        GetSeqNumber() const    { return mSeqNumber; }
        std::string      GetInsertionCode() const  { return mInsertionCode; }

        void        SetNumber(u_int16_t inNumber)  { mNumber = inNumber; }
        u_int16_t        GetNumber() const      { return mNumber; }

        static double CalculateHBondEnergy(MResidue &inDonor, MResidue &inAcceptor);

        std::vector<MAtom> &GetSideChain()        { return mSideChain; }
        const std::vector<MAtom> &
        GetSideChain() const    { return mSideChain; }

        void        CalculateSurface(const std::vector<MResidue*> &inResidues);

        void        GetCenterAndRadius(gmx::RVec &outCenter, double &outRadius) const
        { outCenter = mCenter; outRadius = mRadius; }

        static bool      NoChainBreak(const MResidue* from, const MResidue* to);

    protected:

        double        CalculateSurface(
            const MAtom &inAtom, double inRadius,
            const std::vector<MResidue*> &inResidues);

        bool        TestBond(const MResidue* other) const;

        void        ExtendBox(const MAtom &atom, double inRadius);
        bool        AtomIntersectsBox(const MAtom &atom, double inRadius) const;

        std::string           mChainID;
        MResidue       *      mPrev;
        MResidue       *      mNext;
        int32_t               mSeqNumber, mNumber;
        std::string           mInsertionCode;
        ResidueType           mType;
        u_int8_t              mSSBridgeNr;
        double                mAccessibility;
        SecondaryStructure    mSecondaryStructure;
        MAtom                 mC, mN, mCA, mO, mH;
        HBond                 mHBondDonor[2], mHBondAcceptor[2];
        std::vector<MAtom>    mSideChain;
        MBridgeParner         mBetaPartner[2];
        u_int32_t             mSheet;
        MHelixFlag            mHelixFlags[3]; //
        bool                  mBend;
        gmx::RVec             mBox[2];        // The 3D box containing all atoms
        gmx::RVec             mCenter;        // and the 3d Sphere containing all atoms
        double                mRadius;

    private:
        MResidue &operator=(const MResidue &residue);
};

class MChain
{
    public:

        MChain(const MChain &chain);
        MChain(const std::string &inChainID) : mChainID(inChainID) {}
        ~MChain();

        MChain &operator=(const MChain &chain);

        std::string GetChainID() const { return mChainID; }
        void        SetChainID(const std::string &inChainID);

        std::string GetAuthChainID(void) const;
        void SetAuthChainID(const std::string &inAuthChainID);

        const MResidue* GetResidueBySeqNumber(u_int16_t inSeqNumber, const std::string &inInsertionCode) const;

        std::vector<MResidue*> &
        GetResidues() { return mResidues; }
        const std::vector<MResidue*> &
        GetResidues() const { return mResidues; }

        bool Empty() const { return mResidues.empty(); }

    private:
        std::string      mChainID,
                         mAuthChainID;
        std::vector<MResidue*>
        mResidues;
};

class MProtein
{
    public:
        MProtein();
        MProtein(const std::string &inID, MChain* inChain);
        ~MProtein();

        void        ReadPDB(const std::string &pdbString, bool inCAlphaOnly = false);

        const std::string &GetID() const          { return mID; }
        const std::string &GetHeader() const        { return mHeader; }
        std::string      GetCompound() const;
        std::string      GetSource() const;
        std::string      GetAuthor() const;
        const std::vector<std::string> &GetDbRef() const { return mDbRef; }

        void CalculateSecondaryStructure(bool inPreferPiHelices = true);

        void GetStatistics(u_int32_t& outNrOfResidues, u_int32_t& outNrOfChains, u_int32_t& outNrOfIntraChainSSBridges, u_int32_t& outNrOfHBonds, u_int32_t outNrOfHBondsPerDistance[11]) const;

        void GetCAlphaLocations(const std::string      &inChainID,
                                std::vector<gmx::RVec> &outPoints) const;
        gmx::RVec GetCAlphaPosition(const std::string &inChainID, int16_t inPDBResSeq) const;

        std::string GetFirstChainID() const
        {
            return mChains.front()->GetChainID();
        }

        void        SetChain(const std::string &inChainID, const MChain &inChain);

        MChain       &GetChain(const std::string &inChainID);
        const MChain &GetChain(const std::string &inChainID) const;

        const std::vector<MChain*> &GetChains() const { return mChains; }

        MResidue* GetResidue(const std::string &inChainID, u_int16_t inSeqNumber,
                             const std::string &inInsertionCode);

        const MResidue* GetResidue(const std::string &inChainID, u_int16_t inSeqNumber,
                                   const std::string &inInsertionCode) const;
        std::vector<std::string> writeSecondaryStructure();
    private:

        void AddResidue(const std::vector<MAtom> &inAtoms);

        void CalculateHBondEnergies(const std::vector<MResidue*> &inResidues);
        void CalculateAlphaHelices(const std::vector<MResidue*> &inResidues,
                                   bool                          inPreferPiHelices);
        void CalculateBetaSheets(const std::vector<MResidue*> &inResidues);
        void CalculateAccessibilities(const std::vector<MResidue*> &inResidues);

        std::string                                   mID, mHeader;

        std::vector<std::string>                      mDbRef;
        std::string                                   mCompound, mSource, mAuthor;
        std::vector<MChain*>                          mChains;
        u_int32_t                                     mResidueCount, mChainBreaks;

        std::vector<std::pair<MResidue*, MResidue*> > mSSBonds;
};

// inlines

// Write the DSSP line for a single residue
std::string ResidueToDSSPLine(const MResidue &residue);

// Write a complete DSSP file for a protein
std::string WriteDSSP(MProtein &protein);

#endif
