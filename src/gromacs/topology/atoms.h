/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2018, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_ATOMS_H
#define GMX_TOPOLOGY_ATOMS_H

#include <vector>
#include <algorithm>
#include <stdio.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

struct t_symtab;

/* The particle type */
enum {
    eptAtom, eptNucleus, eptShell, eptBond, eptVSite, eptNR
};

/* The particle type names */
extern const char *ptype_str[eptNR+1];

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file
 */
enum PDB_record {
    epdbATOM,   epdbHETATM, epdbANISOU, epdbCRYST1, epdbCOMPND,
    epdbMODEL,  epdbENDMDL, epdbTER,    epdbHEADER, epdbTITLE, epdbREMARK,
    epdbCONECT, epdbNR
};

/*! \brief
 * Data structure describing a single residue.
 */
struct Residue
{
    public:

    //! Constructor to explicitly initialize all fields.
    explicit Residue(char **name, int nr, unsigned char ic, int chainnum,
           char chainid, const char **rtp) :
           name_(name), nr_(nr), ic_(ic), chainnum_(chainnum), chainid_(chainid), rtp_(rtp) {}
    Residue() {}

    //! Residue name in SymbolTable.
    char **name_ = nullptr;
    //! Residue number.
    int              nr_ = -1;
    //! Code for insertion of residues.
    unsigned char    ic_ = ' ';
    //! Chain number, incremented at TER or new chain id.
    int              chainnum_ = -1;
    //! Chain identifier written/read to pdb.
    char             chainid_ = ' ';
    //! Optional rtp building block name.
    const char **rtp_ = nullptr;
};

/*! \brief
 * Data structure describing a single PDB entry.
 */
struct PdbEntry
{
    //! Constructor to initialize all fields.
    explicit PdbEntry(int type, int atomnr, char altloc, std::string atomnm, real occup,
            real bfac) :
        type_(type), atomnr_(atomnr), altloc_(altloc), atomnm_(atomnm), occup_(occup), bfac_(bfac),
        haveAnisotropic_(false), isSet_(true) {}
    //! Constructor for setting empty object that wont be used.
    explicit PdbEntry() :
        type_(-1), atomnr_(-1), altloc_(' '), occup_(0.0), bfac_(0.0), haveAnisotropic_(false),
        uij_({0}), isSet_(false) {}

    //! PDB record name.
    int      type_;
    //! PDB atom number.
    int      atomnr_;
    //! Alternate location indicator.
    char     altloc_;
    //! True atom name including leading spaces.
    std::string atomnm_;
    //! Occupancy. Can be fudged to contain other entries for special file formats.
    real     occup_;
    //! B-factor. Can be fudged to contain other entries for special file formats.
    real     bfac_;
    //! (an)isotropic switch.
    bool haveAnisotropic_;
    //! Anisotropic B-factor
    std::vector<int>      uij_;
    //! Has the information been set?
    bool isSet_;
};

/*! \brief
 * Structure for atomic data.
 *
 * Bundles together names and atom entries.
 * Blabla add more info here.
 */
struct AtomInfo
{
    //! Atom mass.
    real m_ = 0.0;
    //! Atom charge.
    real q_ = 0.0;
    //! Atom mass for free energy calculation.
    real mB_ = 0.0;
    //! Atom charge for free energy calculation.
    real qB_ = 0.0;
    //! Atom type. 
    unsigned short type_ = 0;
    //! Atom type for free energy calculation.
    unsigned short typeB_ = 0;
    //! Particle type.
    int            ptype_ = -1;
    //! Index into ResidueInformation.
    int            resind_ = -1;
    //! Atomic number or zero. 
    int            atomnumber_ = 0;
    //! Element name.
    std::string elem_;
    //! Has the atom charge been set?
    bool haveCharge_ = false;
    //! Has the atom mass been set?
    bool haveMass_ = false;
    //! Has the atom type been set?
    bool haveType_ = false;
    //! Have the B state parameters been set?
    bool haveBstate_ = false;
    
    //! Atom name for entry.
    char **atomname = nullptr;
    //! Atom type for entry.
    char **atomtype = nullptr;
    //! Atom B state type for entry.
    char **atomtypeB = nullptr;

    //! Mass available
    bool haveMass() const { return haveMass_; }
    //! Charge available
    bool haveCharge() const { return haveCharge_; }
    //! Atomname is set.
    bool haveAtomname() const { return atomname != nullptr; }
    //! Atom type available
    bool haveType() const { return atomtype != nullptr && haveType_; }
    //! B-state parameters available
    bool haveBState() const { return atomtypeB != nullptr && haveMass() && haveCharge() && haveType(); }
};

//! Convenience function to check if all entries have mass.
inline bool allAtomsHaveMass(gmx::ArrayRef<const AtomInfo> atoms)
{
    if (atoms.empty())
    {
        return false;
    }
    return std::all_of(atoms.begin(), atoms.end(),
                       [](AtomInfo atom)
                       {return atom.haveMass();});
}

//! Convenience function to check if all entries have charge.
inline bool allAtomsHaveCharge(gmx::ArrayRef<const AtomInfo> atoms)
{
    if (atoms.empty())
    {
        return false;
    }
    return std::all_of(atoms.begin(), atoms.end(),
                       [](AtomInfo atom)
                       {return atom.haveCharge();});
}

//! Convenience function to check if all entries have atomnames set.
inline bool allAtomsHaveAtomname(gmx::ArrayRef<const AtomInfo> atoms)
{
    if (atoms.empty())
    {
        return false;
    }
    return std::all_of(atoms.begin(), atoms.end(),
                       [](AtomInfo atom)
                       {return atom.haveAtomname();});
}


//! Convenience function to check if all entries have a type.
inline bool allAtomsHaveType(gmx::ArrayRef<const AtomInfo> atoms)
{
    if (atoms.empty())
    {
        return false;
    }
    return std::all_of(atoms.begin(), atoms.end(),
                       [](AtomInfo atom)
                       {return atom.haveType();});
}

//! Convenience function to check if all entries have B state.
inline bool allAtomsHaveBstate(gmx::ArrayRef<const AtomInfo> atoms)
{
    if (atoms.empty())
    {
        return false;
    }
    return std::all_of(atoms.begin(), atoms.end(),
                       [](AtomInfo atom)
                       {return atom.haveBState();});
}

//! Convenience function to check if all entries have PDB information.
inline bool allAtomsHavePdbInfo(gmx::ArrayRef<const PdbEntry> pdb)
{
    if (pdb.empty())
    {
        return false;
    }
    return std::all_of(pdb.begin(), pdb.end(),
                       [](PdbEntry entry)
                       {return entry.isSet_;});
}

typedef struct t_grps
{
    int   nr;                   /* Number of different groups           */
    int  *nm_ind;               /* Index in the group names             */
} t_grps;

typedef struct t_atomtypes
{
    int           nr;           /* number of atomtypes                          */
    int          *atomnumber;   /* Atomic number, used for QM/MM                */
} t_atomtypes;

#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))

void init_atomtypes(t_atomtypes *at);
void done_atomtypes(t_atomtypes *at);

void printAtoms(FILE *fp,
                int indent,
                const char *title,
                gmx::ArrayRef<const AtomInfo> atoms,
                gmx_bool bShownumbers);

void printResidues(FILE *fp,
                   int indent,
                   const char *title,
                   gmx::ArrayRef<const Residue> resinfo,
                   gmx_bool bShowNumbers);

void pr_atomtypes(FILE *fp,
                  int indent,
                  const char *title,
                  const t_atomtypes *atomtypes,
                  gmx_bool bShowNumbers);

void compareAtomInfo(FILE *fp,
                     gmx::ArrayRef<const AtomInfo> a1,
                     gmx::ArrayRef<const AtomInfo> a2,
                     real ftol,
                     real abstol);

void compareAtomFEPData(FILE *fp,
                        gmx::ArrayRef<const AtomInfo> atoms,
                        real ftol,
                        real abstol);

/*! \brief
 * Set mass for each atom using the atom and residue names using a database
 *
 * If all atoms already have masses, nothing is done.
 * If \p printMissingMasses is true, prints details for first 10 missing entries.
 *
 * \param[in] atoms The atom information.
 * \param[in] resinfo The residue information.
 * \param[in] printMissingMasses Whever we will print messages about missing masses or not.
 */
void atomsSetMassesBasedOnNames(gmx::ArrayRef<AtomInfo> atoms,
                                gmx::ArrayRef<const Residue> resinfo,
                                bool printMissingMasses);

/*! \brief
 * Convenience struct that contains all information for atoms.
 */
struct AtomResiduePdb {
    //! The normal atoms information.
    std::vector<AtomInfo> atoms;
    //! The separate residue information.
    std::vector<Residue> resinfo;
    //! The optional pdb information.
    std::vector<PdbEntry> pdb;
};



//! Convenience type to get new handle to AtomInfo.
using AtomsDataPtr = std::unique_ptr<std::vector<AtomInfo>>;
//! Convenience type to get new handle to all atom information.
using AtomResiduePdbDataPtr = std::unique_ptr<AtomResiduePdb>;


#endif
