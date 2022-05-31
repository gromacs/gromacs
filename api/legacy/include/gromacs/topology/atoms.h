/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_TOPOLOGY_ATOMS_H
#define GMX_TOPOLOGY_ATOMS_H

#include <stdio.h>

#include <optional>
#include <vector>

#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
class ISerializer;
} // namespace gmx

/*! \brief
 * Contains information for a single particle in a PDB file.
 *
 * Currently only supports ATOM/HETATM lines, as well as anisotropy information.
 */
class PdbAtomEntry
{
public:
    //! Construct full structure without anisotropy information, bfactor or occupancy.
    PdbAtomEntry(PdbRecordType type, int pdbAtomNumber, char alternativeLocation, const std::string& atomName) :
        PdbAtomEntry(type, pdbAtomNumber, alternativeLocation, atomName, std::nullopt, std::nullopt)
    {
    }

    //! Construct full structure without anisotropy information, but with bfactor and occupancy.
    PdbAtomEntry(PdbRecordType       type,
                 int                 pdbAtomNumber,
                 char                alternativeLocation,
                 const std::string&  atomName,
                 std::optional<real> occupancy,
                 std::optional<real> bFactor) :
        PdbAtomEntry(type, pdbAtomNumber, alternativeLocation, atomName, occupancy, bFactor, std::nullopt)
    {
    }
    //! Construct full structure.
    PdbAtomEntry(PdbRecordType                      type,
                 int                                atomSerialNumber,
                 char                               alternativeLocation,
                 const std::string&                 atomName,
                 std::optional<real>                occupancy,
                 std::optional<real>                bFactor,
                 std::optional<std::array<real, 6>> anisotropy) :
        type_(type),
        atomSerialNumber_(atomSerialNumber),
        alternativeLocation_(alternativeLocation),
        atomName_(atomName),
        occupancy_(occupancy),
        bFactor_(bFactor),
        anisotropyTensor_(anisotropy)
    {
        if (atomName.size() > 6)
        {
            GMX_THROW(gmx::InconsistentInputError(
                    "Cannot have atom name with more than 6 characters"));
        }
    }
    //! Get PDB record type
    PdbRecordType type() const { return type_; }
    //! Get atom number.
    int atomSerialNumber() const { return atomSerialNumber_; }
    //! Get access to alternative location identifier.
    char altloc() const { return alternativeLocation_; }
    //! Get access to real atom name.
    const std::string& atomName() const { return atomName_; }
    //! Get access to occupancy.
    std::optional<real> occupancy() const { return occupancy_; }
    //! Get access to b factor.
    std::optional<real> bFactor() const { return bFactor_; }
    //! Get access to anisotropy values.
    std::optional<gmx::ArrayRef<const real>> anisotropy() const
    {
        return anisotropyTensor_.has_value()
                       ? std::make_optional(gmx::makeConstArrayRef(anisotropyTensor_.value()))
                       : std::nullopt;
    }

private:
    //! PDB record type
    PdbRecordType type_;
    //! PDB atom number.
    int atomSerialNumber_;
    //! Defines alternative location in PDB.
    char alternativeLocation_;
    //! The actual atom name from the pdb file.
    std::string atomName_;
    //! Occupancy field, abused for other things.
    std::optional<real> occupancy_;
    //! B-Factor field, abused for other things.
    std::optional<real> bFactor_;
    //! Tensor of anisotropy values.
    std::optional<std::array<real, 6>> anisotropyTensor_;
};

// Legacy types begin here
typedef struct t_atom
{
    real           m, q;       /* Mass and charge                      */
    real           mB, qB;     /* Mass and charge for Free Energy calc */
    unsigned short type;       /* Atom type                            */
    unsigned short typeB;      /* Atom type for Free Energy calc       */
    ParticleType   ptype;      /* Particle type                        */
    int            resind;     /* Index into resinfo (in t_atoms)      */
    int            atomnumber; /* Atomic Number or 0                   */
    char           elem[4];    /* Element name                         */
} t_atom;

typedef struct t_resinfo
{
    char**        name;     /* Pointer to the residue name          */
    int           nr;       /* Residue number                       */
    unsigned char ic;       /* Code for insertion of residues       */
    int           chainnum; /* Iincremented at TER or new chain id  */
    char          chainid;  /* Chain identifier written/read to pdb */
    char**        rtp;      /* rtp building block name (optional)   */
} t_resinfo;

typedef struct t_pdbinfo
{
    PdbRecordType type;         /* PDB record name                      */
    int           atomnr;       /* PDB atom number                      */
    char          altloc;       /* Alternate location indicator         */
    char          atomnm[6];    /* True atom name including leading spaces */
    real          occup;        /* Occupancy                            */
    real          bfac;         /* B-factor                             */
    bool          bAnisotropic; /* (an)isotropic switch                 */
    int           uij[6];       /* Anisotropic B-factor                 */
} t_pdbinfo;

//! Contains indices into group names for different groups.
using AtomGroupIndices = std::vector<int>;

typedef struct t_atoms
{
    int     nr;         /* Nr of atoms                          */
    t_atom* atom;       /* Array of atoms (dim: nr)             */
                        /* The following entries will not       */
                        /* always be used (nres==0)             */
    char*** atomname;   /* Array of pointers to atom name       */
                        /* use: (*(atomname[i]))                */
    char*** atomtype;   /* Array of pointers to atom types      */
                        /* use: (*(atomtype[i]))                */
    char*** atomtypeB;  /* Array of pointers to B atom types    */
                        /* use: (*(atomtypeB[i]))               */
    int        nres;    /* The number of resinfo entries        */
    t_resinfo* resinfo; /* Array of residue names and numbers   */
    t_pdbinfo* pdbinfo; /* PDB Information, such as aniso. Bfac */

    /* Flags that tell if properties are set for all nr atoms.
     * For B-state parameters, both haveBState and the mass/charge/type
     * flag should be TRUE.
     */
    bool haveMass;    /* Mass available                       */
    bool haveCharge;  /* Charge available                     */
    bool haveType;    /* Atom type available                  */
    bool haveBState;  /* B-state parameters available         */
    bool havePdbInfo; /* pdbinfo available                    */
} t_atoms;

#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))

void init_atom(t_atoms* at);
void done_atom(t_atoms* at);
void done_and_delete_atoms(t_atoms* atoms);

void init_t_atoms(t_atoms* atoms, int natoms, bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */

void gmx_pdbinfo_init_default(t_pdbinfo* pdbinfo);

t_atoms* copy_t_atoms(const t_atoms* src);
/* copy an atoms struct from src to a new one */

void add_t_atoms(t_atoms* atoms, int natom_extra, int nres_extra);
/* allocate extra space for more atoms and or residues */

void t_atoms_set_resinfo(t_atoms*         atoms,
                         int              atom_ind,
                         struct t_symtab* symtab,
                         const char*      resname,
                         int              resnr,
                         unsigned char    ic,
                         int              chainnum,
                         char             chainid);
/* Set the residue name, number, insertion code and chain identifier
 * of atom index atom_ind.
 */

void pr_atoms(FILE* fp, int indent, const char* title, const t_atoms* atoms, bool bShownumbers);

/*! \brief Compare information in the t_atoms data structure.
 *
 * \param[in] fp Pointer to file to write to.
 * \param[in] a1 Pointer to first data structure to compare.
 * \param[in] a2 Pointer to second data structure or nullptr.
 * \param[in] relativeTolerance Relative floating point comparison tolerance.
 * \param[in] absoluteTolerance Absolute floating point comparison tolerance.
 */
void compareAtoms(FILE* fp, const t_atoms* a1, const t_atoms* a2, real relativeTolerance, real absoluteTolerance);

/*! \brief Set mass for each atom using the atom and residue names using a database
 *
 * If atoms->haveMass = TRUE does nothing.
 * If printMissingMasss = TRUE, prints details for first 10 missing masses
 * to stderr.
 */
void atomsSetMassesBasedOnNames(t_atoms* atoms, bool printMissingMasses);

//! Deleter for t_atoms, needed until it has a proper destructor.
using AtomsDataPtr = gmx::unique_cptr<t_atoms, done_and_delete_atoms>;


#endif
