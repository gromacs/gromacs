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

#include <cstdio>

#include <optional>
#include <vector>

#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

/* The particle type */
enum class ParticleType : int;
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


//! Typedef for names that might be set.
using NameHolder = std::optional<StringTableEntry>;

//! Template wrapper struct for particle FEP state values
template<typename T>
struct FEPStateValue
{
    //! Empty object.
    FEPStateValue() = default;
    //! Construct without FEP state changes.
    explicit FEPStateValue(T value) : storage_({ value, T() }), haveBState_{ false } {}
    //! Construct with FEP state changes.
    FEPStateValue(T valueA, T valueB) : storage_({ valueA, valueB }), haveBState_{ true } {}
    //! Internal storage.
    std::array<T, 2> storage_ = {};
    //! Whether this value has B state set or not.
    bool haveBState_ = false;
};

//! Typedef for particle mass.
using ParticleMass = FEPStateValue<real>;
//! Typedef for particle charge.
using ParticleCharge = FEPStateValue<real>;
//! Typedef for particle type value.
using ParticleTypeValue = FEPStateValue<unsigned short>;
//! Typedef for particle type name.
using ParticleTypeName = FEPStateValue<NameHolder>;

//! Single particle in a simulation.
class SimulationParticle
{
public:
    //! Write info to serializer.
    void serializeParticle(gmx::ISerializer* serializer);
    //! Access mass. A state.
    real m() const { return mass_.storage_[0]; }
    //! Access charge. A state.
    real q() const { return charge_.storage_[0]; }
    //! Access atom type. A state.
    unsigned short type() const { return particleTypeValue_.storage_[0]; }
    //! Access mass. B state.
    real mB() const { return haveBStateForAll() ? mass_.storage_[1] : mass_.storage_[0]; }
    //! Access charge. B state.
    real qB() const { return haveBStateForAll() ? charge_.storage_[1] : charge_.storage_[0]; }
    //! Access atom type. B state.
    unsigned short typeB() const
    {
        return haveBStateForAll() ? particleTypeValue_.storage_[1] : particleTypeValue_.storage_[0];
    }
    //! Access particle name.
    std::string particleName() const
    {
        GMX_ASSERT(haveParticleName(), "Can not access uninitialized element");
        return *particleName_.value();
    }
    //! Access type name for state A.
    std::string particleTypeNameA() const
    {
        GMX_ASSERT(haveParticleTypeName(), "Can not access uninitialized element");
        return *particleTypeName_.storage_[0].value();
    }
    //! Access type name for state B if it exists.
    std::string particleTypeNameB() const
    {
        GMX_ASSERT(haveParticleTypeName(), "Can not access uninitialized element");
        const auto entry =
                haveBStateForAll() ? particleTypeName_.storage_[1] : particleTypeName_.storage_[0];
        GMX_ASSERT(entry.has_value(), "Can not access uninitialized element");
        return *entry.value();
    }

    //! Access particle type.
    ParticleType ptype() const { return particleType_; }
    //! Access residue index.
    gmx::index resind() const { return residueIndex_; }
    //! Access atomic number.
    int atomnumber() const { return atomicNumber_; }
    //! Access element name.
    std::string elem() const { return elementName_; }
    //! Do we have mass?
    bool haveMass() const { return haveMass_; }
    //! Do we have charge?
    bool haveCharge() const { return haveCharge_; }
    //! Do we have type?
    bool haveType() const { return haveType_; }
    //! Do we have particle name set?
    bool haveParticleName() const { return haveParticleName_; }
    //! Do we have the particle type name set.
    bool haveParticleTypeName() const { return haveParticleTypeName_; }
    //! Do we have B state?
    bool haveBStateForAll() const { return haveBStateForAll_; }

    //! Constructor with complete information. A and B states are equivalent.
    SimulationParticle(const std::optional<ParticleMass>      mass,
                       const std::optional<ParticleCharge>    charge,
                       const std::optional<ParticleTypeValue> particleTypeValue,
                       const std::optional<ParticleTypeName>  particleTypeName,
                       NameHolder                             particleName,
                       ParticleType                           particleType,
                       gmx::index                             residueIndex,
                       int                                    atomicNumber,
                       const std::string&                     elementName) :
        mass_(mass.has_value() ? *mass : ParticleMass()),
        charge_(charge.has_value() ? *charge : ParticleCharge()),
        particleTypeValue_(particleTypeValue.has_value() ? *particleTypeValue : ParticleTypeValue()),
        particleTypeName_(particleTypeName.has_value() ? *particleTypeName : ParticleTypeName()),
        particleName_(particleName),
        particleType_(particleType),
        residueIndex_(residueIndex),
        atomicNumber_(atomicNumber),
        elementName_(elementName),
        haveMass_(mass.has_value()),
        haveCharge_(charge.has_value()),
        haveType_(particleTypeValue.has_value()),
        haveParticleName_(particleName.has_value()),
        haveParticleTypeName_(particleTypeName.has_value()),
        haveBStateForAll_(mass_.haveBState_ && charge_.haveBState_ && particleTypeValue_.haveBState_
                          && particleTypeName_.haveBState_)
    {
        GMX_ASSERT(elementName.length() <= 4, "Element name can only be three characters");
    }
    //! Construct new datastructure from deserialization.
    SimulationParticle(gmx::ISerializer* serializer, const StringTable& table);

    //! Copy constructor.
    SimulationParticle(const SimulationParticle&) = default;
    //! Copy assignment.
    SimulationParticle& operator=(const SimulationParticle&) = default;
    //! Default move constructor.
    SimulationParticle(SimulationParticle&&) = default;
    //! Default move assignment.
    SimulationParticle& operator=(SimulationParticle&&) = default;

private:
    //! Mass of the particle. A and B state.
    ParticleMass mass_;
    //! Charge of the particle. A and B state.
    ParticleCharge charge_;
    //! Atom type. A and B state.
    ParticleTypeValue particleTypeValue_;
    //! Atom type name. A and B state.
    ParticleTypeName particleTypeName_;
    //! Atom name.
    NameHolder particleName_;
    //! Type of the particle.
    ParticleType particleType_;
    //! Residue this atoms is part of.
    int64_t residueIndex_;
    //! Atomic Number or 0.
    int atomicNumber_;
    //! Name of the element if applicable.
    std::string elementName_;
    //! If we have mass for the particle.
    bool haveMass_;
    //! If we have charge for the particle.
    bool haveCharge_;
    //! If the particle type is set.
    bool haveType_;
    //! If the particle name is set.
    bool haveParticleName_;
    //! If the particle type name is set.
    bool haveParticleTypeName_;
    //! If all fields have B state set.
    bool haveBStateForAll_;
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

//! Returns whether the charge of the atom is perturbed
inline bool atomHasPerturbedCharge(const t_atom& atom)
{
    return atom.qB != atom.q;
}

//! Returns whether the atom is perturbed, i.e. mass, charge or atom type
inline bool PERTURBED(const t_atom& atom)
{
    return atom.mB != atom.m || atomHasPerturbedCharge(atom) || atom.typeB != atom.type;
}

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
 * If printMissingMasses = TRUE, prints details for first 10 missing masses
 * to stderr.
 */
void atomsSetMassesBasedOnNames(t_atoms* atoms, bool printMissingMasses);

//! Deleter for t_atoms, needed until it has a proper destructor.
using AtomsDataPtr = gmx::unique_cptr<t_atoms, done_and_delete_atoms>;


#endif
