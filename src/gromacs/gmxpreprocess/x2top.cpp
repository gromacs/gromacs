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
#include "gmxpre.h"

#include "x2top.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/gmxpreprocess/extendableforcefield.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/nm2type.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vecdump.h"
#include "gromacs/utility/vectypes.h"

#include "hackblock.h"

struct gmx_output_env_t;

namespace gmx
{

static bool is_bond(int nnm, t_nm2type nmt[], const char* ai, const char* aj, real blen)
{
    int i, j;

    for (i = 0; (i < nnm); i++)
    {
        for (j = 0; (j < nmt[i].nbonds); j++)
        {
            if ((((equalCaseInsensitive(ai, nmt[i].elem, 1)) && (equalCaseInsensitive(aj, nmt[i].bond[j], 1)))
                 || ((equalCaseInsensitive(ai, nmt[i].bond[j], 1))
                     && (equalCaseInsensitive(aj, nmt[i].elem, 1))))
                && (std::fabs(blen - nmt[i].blen[j]) <= 0.1 * nmt[i].blen[j]))
            {
                return TRUE;
            }
        }
    }
    return FALSE;
}

/*! \brief Return the PBC type to use for distance calculations on the input structure.
 *
 * Structures without a box (PbcType::No) get the automatic detection of set_pbc(). */
static PbcType pbcTypeForCalculation(PbcType pbcType)
{
    return (pbcType == PbcType::No) ? PbcType::Unset : pbcType;
}

static void mk_bonds(int                 nnm,
                     t_nm2type           nmt[],
                     t_atoms*            atoms,
                     const rvec          x[],
                     InteractionsOfType* bond,
                     int                 nbond[],
                     bool                bPBC,
                     PbcType             pbcType,
                     matrix              box)
{
    int   i, j;
    t_pbc pbc;
    rvec  dx;
    real  dx2;

    std::array<real, MAXFORCEPARAM> forceParam = { 0.0 };
    if (bPBC)
    {
        set_pbc(&pbc, pbcTypeForCalculation(pbcType), box);
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        for (j = i + 1; (j < atoms->nr); j++)
        {
            if (bPBC)
            {
                pbc_dx(&pbc, x[i], x[j], dx);
            }
            else
            {
                rvec_sub(x[i], x[j], dx);
            }

            dx2 = iprod(dx, dx);
            if (is_bond(nnm, nmt, *atoms->atomname[i], *atoms->atomname[j], std::sqrt(dx2)))
            {
                forceParam[0]             = std::sqrt(dx2);
                std::vector<int> atomPair = { i, j };
                add_param_to_list(bond, InteractionOfType(atomPair, forceParam));
                nbond[i]++;
                nbond[j]++;
            }
        }
    }
}

static int* set_cgnr(t_atoms* atoms, bool bUsePDBcharge, real* qtot, real* mtot)
{
    int    i, n = 1;
    int*   cgnr;
    double qt = 0;

    *qtot = *mtot = 0;
    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->pdbinfo && bUsePDBcharge)
        {
            atoms->atom[i].q = atoms->pdbinfo[i].bfac;
        }
        qt += atoms->atom[i].q;
        *qtot += atoms->atom[i].q;
        *mtot += atoms->atom[i].m;
        cgnr[i] = n;
        if (is_int(qt))
        {
            n++;
            qt = 0;
        }
    }
    return cgnr;
}

static void set_atom_type(PreprocessingAtomTypes* atypes,
                          t_atoms*                atoms,
                          InteractionsOfType*     bonds,
                          int*                    nbonds,
                          int                     nnm,
                          t_nm2type               nm2t[],
                          const MDLogger&         logger)
{
    int nresolved;

    snew(atoms->atomtype, atoms->nr);
    nresolved = nm2type(nnm, nm2t, atoms, atypes, nbonds, bonds);
    if (nresolved != atoms->nr)
    {
        gmx_fatal(FARGS, "Could only find a forcefield type for %d out of %d atoms", nresolved, atoms->nr);
    }

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("There are %zu different atom types in your sample", atypes->size());
}

static void lo_set_force_const(InteractionsOfType* plist, real c[], int nrfp, bool bRound, bool bDih, bool bParam)
{
    double cc;
    char   buf[32];

    for (auto& param : plist->interactionTypes)
    {
        if (!bParam)
        {
            for (int j = 0; j < nrfp; j++)
            {
                c[j] = NOTSET;
            }
        }
        else
        {
            if (bRound)
            {
                sprintf(buf, "%.2e", param.c0());
                sscanf(buf, "%lf", &cc);
                c[0] = cc;
            }
            else
            {
                c[0] = param.c0();
            }
            if (bDih)
            {
                c[0] *= c[2];
                c[0] = (static_cast<int>(c[0] + 3600)) % 360;
                if (c[0] > 180)
                {
                    c[0] -= 360;
                }
                /* To put the minimum at the current angle rather than the maximum */
                c[0] += 180;
            }
        }
        GMX_ASSERT(nrfp <= MAXFORCEPARAM / 2, "Only 6 parameters may be used for an interaction");
        std::array<real, MAXFORCEPARAM> forceParam;
        for (int j = 0; (j < nrfp); j++)
        {
            forceParam[j]        = c[j];
            forceParam[nrfp + j] = c[j];
        }
        param = InteractionOfType(param.atoms(), forceParam);
    }
}

/*! \brief Return whether the atom is a hydrogen (or deuterium).
 *
 * The mass that the n2t table gave the atom is used first, since that table is made for the
 * structure at hand (an atom named "Ho" is then a hydroxyl hydrogen when the table says 1.008 u,
 * not holmium). Without a mass the element column of a PDB file is used, and as a last resort
 * the atom name: its first letter (skipping leading blanks and digits, as in "1HB") has to be an H. */
static bool isHydrogenAtom(const t_atoms* atoms, int atomIndex)
{
    const t_atom& atom = atoms->atom[atomIndex];
    // Hydrogen is 1.008 u and deuterium 2.014 u, all other elements are heavier than 4 u
    constexpr real c_maxHydrogenMass = 2.5;
    if (atom.m > 0)
    {
        return (atom.m < c_maxHydrogenMass);
    }
    if (atom.elem[0] != '\0')
    {
        const int element = std::toupper(static_cast<unsigned char>(atom.elem[0]));
        return (atom.elem[1] == '\0' && (element == 'H' || element == 'D'));
    }
    const char* name = *atoms->atomname[atomIndex];
    while (*name != '\0' && !std::isalpha(static_cast<unsigned char>(*name)))
    {
        ++name;
    }
    return (*name == 'H' || *name == 'h');
}

static void set_force_const(EnumerationArray<InteractionFunction, InteractionsOfType>& plist,
                            const t_atoms*                                             atoms,
                            real                                                       kb,
                            real                                                       kt,
                            real                                                       ktHydrogen,
                            real                                                       kp,
                            bool                                                       bRound,
                            bool                                                       bParam)
{
    real c[MAXFORCEPARAM];

    c[0] = 0;
    c[1] = kb;
    lo_set_force_const(&plist[InteractionFunction::Bonds], c, 2, bRound, FALSE, bParam);
    c[1] = kt;
    lo_set_force_const(&plist[InteractionFunction::Angles], c, 2, bRound, FALSE, bParam);
    if (bParam && ktHydrogen >= 0)
    {
        // Angles that involve a hydrogen get their own force constant
        const int nrfp = interaction_function[InteractionFunction::Angles].nrfpA;
        for (auto& angle : plist[InteractionFunction::Angles].interactionTypes)
        {
            const auto angleAtoms = angle.atoms();
            if (std::any_of(angleAtoms.begin(),
                            angleAtoms.end(),
                            [atoms](int atomIndex) { return isHydrogenAtom(atoms, atomIndex); }))
            {
                angle.setForceParameter(1, ktHydrogen);
                if (nrfp > 1)
                {
                    angle.setForceParameter(nrfp + 1, ktHydrogen);
                }
            }
        }
    }
    c[1] = kp;
    c[2] = 3;
    lo_set_force_const(&plist[InteractionFunction::ProperDihedrals], c, 3, bRound, TRUE, bParam);
}

static void calc_angles_dihs(InteractionsOfType* ang,
                             InteractionsOfType* dih,
                             const rvec          x[],
                             bool                bPBC,
                             PbcType             pbcType,
                             matrix              box)
{
    int   t1, t2, t3;
    rvec  r_ij, r_kj, r_kl, m, n;
    real  costh;
    t_pbc pbc;

    if (bPBC)
    {
        set_pbc(&pbc, pbcTypeForCalculation(pbcType), box);
    }
    for (auto& angle : ang->interactionTypes)
    {
        int  ai = angle.ai();
        int  aj = angle.aj();
        int  ak = angle.ak();
        real th = c_rad2Deg
                  * bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : nullptr, r_ij, r_kj, &costh, &t1, &t2);
        angle.setForceParameter(0, th);
    }
    for (auto dihedral : dih->interactionTypes)
    {
        int  ai = dihedral.ai();
        int  aj = dihedral.aj();
        int  ak = dihedral.ak();
        int  al = dihedral.al();
        real ph =
                c_rad2Deg
                * dih_angle(
                        x[ai], x[aj], x[ak], x[al], bPBC ? &pbc : nullptr, r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);
        dihedral.setForceParameter(0, ph);
    }
}

static void dump_hybridization(FILE* fp, t_atoms* atoms, int nbonds[])
{
    for (int i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "Atom %5s has %1d bonds\n", *atoms->atomname[i], nbonds[i]);
    }
}

static void print_pl(FILE*                                                            fp,
                     const EnumerationArray<InteractionFunction, InteractionsOfType>& plist,
                     InteractionFunction                                              ftp,
                     const char*                                                      name,
                     char***                                                          atomname)
{
    if (!plist[ftp].interactionTypes.empty())
    {
        fprintf(fp, "\n");
        fprintf(fp, "[ %s ]\n", name);
        int nrfp = interaction_function[ftp].nrfpA;
        for (const auto& param : plist[ftp].interactionTypes)
        {
            ArrayRef<const int>  atoms      = param.atoms();
            ArrayRef<const real> forceParam = param.forceParam();
            for (const auto& atom : atoms)
            {
                fprintf(fp, "  %5s", *atomname[atom]);
            }
            for (int j = 0; (j < nrfp); j++)
            {
                if (forceParam[j] != NOTSET)
                {
                    fprintf(fp, "  %10.3e", forceParam[j]);
                }
            }
            fprintf(fp, "\n");
        }
    }
}

/*! \brief Return the name of the atom type of atom \p index, or the atom name if it has none. */
static std::string atomTypeName(PreprocessingAtomTypes* atypes, const t_atoms* systemAtoms, int index)
{
    auto typeName = atypes->atomNameFromAtomType(systemAtoms->atom[index].type);
    if (typeName.has_value())
    {
        return *typeName;
    }
    return std::string(*systemAtoms->atomname[index]);
}

/*! \brief Return the force constant of a bond: parameter 1, else 3, else the first one that is set. */
static real minffBondForceConstant(const InteractionOfType& bond, real fallback)
{
    const auto params = bond.forceParam();
    if (params.size() > 1 && params[1] != NOTSET && params[1] != 0)
    {
        return params[1];
    }
    if (params.size() > 3 && params[3] != NOTSET && params[3] != 0)
    {
        return params[3];
    }
    for (std::size_t i = 1; i < params.size(); i++)
    {
        if (params[i] != NOTSET && params[i] != 0)
        {
            return params[i];
        }
    }
    return fallback;
}

/*! \brief Return the force constant of an angle: parameter 1, else 3, else the first one that differs from the angle. */
static real minffAngleForceConstant(const InteractionOfType& angle)
{
    const auto params = angle.forceParam();
    if (params.size() > 1 && params[1] != NOTSET && params[1] != 0)
    {
        return params[1];
    }
    if (params.size() > 3 && params[3] != NOTSET && params[3] != 0)
    {
        return params[3];
    }
    for (const real value : params)
    {
        if (value != NOTSET && value != 0 && value != params[0])
        {
            return value;
        }
    }
    return 0;
}

static void writeMinffBondsSection(FILE*                     fp,
                                   const InteractionsOfType& bonds,
                                   const rvec                x[],
                                   bool                      bPBC,
                                   PbcType                   pbcType,
                                   const matrix              box,
                                   const t_atoms*            systemAtoms,
                                   PreprocessingAtomTypes*   atypes,
                                   real                      defaultBondForceConstant,
                                   bool                      bExplicitForceConstant,
                                   std::optional<real>       hydrogenBondOverrideDistance,
                                   std::optional<real>       bondLengthOverrideDistance)
{
    if (bonds.size() == 0)
    {
        return;
    }

    t_pbc pbc;
    if (bPBC)
    {
        set_pbc(&pbc, pbcTypeForCalculation(pbcType), box);
    }

    fprintf(fp, "[ bonds ]\n");
    fprintf(fp, "; i     j       funct");
    if (bExplicitForceConstant)
    {
        fprintf(fp, "     length        k_b");
    }
    fprintf(fp, "\n");
    for (const auto& bond : bonds.interactionTypes)
    {
        const int ai = bond.ai();
        const int aj = bond.aj();

        if (!isHydrogenAtom(systemAtoms, ai) && !isHydrogenAtom(systemAtoms, aj))
        {
            continue;
        }

        const auto params           = bond.forceParam();
        real       measuredDistance = 0;
        {
            rvec dx;
            if (bPBC)
            {
                pbc_dx(&pbc, x[ai], x[aj], dx);
            }
            else
            {
                rvec_sub(x[ai], x[aj], dx);
            }
            measuredDistance = std::sqrt(iprod(dx, dx));
        }

        real eqDistance = 0;
        if (hydrogenBondOverrideDistance.has_value()
            && (isHydrogenAtom(systemAtoms, ai) || isHydrogenAtom(systemAtoms, aj)))
        {
            eqDistance = hydrogenBondOverrideDistance.value();
        }
        else if (bondLengthOverrideDistance.has_value())
        {
            eqDistance = bondLengthOverrideDistance.value();
        }
        else if (!params.empty() && params[0] != NOTSET && params[0] != 0)
        {
            eqDistance = params[0];
        }
        if (eqDistance <= 0 && bExplicitForceConstant && measuredDistance > 0)
        {
            eqDistance = measuredDistance;
        }

        const std::string typeNameA = atomTypeName(atypes, systemAtoms, ai);
        const std::string typeNameB = atomTypeName(atypes, systemAtoms, aj);

        if (bExplicitForceConstant)
        {
            fprintf(fp,
                    "%6d %6d %6d %10.4f %12.1f     ; %s-%s\n",
                    ai + 1,
                    aj + 1,
                    1,
                    eqDistance,
                    minffBondForceConstant(bond, defaultBondForceConstant),
                    typeNameA.c_str(),
                    typeNameB.c_str());
        }
        else
        {
            fprintf(fp, "%6d %6d %6d     ; %s-%s\n", ai + 1, aj + 1, 1, typeNameA.c_str(), typeNameB.c_str());
        }
    }
    fprintf(fp, "\n");
}

static void writeMinffAnglesSection(FILE*                     fp,
                                    const InteractionsOfType& angles,
                                    const rvec                x[],
                                    bool                      bPBC,
                                    PbcType                   pbcType,
                                    const matrix              box,
                                    const t_atoms*            systemAtoms,
                                    PreprocessingAtomTypes*   atypes,
                                    std::optional<real>       hydrogenAngleOverrideDegrees)
{
    if (angles.size() == 0)
    {
        return;
    }

    t_pbc pbc;
    if (bPBC)
    {
        set_pbc(&pbc, pbcTypeForCalculation(pbcType), box);
    }

    fprintf(fp, "[ angles ]\n");
    fprintf(fp, "; i    j   k   type   theta   force.c.\n");
    for (const auto& angle : angles.interactionTypes)
    {
        const int ai = angle.ai();
        const int aj = angle.aj();
        const int ak = angle.ak();

        rvec r_ij, r_kj;
        real costh;
        int  t1, t2;
        real thetaDegrees =
                c_rad2Deg
                * bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : nullptr, r_ij, r_kj, &costh, &t1, &t2);
        if (hydrogenAngleOverrideDegrees.has_value()
            && (isHydrogenAtom(systemAtoms, ai) || isHydrogenAtom(systemAtoms, aj)
                || isHydrogenAtom(systemAtoms, ak)))
        {
            thetaDegrees = hydrogenAngleOverrideDegrees.value();
        }

        const std::string typeNameA = atomTypeName(atypes, systemAtoms, ai);
        const std::string typeNameB = atomTypeName(atypes, systemAtoms, aj);
        const std::string typeNameC = atomTypeName(atypes, systemAtoms, ak);

        fprintf(fp,
                "%6d %6d %6d %6d %10.2f %10.2f     ; %s-%s-%s\n",
                ai + 1,
                aj + 1,
                ak + 1,
                1,
                thetaDegrees,
                minffAngleForceConstant(angle),
                typeNameA.c_str(),
                typeNameB.c_str(),
                typeNameC.c_str());
    }
    fprintf(fp, "\n");
}

static void writeMinffTopology(FILE*                   fp,
                               const char*             molname,
                               int                     nrexcl,
                               t_atoms*                atoms,
                               PreprocessingAtomTypes* atypes,
                               const EnumerationArray<InteractionFunction, InteractionsOfType>& plist,
                               t_excls*            excls,
                               const rvec          x[],
                               bool                bPBC,
                               PbcType             pbcType,
                               const matrix        box,
                               real                kb,
                               bool                bExplicitBondConstant,
                               std::optional<real> hydrogenBondOverrideDistance,
                               std::optional<real> bondLengthOverrideDistance,
                               std::optional<real> hydrogenAngleOverrideDegrees)
{
    fprintf(fp, "[ moleculetype ]\n");
    fprintf(fp, "; %-15s %5s\n", "Name", "nrexcl");
    fprintf(fp, "%-15s %5d\n\n", molname, nrexcl);

    print_atoms(fp, atypes, atoms, FALSE);

    writeMinffBondsSection(fp,
                           plist[InteractionFunction::Bonds],
                           x,
                           bPBC,
                           pbcType,
                           box,
                           atoms,
                           atypes,
                           kb,
                           bExplicitBondConstant,
                           hydrogenBondOverrideDistance,
                           bondLengthOverrideDistance);

    print_excl(fp, atoms->nr, excls);

    writeMinffAnglesSection(
            fp, plist[InteractionFunction::Angles], x, bPBC, pbcType, box, atoms, atypes, hydrogenAngleOverrideDegrees);
}

/*! \brief Return whether the name of an atom type starts with \p prefix, ignoring case. */
static bool typeNameStartsWith(const std::string& typeName, const char* prefix)
{
    return gmx_strncasecmp(typeName.c_str(), prefix, std::strlen(prefix)) == 0;
}

/*! \brief Return the formal valence of a neighbor of an oxygen from its MINFF atom type name.
 *
 * Types that do not take part in the charge smearing give 0. */
static int formalValenceOfType(const std::string& typeName)
{
    // The first matching prefix counts, so the specific iron types come before plain Fe
    static const std::array<std::pair<const char*, int>, 16> c_valences = { { { "si", 4 },
                                                                              { "al", 3 },
                                                                              { "fe2", 2 },
                                                                              { "feo2", 2 },
                                                                              { "fet2", 2 },
                                                                              { "fee2", 2 },
                                                                              { "fe3", 3 },
                                                                              { "feo3", 3 },
                                                                              { "fet3", 3 },
                                                                              { "fee3", 3 },
                                                                              { "fe", 3 },
                                                                              { "ti", 4 },
                                                                              { "li", 1 },
                                                                              { "mg", 2 },
                                                                              { "ca", 2 },
                                                                              { "h", 1 } } };
    for (const auto& entry : c_valences)
    {
        if (typeNameStartsWith(typeName, entry.first))
        {
            return entry.second;
        }
    }
    return 0;
}

/*! \brief Smear the charge of isomorphic substitutions over the oxygens and fluorines.
 *
 * The charges of all other atoms are kept as they are given by the n2t table. The charge of
 * an oxygen (type name starting with O) becomes
 *
 *   q_O = -2 + sum_j (Z_j - q_j) / CN_j
 *
 * where the sum runs over the bonded neighbors j, Z_j is the formal valence of neighbor j,
 * q_j its charge, and CN_j its number of bonds. A structural fluorine (type name starting
 * with Fs) gets the same with -1 instead of -2. This is the charge smearing described in
 * Tournassat et al., Clays and Clay Minerals 64, 374-388 (2016), doi:10.1346/CCMN.2016.0640403
 * and Lammers et al., J. Colloid Interface Sci. 2017, doi:10.1016/j.jcis.2016.11.084, as also
 * used in the atom and atomipy libraries. Finally, if the total charge is not an integer, all
 * oxygens are shifted by the same amount so that it becomes the nearest integer. */
static void smearCharges(t_atoms*                  atoms,
                         const InteractionsOfType& bonds,
                         PreprocessingAtomTypes*   atypes,
                         const MDLogger&           logger)
{
    const int numAtoms = atoms->nr;

    // Round to 1e-6 so that the printed charges do not depend on the last bit of the
    // floating-point arithmetic, which varies between compilers and platforms.
    const auto roundCharge = [](const double q)
    { return static_cast<real>(std::round(q * 1e6) / 1e6); };

    std::vector<std::vector<int>> neighbors(numAtoms);
    for (const auto& bond : bonds.interactionTypes)
    {
        neighbors[bond.ai()].push_back(bond.aj());
        neighbors[bond.aj()].push_back(bond.ai());
    }
    std::vector<std::string> typeNames(numAtoms);
    std::vector<real>        tableCharges(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        typeNames[i]    = atomTypeName(atypes, atoms, i);
        tableCharges[i] = atoms->atom[i].q;
    }

    int numSmeared = 0;
    for (int i = 0; i < numAtoms; i++)
    {
        real baseCharge;
        if (typeNameStartsWith(typeNames[i], "o"))
        {
            baseCharge = -2;
        }
        else if (typeNameStartsWith(typeNames[i], "fs"))
        {
            baseCharge = -1;
        }
        else
        {
            continue;
        }
        if (neighbors[i].empty())
        {
            continue;
        }
        double sum = 0;
        for (const int j : neighbors[i])
        {
            const int coordination = std::max<int>(1, neighbors[j].size());
            sum += (formalValenceOfType(typeNames[j]) - tableCharges[j]) / coordination;
        }
        atoms->atom[i].q = atoms->atom[i].qB = roundCharge(baseCharge + sum);
        numSmeared++;
    }

    double totalCharge = 0;
    int    numOxygens  = 0;
    for (int i = 0; i < numAtoms; i++)
    {
        totalCharge += atoms->atom[i].q;
        if (typeNameStartsWith(typeNames[i], "o"))
        {
            numOxygens++;
        }
    }
    const double integerCharge = std::round(totalCharge);
    if (std::abs(integerCharge - totalCharge) > 1e-5 && numOxygens > 0)
    {
        const double shift = (integerCharge - totalCharge) / numOxygens;
        for (int i = 0; i < numAtoms; i++)
        {
            if (typeNameStartsWith(typeNames[i], "o"))
            {
                atoms->atom[i].q  = roundCharge(atoms->atom[i].q + shift);
                atoms->atom[i].qB = atoms->atom[i].q;
            }
        }
    }
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "Smeared the charge over %d oxygen/fluorine atoms, the total charge before "
                    "balancing was %g",
                    numSmeared,
                    totalCharge);
}

static void print_rtp(const char*                                                      filenm,
                      const char*                                                      title,
                      t_atoms*                                                         atoms,
                      const EnumerationArray<InteractionFunction, InteractionsOfType>& plist,
                      PreprocessingAtomTypes*                                          atypes,
                      int                                                              cgnr[])
{
    FILE* fp;
    int   i, tp;

    fp = gmx_fio_fopen(filenm, "w");
    fprintf(fp, "; %s\n", title);
    fprintf(fp, "\n");
    fprintf(fp, "[ %s ]\n", *atoms->resinfo[0].name);
    fprintf(fp, "\n");
    fprintf(fp, "[ atoms ]\n");
    for (i = 0; (i < atoms->nr); i++)
    {
        tp        = atoms->atom[i].type;
        auto tpnm = atypes->atomNameFromAtomType(tp);
        if (!tpnm.has_value())
        {
            gmx_fatal(FARGS, "tp = %d, i = %d in print_rtp", tp, i);
        }
        fprintf(fp, "%-8s  %12s  %8.4f  %5d\n", *atoms->atomname[i], tpnm->c_str(), atoms->atom[i].q, cgnr[i]);
    }
    print_pl(fp, plist, InteractionFunction::Bonds, "bonds", atoms->atomname);
    print_pl(fp, plist, InteractionFunction::Angles, "angles", atoms->atomname);
    print_pl(fp, plist, InteractionFunction::ProperDihedrals, "dihedrals", atoms->atomname);
    print_pl(fp, plist, InteractionFunction::ImproperDihedrals, "impropers", atoms->atomname);

    gmx_fio_fclose(fp);
}

int gmx_x2top(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] generates a primitive topology from a coordinate file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.",
        "The program can also make an [REF].rtp[ref] entry, which you can then add",
        "to the [REF].rtp[ref] database.[PAR]",
        "When [TT]-param[tt] is set, equilibrium distances and angles",
        "and force constants will be printed in the topology for all",
        "interactions. The equilibrium distances and angles are taken",
        "from the input coordinates, the force constant are set with",
        "command line options.",
        "The force fields somewhat supported currently are:[PAR]",
        "G53a5  GROMOS96 53a5 Forcefield (official distribution)[PAR]",
        "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
        "min    MINFF force field for minerals (alias: minff)[PAR]",
        "The corresponding data files can be found in the library directory",
        "with name [TT]atomname2type.n2t[tt]. Check Chapter 5 of the manual for more",
        "information about file formats. By default, the force field selection",
        "is interactive, but you can use the [TT]-ff[tt] option to specify",
        "one of the short names above on the command line instead. In that",
        "case [THISMODULE] just looks for the corresponding file.[PAR]",
        "The MINFF force field is described in Holmboe, Lavauzelle and Swai, MINFF: A Family of",
        "Classical Force Fields for Minerals with Tunable Angle Terms, arXiv preprint (2026),",
        "doi:10.48550/arXiv.2610.18220, and its parameters are available at",
        "https://github.com/mholmboe/minff.[PAR]",
        "With [TT]-ff min[tt] (or [TT]-ff minff[tt]) the topology is written in the compact",
        "MINFF format. Only bonds involving hydrogen are listed in [TT][ bonds ][tt],",
        "without parameters unless [TT]-kb[tt] is given, and the [TT][ angles ][tt]",
        "get the measured angle and the force constant, annotated with the atom types.",
        "No pairs or dihedrals are written. For these force fields the defaults are",
        "441050 for [TT]-kb[tt], 500 for [TT]-kt[tt], 110 for [TT]-ktH[tt], 1 for",
        "[TT]-nexcl[tt] and MIN for [TT]-name[tt]. The bonded parameters are measured",
        "from the structure, with periodic boundary conditions taken into account when",
        "[TT]-pbc[tt] is set. With [TT]-dM[tt] and [TT]-dH[tt] the printed bond length",
        "can be set for all bonds or for bonds to hydrogen only, and with [TT]-aH[tt] the",
        "angles involving hydrogen are printed with a fixed value. Any of these, as",
        "well as [TT]-kb[tt], also adds the bond length and force constant columns to",
        "[TT][ bonds ][tt].[PAR]",
        "Hydrogen atoms (for [TT]-ktH[tt]) are recognized from the atomic mass given in the",
        "[TT]atomname2type.n2t[tt] file (below 2.5 u), otherwise from the element column of a",
        "[REF].pdb[ref] file, and as a last resort from the atom name starting with an H.[PAR]",
        "With [TT]-smear[tt] ([TT]-ff min[tt] only) the charges of oxygens (and structural",
        "fluorines, type Fs) are not read from the table but computed from their",
        "neighbors: q = -2 + sum((Z - q_j)/CN_j) over the bonded neighbors, with Z the formal",
        "valence of the neighbor, q_j its charge from the table and CN_j its number of",
        "bonds. This smears out the charge of isomorphic substitutions, see Tournassat et al.,",
        "Clays and Clay Minerals 64, 374-388 (2016), doi:10.1346/CCMN.2016.0640403 and",
        "Lammers et al., J. Colloid Interface Sci. 2017, doi:10.1016/j.jcis.2016.11.084.",
        "The formal valence is derived from the MINFF atom type name (Si 4, Al 3, Fe2/Fe3,",
        "Ti 4, Li 1, Mg 2, Ca 2, H 1). A non-integer total charge is then rounded to the",
        "nearest integer by shifting all oxygens equally.[PAR]",
    };
    const char* bugs[] = {
        "The atom type selection is primitive. Virtually no chemical knowledge is used",
        "Periodic boundary conditions screw up the bonding",
        "No improper dihedrals are generated",
        ("The atoms to atomtype translation table is incomplete ([TT]atomname2type.n2t[tt] file in "
         "the data directory). Please extend it and send the results back to the GROMACS crew.")
    };
    FILE*                                                     fp;
    EnumerationArray<InteractionFunction, InteractionsOfType> plist;
    t_excls*                                                  excls;
    t_nm2type*                                                nm2t;
    t_mols                                                    mymol;
    int                                                       nnm;
    char                                                      forcefield[32];
    rvec*                                                     x; /* coordinates? */
    int *                                                     nbonds, *cgnr;
    int               bts[static_cast<int>(BondedTypes::Count)] = { 1, 1, 1, 2, 1, 1 };
    matrix            box;    /* box length matrix */
    int               natoms; /* number of atoms in one molecule  */
    PbcType           pbcType;
    bool              bRTP, bTOP, bOPLS, bMinff;
    real              qtot, mtot;
    gmx_output_env_t* oenv;

    t_filenm fnm[] = { { efSTX, "-f", "conf", ffREAD },
                       { efTOP, "-o", "out", ffOPTWR },
                       { efRTP, "-r", "out", ffOPTWR } };
#define NFILE asize(fnm)
    real              kb = 4e5, kt = 400, ktHydrogen = -1, kp = 5;
    real              bondLengthOverride = -1, hydrogenBondDistance = 0, hydrogenAngleDegrees = 0;
    PreprocessResidue rtp_header_settings;
    bool              bRemoveDihedralIfWithImproper = FALSE;
    bool              bGenerateHH14Interactions     = TRUE;
    bool              bKeepAllGeneratedDihedrals    = FALSE;
    int               nrexcl                        = 3;
    bool              bParam = TRUE, bRound = TRUE;
    bool              bPairs = TRUE, bPBC = TRUE;
    bool              bUsePDBcharge = FALSE, bVerbose = FALSE;
    const char*       molnm   = "ICE";
    const char*       ff      = "oplsaa";
    const char*       n2tFile = nullptr;
    bool              bSmear  = FALSE;
    t_pargs           pa[]    = {
        { "-ff",
                       FALSE,
                       etSTR,
                       { &ff },
                       "Force field for your simulation. Type \"select\" for interactive selection." },
        { "-v", FALSE, etBOOL, { &bVerbose }, "Generate verbose output in the top file." },
        { "-nexcl", FALSE, etINT, { &nrexcl }, "Number of exclusions" },
        { "-H14",
                       FALSE,
                       etBOOL,
                       { &bGenerateHH14Interactions },
                       "Use 3rd neighbour interactions for hydrogen atoms" },
        { "-alldih",
                       FALSE,
                       etBOOL,
                       { &bKeepAllGeneratedDihedrals },
                       "Generate all proper dihedrals" },
        { "-remdih",
                       FALSE,
                       etBOOL,
                       { &bRemoveDihedralIfWithImproper },
                       "Remove dihedrals on the same bond as an improper" },
        { "-pairs", FALSE, etBOOL, { &bPairs }, "Output 1-4 interactions (pairs) in topology file" },
        { "-name", FALSE, etSTR, { &molnm }, "Name of your molecule" },
        { "-pbc", FALSE, etBOOL, { &bPBC }, "Use periodic boundary conditions." },
        { "-pdbq",
                       FALSE,
                       etBOOL,
                       { &bUsePDBcharge },
                       "Use the B-factor supplied in a [REF].pdb[ref] file for the atomic charges" },
        { "-n2t",
                       FALSE,
                       etSTR,
                       { &n2tFile },
                       "Read the atom name to type translation from this file instead of the "
                                    "[TT]atomname2type.n2t[tt] files of the force field" },
        { "-smear",
                       FALSE,
                       etBOOL,
                       { &bSmear },
                       "-ff min only: compute the charges of oxygens and fluorines from the "
                                    "charges of their neighbors (charge smearing) instead of reading them from the "
                                    "n2t table" },
        { "-param", FALSE, etBOOL, { &bParam }, "Print parameters in the output" },
        { "-round", FALSE, etBOOL, { &bRound }, "Round off measured values" },
        { "-kb", FALSE, etREAL, { &kb }, "Bonded force constant (kJ/mol/nm^2)" },
        { "-kt", FALSE, etREAL, { &kt }, "Angle force constant (kJ/mol/rad^2)" },
        { "-ktH",
                       FALSE,
                       etREAL,
                       { &ktHydrogen },
                       "Angle force constant for angles involving hydrogen atoms (kJ/mol/rad^2). "
                                    "Negative means use [TT]-kt[tt] for all angles" },
        { "-dM",
                       FALSE,
                       etREAL,
                       { &bondLengthOverride },
                       "MINFF format only: print this bond length (nm) for all bonds instead of the "
                                    "default" },
        { "-dH",
                       FALSE,
                       etREAL,
                       { &hydrogenBondDistance },
                       "MINFF format only: print this bond length (nm) for bonds involving hydrogen" },
        { "-aH",
                       FALSE,
                       etREAL,
                       { &hydrogenAngleDegrees },
                       "MINFF format only: print this angle (degrees) for angles involving hydrogen "
                                    "instead of the measured one" },
        { "-kp", FALSE, etREAL, { &kp }, "Dihedral angle force constant (kJ/mol/rad^2)" }
    };

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv, nullptr))
    {
        return 0;
    }
    bRTP = opt2bSet("-r", NFILE, fnm);
    bTOP = opt2bSet("-o", NFILE, fnm);
    /* C89 requirements mean that these struct members cannot be used in
     * the declaration of pa. So some temporary variables are needed. */
    rtp_header_settings.bRemoveDihedralIfWithImproper = bRemoveDihedralIfWithImproper;
    rtp_header_settings.bGenerateHH14Interactions     = bGenerateHH14Interactions;
    rtp_header_settings.bKeepAllGeneratedDihedrals    = bKeepAllGeneratedDihedrals;
    rtp_header_settings.nrexcl                        = nrexcl;

    if (!bRTP && !bTOP)
    {
        gmx_fatal(FARGS, "Specify at least one output file");
    }

    LoggerBuilder builder;
    builder.addTargetStream(MDLogger::LogLevel::Info, &TextOutputFile::standardOutput());
    builder.addTargetStream(MDLogger::LogLevel::Warning, &TextOutputFile::standardError());
    LoggerOwner logOwner(builder.build());
    MDLogger    logger(logOwner.logger());


    /* Force field selection, interactive or direct */
    if (std::strcmp(ff, "minff") == 0)
    {
        ff = "min";
    }
    auto ffdir = choose_ff(
            std::strcmp(ff, "select") == 0 ? nullptr : ff, forcefield, sizeof(forcefield), logger);

    bOPLS  = (std::strcmp(forcefield, "oplsaa") == 0);
    bMinff = (std::strcmp(forcefield, "min") == 0);

    if (bSmear && !bMinff)
    {
        gmx_fatal(FARGS, "-smear can only be used with -ff min or minff.");
    }
    if (bSmear && bUsePDBcharge)
    {
        gmx_fatal(FARGS, "-smear cannot be combined with -pdbq.");
    }

    const bool bKbSet = opt2parg_bSet("-kb", asize(pa), pa);

    std::optional<real> hydrogenBondDistanceOpt;
    std::optional<real> bondLengthOverrideOpt;
    std::optional<real> hydrogenAngleOpt;
    if (opt2parg_bSet("-dH", asize(pa), pa))
    {
        if (hydrogenBondDistance <= 0)
        {
            gmx_fatal(FARGS, "-dH requires a positive bond length in nm.");
        }
        hydrogenBondDistanceOpt = hydrogenBondDistance;
    }
    if (opt2parg_bSet("-dM", asize(pa), pa) && bondLengthOverride > 0)
    {
        bondLengthOverrideOpt = bondLengthOverride;
    }
    if (opt2parg_bSet("-aH", asize(pa), pa))
    {
        if (hydrogenAngleDegrees <= 0)
        {
            gmx_fatal(FARGS, "-aH requires a positive angle in degrees.");
        }
        hydrogenAngleOpt = hydrogenAngleDegrees;
    }
    if (bMinff)
    {
        if (!opt2parg_bSet("-nexcl", asize(pa), pa))
        {
            nrexcl                     = 1;
            rtp_header_settings.nrexcl = nrexcl;
        }
        if (!bKbSet)
        {
            kb = 441050;
        }
        if (!opt2parg_bSet("-kt", asize(pa), pa))
        {
            kt = 500;
        }
        if (!opt2parg_bSet("-ktH", asize(pa), pa))
        {
            ktHydrogen = 110;
        }
        if (!opt2parg_bSet("-name", asize(pa), pa))
        {
            molnm = "MIN";
        }
    }


    mymol.name = molnm;
    mymol.nr   = 1;

    /* Read coordinates */
    t_topology* top;
    snew(top, 1);
    read_tps_conf(opt2fn("-f", NFILE, fnm), top, &pbcType, &x, nullptr, box, FALSE);
    t_atoms* atoms = &top->atoms;
    natoms         = atoms->nr;
    if (atoms->pdbinfo == nullptr)
    {
        snew(atoms->pdbinfo, natoms);
    }

    const ExtendableForceField forceField(ffdir);
    if (n2tFile != nullptr)
    {
        nm2t = rd_nm2type_from_file(n2tFile, &nnm);
        if (nnm == 0)
        {
            gmx_fatal(FARGS, "No atomname2type entries could be read from %s", n2tFile);
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("There are %d name to type translations in file %s", nnm, n2tFile);
    }
    else
    {
        nm2t = rd_nm2type(forceField, &nnm);
        if (nnm == 0)
        {
            gmx_fatal(FARGS,
                      "No or incorrect atomname2type.n2t file found (looking for %s)",
                      ffdir.string().c_str());
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "There are %d name to type translations in file %s", nnm, ffdir.string().c_str());
    }
    if (debug)
    {
        dump_nm2type(debug, nnm, nm2t);
    }
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Generating bonds from distances...");
    snew(nbonds, atoms->nr);
    mk_bonds(nnm, nm2t, atoms, x, &(plist[InteractionFunction::Bonds]), nbonds, bPBC, pbcType, box);

    PreprocessingAtomTypes atypes;
    set_atom_type(&atypes, atoms, &(plist[InteractionFunction::Bonds]), nbonds, nnm, nm2t, logger);
    if (bSmear)
    {
        smearCharges(atoms, plist[InteractionFunction::Bonds], &atypes, logger);
    }

    /* Make Angles and Dihedrals */
    snew(excls, atoms->nr);
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("Generating angles and dihedrals from bonds...");
    gen_pad(atoms, arrayRefFromArray(&rtp_header_settings, 1), plist, excls, {}, TRUE, {}, {});

    if (!bPairs || bMinff)
    {
        plist[InteractionFunction::LennardJones14].interactionTypes.clear();
    }
    if (bMinff)
    {
        plist[InteractionFunction::ProperDihedrals].interactionTypes.clear();
        plist[InteractionFunction::ImproperDihedrals].interactionTypes.clear();
    }
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "There are %4zu %s dihedrals, %4zu impropers, %4zu angles\n"
                    "          %4zu pairs,     %4zu bonds and  %4d atoms\n",
                    plist[InteractionFunction::ProperDihedrals].size(),
                    bOPLS ? "Ryckaert-Bellemans" : "proper",
                    plist[InteractionFunction::ImproperDihedrals].size(),
                    plist[InteractionFunction::Angles].size(),
                    plist[InteractionFunction::LennardJones14].size(),
                    plist[InteractionFunction::Bonds].size(),
                    atoms->nr);

    calc_angles_dihs(
            &plist[InteractionFunction::Angles], &plist[InteractionFunction::ProperDihedrals], x, bPBC, pbcType, box);

    set_force_const(plist, atoms, kb, kt, ktHydrogen, kp, bRound, bParam);

    cgnr = set_cgnr(atoms, bUsePDBcharge, &qtot, &mtot);
    if (bMinff)
    {
        // MINFF molecules are a single charge group
        std::fill(cgnr, cgnr + atoms->nr, 1);
    }
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Total charge is %g, total mass is %g", qtot, mtot);
    if (bOPLS)
    {
        bts[2] = 3;
        bts[3] = 1;
    }

    if (bTOP)
    {
        fp = ftp2FILE(efTOP, NFILE, fnm, "w");
        print_top_header(fp, ftp2fn(efTOP, NFILE, fnm), TRUE, ffdir, 1.0);

        if (bMinff)
        {
            const bool bExplicitBondConstant = bKbSet || hydrogenBondDistanceOpt.has_value()
                                               || bondLengthOverrideOpt.has_value();
            writeMinffTopology(fp,
                               mymol.name.c_str(),
                               rtp_header_settings.nrexcl,
                               atoms,
                               &atypes,
                               plist,
                               excls,
                               x,
                               bPBC,
                               pbcType,
                               box,
                               kb,
                               bExplicitBondConstant,
                               hydrogenBondDistanceOpt,
                               bondLengthOverrideOpt,
                               hydrogenAngleOpt);
        }
        else
        {
            write_top(fp,
                      {},
                      mymol.name.c_str(),
                      atoms,
                      FALSE,
                      bts,
                      plist,
                      excls,
                      &atypes,
                      rtp_header_settings.nrexcl);
        }
        print_top_mols(fp, mymol.name.c_str(), forceField, nullptr, {}, arrayRefFromArray(&mymol, 1));

        gmx_ffclose(fp);
    }
    if (bRTP)
    {
        print_rtp(ftp2fn(efRTP, NFILE, fnm), "Generated by x2top", atoms, plist, &atypes, cgnr);
    }

    if (debug)
    {
        dump_hybridization(debug, atoms, nbonds);
    }

    GMX_LOG(logger.warning)
            .asParagraph()
            .appendTextFormatted(
                    "Topologies generated by %s can not be trusted at face value. "
                    "Please verify atomtypes and charges by comparison to other topologies.",
                    output_env_get_program_display_name(oenv));

    for (int i = 0; i < atoms->nr; i++)
    {
        sfree(excls[i].e);
    }
    sfree(excls);
    sfree(cgnr);
    sfree(nbonds);
    free_nm2type(nnm, nm2t);
    sfree(x);
    done_top(top);
    sfree(top);
    output_env_done(oenv);

    return 0;
}

} // namespace gmx
