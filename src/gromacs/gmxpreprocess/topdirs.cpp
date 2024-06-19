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

#include "topdirs.h"

#include <cstdarg>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringtoenumvalueconverter.h"

const char* enumValueToString(Directive d)
{
    /* Must correspond to the Directive enum in topdirs.h */
    static constexpr gmx::EnumerationArray<Directive, const char*> directiveNames = {
        "defaults",
        "atomtypes",
        "bondtypes",
        "constrainttypes",
        "pairtypes",
        "angletypes",
        "dihedraltypes",
        "nonbond_params",
        "implicit_genborn_params",
        "implicit_surface_params",
        "cmaptypes",
        /* All the directives above can not appear after moleculetype */
        "moleculetype",
        "atoms",
        "virtual_sites1",
        "virtual_sites2",
        "virtual_sites3",
        "virtual_sites4",
        "virtual_sitesn",
        "bonds",
        "exclusions",
        "pairs",
        "pairs_nb",
        "angles",
        "dihedrals",
        "constraints",
        "settles",
        "polarization",
        "water_polarization",
        "thole_polarization",
        "system",
        "molecules",
        "position_restraints",
        "angle_restraints",
        "angle_restraints_z",
        "distance_restraints",
        "orientation_restraints",
        "dihedral_restraints",
        "cmap",
        "intermolecular_interactions",
        "maxdirs",
        "invalid",
        "none"
    };
    return directiveNames[d];
}

int ifunc_index(Directive d, int type)
{
    switch (d)
    {
        case Directive::d_bondtypes:
        case Directive::d_bonds:
            switch (type)
            {
                case 1: return F_BONDS;
                case 2: return F_G96BONDS;
                case 3: return F_MORSE;
                case 4: return F_CUBICBONDS;
                case 5: return F_CONNBONDS;
                case 6: return F_HARMONIC;
                case 7: return F_FENEBONDS;
                case 8: return F_TABBONDS;
                case 9: return F_TABBONDSNC;
                case 10: return F_RESTRBONDS;
                default: gmx_fatal(FARGS, "Invalid bond type %d", type);
            }
        case Directive::d_angles:
        case Directive::d_angletypes:
            switch (type)
            {
                case 1: return F_ANGLES;
                case 2: return F_G96ANGLES;
                case 3: return F_CROSS_BOND_BONDS;
                case 4: return F_CROSS_BOND_ANGLES;
                case 5: return F_UREY_BRADLEY;
                case 6: return F_QUARTIC_ANGLES;
                case 8: return F_TABANGLES;
                case 9: return F_LINEAR_ANGLES;
                case 10: return F_RESTRANGLES;
                default: gmx_fatal(FARGS, "Invalid angle type %d", type);
            }
        case Directive::d_pairs:
        case Directive::d_pairtypes:
            if (type == 1 || (d == Directive::d_pairtypes && type == 2))
            {
                return F_LJ14;
            }
            else if (type == 2)
            {
                return F_LJC14_Q;
            }
            else
            {
                gmx_fatal(FARGS, "Invalid pairs type %d", type);
            }
        case Directive::d_pairs_nb: return F_LJC_PAIRS_NB;
        case Directive::d_dihedrals:
        case Directive::d_dihedraltypes:
            switch (type)
            {
                case 1: return F_PDIHS;
                case 2: return F_IDIHS;
                case 3: return F_RBDIHS;
                case 4: return F_PIDIHS;
                case 5: return F_FOURDIHS;
                case 8: return F_TABDIHS;
                case 9:
                    return F_PDIHS; /* proper dihedrals where we allow multiple terms over single bond */
                case 10: return F_RESTRDIHS;
                case 11: return F_CBTDIHS;
                default: gmx_fatal(FARGS, "Invalid dihedral type %d", type);
            }
        case Directive::d_cmaptypes:
        case Directive::d_cmap: return F_CMAP;

        case Directive::d_nonbond_params:
            if (type == 1)
            {
                return F_LJ;
            }
            else
            {
                return F_BHAM;
            }
        case Directive::d_vsites1:
            if (type == 1)
            {
                return F_VSITE1;
            }
            else
            {
                gmx_fatal(FARGS, "Invalid vsites1 type %d", type);
            }
        case Directive::d_vsites2:
            switch (type)
            {
                case 1: return F_VSITE2;
                case 2: return F_VSITE2FD;
                default: gmx_fatal(FARGS, "Invalid vsites2 type %d", type);
            }
        case Directive::d_vsites3:
            switch (type)
            {
                case 1: return F_VSITE3;
                case 2: return F_VSITE3FD;
                case 3: return F_VSITE3FAD;
                case 4: return F_VSITE3OUT;
                default: gmx_fatal(FARGS, "Invalid vsites3 type %d", type);
            }
        case Directive::d_vsites4:
            switch (type)
            {
                case 1: return F_VSITE4FD;
                case 2: return F_VSITE4FDN;
                default: gmx_fatal(FARGS, "Invalid vsites4 type %d", type);
            }
        case Directive::d_vsitesn: return F_VSITEN;
        case Directive::d_constraints:
        case Directive::d_constrainttypes:
            switch (type)
            {
                case 1: return F_CONSTR;
                case 2: return F_CONSTRNC;
                default: gmx_fatal(FARGS, "Invalid constraints type %d", type);
            }
        case Directive::d_settles: return F_SETTLE;
        case Directive::d_position_restraints:
            switch (type)
            {
                case 1: return F_POSRES;
                case 2: return F_FBPOSRES;
                default: gmx_fatal(FARGS, "Invalid position restraint type %d", type);
            }
        case Directive::d_polarization:
            switch (type)
            {
                case 1: return F_POLARIZATION;
                case 2: return F_ANHARM_POL;
                default: gmx_fatal(FARGS, "Invalid polarization type %d", type);
            }
        case Directive::d_thole_polarization: return F_THOLE_POL;
        case Directive::d_water_polarization: return F_WATER_POL;
        case Directive::d_angle_restraints: return F_ANGRES;
        case Directive::d_angle_restraints_z: return F_ANGRESZ;
        case Directive::d_distance_restraints: return F_DISRES;
        case Directive::d_orientation_restraints: return F_ORIRES;
        case Directive::d_dihedral_restraints: return F_DIHRES;
        default:
            gmx_fatal(FARGS, "invalid directive %s in ifunc_index (%s:%d)", enumValueToString(d), __FILE__, __LINE__);
    }
}

enum class DeprecatedDirectives : int
{
    d_dummies1,
    d_dummies2,
    d_dummies3,
    d_dummies4,
    d_dummiesn,
    Count
};

static const char* enumValueToString(DeprecatedDirectives d)
{
    static constexpr gmx::EnumerationArray<DeprecatedDirectives, const char*> directiveNames = {
        "dummies1", "dummies2", "dummies3", "dummies4", "dummiesn"
    };
    return directiveNames[d];
}

Directive str2dir(char* dstr)
{
    static const gmx::StringToEnumValueConverter<Directive, enumValueToString, gmx::StringCompareType::CaseAndDashInsensitive> s_converter;

    if (std::optional<Directive> d = s_converter.valueFrom(dstr); d.has_value())
    {
        return d.value();
    }
    // Also handle deprecated directives that have modern replacements, like
    // "dummies*" -> "virtual_sites*"

    static const gmx::StringToEnumValueConverter<DeprecatedDirectives, enumValueToString, gmx::StringCompareType::CaseAndDashInsensitive>
            s_converterForDeprecated;

    if (std::optional<DeprecatedDirectives> d = s_converterForDeprecated.valueFrom(dstr); d.has_value())
    {
        static constexpr gmx::EnumerationArray<DeprecatedDirectives, Directive> s_deprecatedDirectiveToDirective = {
            Directive::d_vsites1, Directive::d_vsites2, Directive::d_vsites3,
            Directive::d_vsites4, Directive::d_vsitesn,
        };
        return s_deprecatedDirectiveToDirective[d.value()];
    }
    return Directive::d_invalid;
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx::EnumerationArray<Directive, Directive*> necessary = { { nullptr } };

static void set_nec(Directive** n, ...)
/* Must always have at least one extra argument */
{
    std::va_list ap;
    int          ind = 0;
    Directive    d;

    va_start(ap, n);
    do
    {
        d = static_cast<Directive>(va_arg(ap, int));
        srenew(*n, ++ind);
        (*n)[ind - 1] = d;
    } while (d != Directive::d_none);
    va_end(ap);
}

void DS_Init(DirStack** DS)
{
    if (necessary[0] == nullptr)
    {
        set_nec(&(necessary[Directive::d_defaults]), Directive::d_none);
        set_nec(&(necessary[Directive::d_atomtypes]), Directive::d_defaults, Directive::d_none);
        set_nec(&(necessary[Directive::d_bondtypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_constrainttypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_pairtypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_angletypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_dihedraltypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_nonbond_params]), Directive::d_atomtypes, Directive::d_none);
        // Note that the content of the next two directives are
        // ignored, but if grompp reads them in old force field files,
        // it still needs to understand that they are in a valid place
        // in the .top structure. It doesn't have to require them to
        // be in the same place that was valid in old versions (ie. child
        // directive of [atomtypes]) but any relevant case will
        // satisfy that.
        set_nec(&(necessary[Directive::d_implicit_genborn_params]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_implicit_surface_params]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_cmaptypes]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_moleculetype]), Directive::d_atomtypes, Directive::d_none);
        set_nec(&(necessary[Directive::d_atoms]), Directive::d_moleculetype, Directive::d_none);
        set_nec(&(necessary[Directive::d_vsites1]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_vsites2]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_vsites3]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_vsites4]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_vsitesn]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_bonds]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_exclusions]),
                Directive::d_bonds,
                Directive::d_constraints,
                Directive::d_settles,
                Directive::d_none);
        set_nec(&(necessary[Directive::d_pairs]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_pairs_nb]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_angles]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_polarization]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_water_polarization]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_thole_polarization]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_dihedrals]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_constraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_settles]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_system]), Directive::d_moleculetype, Directive::d_none);
        set_nec(&(necessary[Directive::d_molecules]), Directive::d_system, Directive::d_none);
        set_nec(&(necessary[Directive::d_position_restraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_angle_restraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_angle_restraints_z]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_distance_restraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_orientation_restraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_dihedral_restraints]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_cmap]), Directive::d_atoms, Directive::d_none);
        set_nec(&(necessary[Directive::d_intermolecular_interactions]),
                Directive::d_molecules,
                Directive::d_none);
    }
    *DS = nullptr;
}

void DS_Done(DirStack** DS)
{
    DirStack* D;

    while (*DS != nullptr)
    {
        D   = *DS;
        *DS = (*DS)->prev;
        sfree(D);
    }
}

void DS_Push(DirStack** DS, Directive d)
{
    DirStack* D;

    snew(D, 1);
    D->d    = d;
    D->prev = *DS;
    *DS     = D;
}

int DS_Search(DirStack* DS, Directive d)
{
    DirStack* D;

    D = DS;
    while ((D != nullptr) && (D->d != d))
    {
        D = D->prev;
    }

    return static_cast<int>(D != nullptr);
}

int DS_Check_Order(DirStack* DS, Directive d)
{
    Directive d0;
    int       i = 0;

    /* Check if parameter definitions appear after a moleculetype directive */
    if (d < Directive::d_moleculetype && DS_Search(DS, Directive::d_moleculetype))
    {
        return FALSE;
    }

    /* Check if all the necessary directives have appeared before directive d */
    if (necessary[d][0] == Directive::d_none)
    {
        return TRUE;
    }
    else
    {
        do
        {
            d0 = necessary[d][i++];
            if (DS_Search(DS, d0))
            {
                return TRUE;
            }
        } while (d0 != Directive::d_none);
    }
    return FALSE;
}
