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

#include "vsite_parm.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"
#include "resall.h"

/*! \internal \brief
 * Data type to store information about bonded interactions for virtual sites.
 */
class VsiteBondedInteraction
{
public:
    //! Constructor initializes datastructure.
    VsiteBondedInteraction(gmx::ArrayRef<const int> atomIndex, real parameterValue) :
        parameterValue_(parameterValue)
    {
        GMX_RELEASE_ASSERT(atomIndex.size() <= atomIndex_.size(),
                           "Cannot add more atom indices than maximum number");
        std::array<int, 4>::iterator atomIndexIt = atomIndex_.begin();
        for (const auto index : atomIndex)
        {
            *atomIndexIt++ = index;
        }
    }
    /*!@{*/
    //! Access the individual elements set for the vsite bonded interaction.
    const int& ai() const { return atomIndex_[0]; }
    const int& aj() const { return atomIndex_[1]; }
    const int& ak() const { return atomIndex_[2]; }
    const int& al() const { return atomIndex_[3]; }

    const real& parameterValue() const { return parameterValue_; }
    /*!@}*/

private:
    //! The distance value for this bonded interaction.
    real parameterValue_;
    //! Array of atom indices
    std::array<int, 4> atomIndex_;
};

/*! \internal \brief
 * Stores information about single virtual site bonded parameter.
 */
struct VsiteBondParameter
{
    //! Constructor initializes datastructure.
    VsiteBondParameter(int ftype, const InteractionOfType& vsiteInteraction) :
        ftype_(ftype), vsiteInteraction_(vsiteInteraction)
    {
    }
    //! Function type for virtual site.
    int ftype_;
    /*! \brief
     * Interaction type data.
     *
     * The datastructure should never be used in a case where the InteractionType
     * used to construct it might go out of scope before this object, as it would cause
     * the reference to type_ to dangle.
     */
    const InteractionOfType& vsiteInteraction_;
};

/*! \internal \brief
 * Helper type for conversion of bonded parameters to virtual sites.
 */
struct Atom2VsiteBond
{
    //! The vsite parameters in a list.
    std::vector<VsiteBondParameter> vSiteBondedParameters;
};

//! Convenience type def for virtual site connections.
using Atom2VsiteConnection = std::vector<int>;

static int vsite_bond_nrcheck(int ftype)
{
    int nrcheck;

    if ((interaction_function[ftype].flags & (IF_BTYPE | IF_CONSTRAINT | IF_ATYPE)) || (ftype == F_IDIHS))
    {
        nrcheck = NRAL(ftype);
    }
    else
    {
        nrcheck = 0;
    }

    return nrcheck;
}

static void enter_bonded(int nratoms, std::vector<VsiteBondedInteraction>* bondeds, const InteractionOfType& type)
{
    GMX_RELEASE_ASSERT(nratoms == type.atoms().ssize(), "Size of atom array must match");
    bondeds->emplace_back(type.atoms(), type.c0());
}

/*! \internal \brief
 * Wraps the datastructures for the different vsite bondeds.
 */
struct AllVsiteBondedInteractions
{
    //! Bond vsites.
    std::vector<VsiteBondedInteraction> bonds;
    //! Angle vsites.
    std::vector<VsiteBondedInteraction> angles;
    //! Dihedral vsites.
    std::vector<VsiteBondedInteraction> dihedrals;
};

static AllVsiteBondedInteractions createVsiteBondedInformation(int                      nrat,
                                                               gmx::ArrayRef<const int> atoms,
                                                               gmx::ArrayRef<const Atom2VsiteBond> at2vb)
{
    AllVsiteBondedInteractions allVsiteBondeds;
    for (int k = 0; k < nrat; k++)
    {
        for (const auto& vsite : at2vb[atoms[k]].vSiteBondedParameters)
        {
            int                      ftype   = vsite.ftype_;
            const InteractionOfType& type    = vsite.vsiteInteraction_;
            int                      nrcheck = vsite_bond_nrcheck(ftype);
            /* abuse nrcheck to see if we're adding bond, angle or idih */
            switch (nrcheck)
            {
                case 2: enter_bonded(nrcheck, &allVsiteBondeds.bonds, type); break;
                case 3: enter_bonded(nrcheck, &allVsiteBondeds.angles, type); break;
                case 4: enter_bonded(nrcheck, &allVsiteBondeds.dihedrals, type); break;
            }
        }
    }
    return allVsiteBondeds;
}

static std::vector<Atom2VsiteBond> make_at2vsitebond(int natoms, gmx::ArrayRef<InteractionsOfType> plist)
{
    bool* bVSI;

    std::vector<Atom2VsiteBond> at2vb(natoms);

    snew(bVSI, natoms);
    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) && ftype != F_VSITEN)
        {
            for (int i = 0; (i < gmx::ssize(plist[ftype])); i++)
            {
                gmx::ArrayRef<const int> atoms = plist[ftype].interactionTypes[i].atoms();
                for (int j = 0; j < NRAL(ftype); j++)
                {
                    bVSI[atoms[j]] = TRUE;
                }
            }
        }
    }

    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        int nrcheck = vsite_bond_nrcheck(ftype);
        if (nrcheck > 0)
        {
            for (int i = 0; (i < gmx::ssize(plist[ftype])); i++)
            {
                gmx::ArrayRef<const int> aa = plist[ftype].interactionTypes[i].atoms();
                for (int j = 0; j < nrcheck; j++)
                {
                    if (bVSI[aa[j]])
                    {
                        at2vb[aa[j]].vSiteBondedParameters.emplace_back(
                                ftype, plist[ftype].interactionTypes[i]);
                    }
                }
            }
        }
    }
    sfree(bVSI);

    return at2vb;
}

static std::vector<Atom2VsiteConnection> make_at2vsitecon(int natoms, gmx::ArrayRef<InteractionsOfType> plist)
{
    std::vector<bool>                 bVSI(natoms);
    std::vector<Atom2VsiteConnection> at2vc(natoms);

    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_VSITE) && ftype != F_VSITEN)
        {
            for (int i = 0; (i < gmx::ssize(plist[ftype])); i++)
            {
                gmx::ArrayRef<const int> atoms = plist[ftype].interactionTypes[i].atoms();
                for (int j = 0; j < NRAL(ftype); j++)
                {
                    bVSI[atoms[j]] = TRUE;
                }
            }
        }
    }

    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_CONSTRAINT)
        {
            for (int i = 0; (i < gmx::ssize(plist[ftype])); i++)
            {
                int ai = plist[ftype].interactionTypes[i].ai();
                int aj = plist[ftype].interactionTypes[i].aj();
                if (bVSI[ai] && bVSI[aj])
                {
                    /* Store forward direction */
                    at2vc[ai].emplace_back(aj);
                    /* Store backward direction */
                    at2vc[aj].emplace_back(ai);
                }
            }
        }
    }
    return at2vc;
}

/* for debug */
static void print_bad(FILE*                                       fp,
                      gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                      gmx::ArrayRef<const VsiteBondedInteraction> angles,
                      gmx::ArrayRef<const VsiteBondedInteraction> idihs)
{
    if (!bonds.empty())
    {
        fprintf(fp, "bonds:");
        for (const auto& bond : bonds)
        {
            fprintf(fp, " %d-%d (%g)", bond.ai() + 1, bond.aj() + 1, bond.parameterValue());
        }
        fprintf(fp, "\n");
    }
    if (!angles.empty())
    {
        fprintf(fp, "angles:");
        for (const auto& angle : angles)
        {
            fprintf(fp, " %d-%d-%d (%g)", angle.ai() + 1, angle.aj() + 1, angle.ak() + 1, angle.parameterValue());
        }
        fprintf(fp, "\n");
    }
    if (!idihs.empty())
    {
        fprintf(fp, "idihs:");
        for (const auto& idih : idihs)
        {
            fprintf(fp,
                    " %d-%d-%d-%d (%g)",
                    idih.ai() + 1,
                    idih.aj() + 1,
                    idih.ak() + 1,
                    idih.al() + 1,
                    idih.parameterValue());
        }
        fprintf(fp, "\n");
    }
}

static void printInteractionOfType(FILE* fp, int ftype, int i, const InteractionOfType& type)
{
    static int pass       = 0;
    static int prev_ftype = NOTSET;
    static int prev_i     = NOTSET;

    if ((ftype != prev_ftype) || (i != prev_i))
    {
        pass       = 0;
        prev_ftype = ftype;
        prev_i     = i;
    }
    fprintf(fp, "(%d) plist[%s].param[%d]", pass, interaction_function[ftype].name, i);
    gmx::ArrayRef<const real> forceParam = type.forceParam();
    for (int j = 0; j < NRFP(ftype); j++)
    {
        fprintf(fp, ".c[%d]=%g ", j, forceParam[j]);
    }
    fprintf(fp, "\n");
    pass++;
}

static real get_bond_length(gmx::ArrayRef<const VsiteBondedInteraction> bonds, t_iatom ai, t_iatom aj)
{
    real bondlen;

    bondlen = NOTSET;
    for (const auto& bond : bonds)
    {
        /* check both ways */
        if (((ai == bond.ai()) && (aj == bond.aj())) || ((ai == bond.aj()) && (aj == bond.ai())))
        {
            bondlen = bond.parameterValue(); /* note: parameterValue might be NOTSET */
            break;
        }
    }
    return bondlen;
}

static real get_angle(gmx::ArrayRef<const VsiteBondedInteraction> angles, t_iatom ai, t_iatom aj, t_iatom ak)
{
    real angle;

    angle = NOTSET;
    for (const auto& ang : angles)
    {
        /* check both ways */
        if (((ai == ang.ai()) && (aj == ang.aj()) && (ak == ang.ak()))
            || ((ai == ang.ak()) && (aj == ang.aj()) && (ak == ang.ai())))
        {
            angle = gmx::c_deg2Rad * ang.parameterValue();
            break;
        }
    }
    return angle;
}

static bool cmp_atomtype_name_AB(t_atom*                 atom,
                                 PreprocessingAtomTypes* atypes,
                                 const std::string&      compare,
                                 int                     compareLength)
{
    /* When using the decoupling option, atom types are changed
     * to decoupled for the non-bonded interactions, but the virtual
     * sites constructions should be based on the "normal" interactions.
     * So we return the state B atom type name if the state A atom
     * type is the decoupled one. We should actually check for the atom
     * type number, but that's not passed here. So we check for
     * the decoupled atom type name. This should not cause trouble
     * as this code is only used for topologies with v-sites without
     * parameters generated by pdb2gmx.
     */
    if (atypes->atomNameFromAtomType(atom->type)->find("decoupled") != std::string::npos)
    {
        return gmx::equalCaseInsensitive(
                atypes->atomNameFromAtomType(atom->typeB).value(), compare, compareLength);
    }
    else
    {
        return gmx::equalCaseInsensitive(
                atypes->atomNameFromAtomType(atom->type).value(), compare, compareLength);
    }
}

static bool calc_vsite3_param(PreprocessingAtomTypes*                     atypes,
                              InteractionOfType*                          vsite,
                              t_atoms*                                    at,
                              gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                              gmx::ArrayRef<const VsiteBondedInteraction> angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     */

    bool bXH3, bError;
    real bjk, bjl, a = -1, b = -1;
    /* check if this is part of a NH3 , NH2-umbrella or CH3 group,
     * i.e. if atom k and l are dummy masses (MNH* or MCH3*) */
    bXH3 = ((cmp_atomtype_name_AB(&at->atom[vsite->ak()], atypes, "MNH", 3)
             && (cmp_atomtype_name_AB(&at->atom[vsite->al()], atypes, "MNH", 3)))
            || ((cmp_atomtype_name_AB(&at->atom[vsite->ak()], atypes, "MCH3", 4)
                 && (cmp_atomtype_name_AB(&at->atom[vsite->al()], atypes, "MCH3", 4)))));

    bjk    = get_bond_length(bonds, vsite->aj(), vsite->ak());
    bjl    = get_bond_length(bonds, vsite->aj(), vsite->al());
    bError = (bjk == NOTSET) || (bjl == NOTSET);
    if (bXH3)
    {
        /* now we get some XH2/XH3 group specific construction */
        /* note: we call the heavy atom 'C' and the X atom 'N' */
        real bMM, bCM, bCN, bNH, aCNH, dH, rH, dM, rM;
        int  aN;

        /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
        bError = bError || (bjk != bjl);

        /* the X atom (C or N) in the XH2/XH3 group is the first after the masses: */
        aN = std::max(vsite->ak(), vsite->al()) + 1;

        /* get common bonds */
        bMM    = get_bond_length(bonds, vsite->ak(), vsite->al());
        bCM    = bjk;
        bCN    = get_bond_length(bonds, vsite->aj(), aN);
        bError = bError || (bMM == NOTSET) || (bCN == NOTSET);

        /* calculate common things */
        rM = 0.5 * bMM;
        dM = std::sqrt(gmx::square(bCM) - gmx::square(rM));

        /* are we dealing with the X atom? */
        if (vsite->ai() == aN)
        {
            /* this is trivial */
            a = b = 0.5 * bCN / dM;
        }
        else
        {
            /* get other bondlengths and angles: */
            bNH    = get_bond_length(bonds, aN, vsite->ai());
            aCNH   = get_angle(angles, vsite->aj(), aN, vsite->ai());
            bError = bError || (bNH == NOTSET) || (aCNH == NOTSET);

            /* calculate */
            dH = bCN - bNH * std::cos(aCNH);
            rH = bNH * std::sin(aCNH);

            a = 0.5 * (dH / dM + rH / rM);
            b = 0.5 * (dH / dM - rH / rM);
        }
    }
    else
    {
        gmx_fatal(FARGS,
                  "calc_vsite3_param not implemented for the general case "
                  "(atom %d)",
                  vsite->ai() + 1);
    }
    vsite->setForceParameter(0, a);
    vsite->setForceParameter(1, b);

    return bError;
}

static bool calc_vsite3fd_param(InteractionOfType*                          vsite,
                                gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                                gmx::ArrayRef<const VsiteBondedInteraction> angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     */

    bool bError;
    real bij, bjk, bjl, aijk, aijl, rk, rl;

    bij  = get_bond_length(bonds, vsite->ai(), vsite->aj());
    bjk  = get_bond_length(bonds, vsite->aj(), vsite->ak());
    bjl  = get_bond_length(bonds, vsite->aj(), vsite->al());
    aijk = get_angle(angles, vsite->ai(), vsite->aj(), vsite->ak());
    aijl = get_angle(angles, vsite->ai(), vsite->aj(), vsite->al());
    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) || (aijk == NOTSET) || (aijl == NOTSET);

    rk = bjk * std::sin(aijk);
    rl = bjl * std::sin(aijl);
    vsite->setForceParameter(0, rk / (rk + rl));
    vsite->setForceParameter(1, -bij);

    return bError;
}

static bool calc_vsite3fad_param(InteractionOfType*                          vsite,
                                 gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                                 gmx::ArrayRef<const VsiteBondedInteraction> angles)
{
    /* i = virtual site          |
     * j = 1st bonded heavy atom | i-j
     * k = 2nd bonded heavy atom |    `k-l
     * l = 3d bonded heavy atom  |
     */

    bool bSwapParity, bError;
    real bij, aijk;

    bSwapParity = (vsite->c1() == -1);

    bij    = get_bond_length(bonds, vsite->ai(), vsite->aj());
    aijk   = get_angle(angles, vsite->ai(), vsite->aj(), vsite->ak());
    bError = (bij == NOTSET) || (aijk == NOTSET);

    vsite->setForceParameter(1, bij);
    vsite->setForceParameter(0, gmx::c_rad2Deg * aijk);

    if (bSwapParity)
    {
        vsite->setForceParameter(0, 360 - vsite->c0());
    }

    return bError;
}

static bool calc_vsite3out_param(PreprocessingAtomTypes*                     atypes,
                                 InteractionOfType*                          vsite,
                                 t_atoms*                                    at,
                                 gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                                 gmx::ArrayRef<const VsiteBondedInteraction> angles)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j
     * k,l = 2nd bonded atoms    |    `l
     * NOTE: i is out of the j-k-l plane!
     */

    bool bXH3, bError, bSwapParity;
    real bij, bjk, bjl, aijk, aijl, akjl, pijk, pijl, a, b, c;

    /* check if this is part of a NH2-umbrella, NH3 or CH3 group,
     * i.e. if atom k and l are dummy masses (MNH* or MCH3*) */
    bXH3 = ((cmp_atomtype_name_AB(&at->atom[vsite->ak()], atypes, "MNH", 3)
             && (cmp_atomtype_name_AB(&at->atom[vsite->al()], atypes, "MNH", 3)))
            || ((cmp_atomtype_name_AB(&at->atom[vsite->ak()], atypes, "MCH3", 4)
                 && (cmp_atomtype_name_AB(&at->atom[vsite->al()], atypes, "MCH3", 4)))));

    /* check if construction parity must be swapped */
    bSwapParity = (vsite->c1() == -1);

    bjk    = get_bond_length(bonds, vsite->aj(), vsite->ak());
    bjl    = get_bond_length(bonds, vsite->aj(), vsite->al());
    bError = (bjk == NOTSET) || (bjl == NOTSET);
    if (bXH3)
    {
        /* now we get some XH3 group specific construction */
        /* note: we call the heavy atom 'C' and the X atom 'N' */
        real bMM, bCM, bCN, bNH, aCNH, dH, rH, rHx, rHy, dM, rM;
        int  aN;

        /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
        bError = bError || (bjk != bjl);

        /* the X atom (C or N) in the XH3 group is the first after the masses: */
        aN = std::max(vsite->ak(), vsite->al()) + 1;

        /* get all bondlengths and angles: */
        bMM    = get_bond_length(bonds, vsite->ak(), vsite->al());
        bCM    = bjk;
        bCN    = get_bond_length(bonds, vsite->aj(), aN);
        bNH    = get_bond_length(bonds, aN, vsite->ai());
        aCNH   = get_angle(angles, vsite->aj(), aN, vsite->ai());
        bError = bError || (bMM == NOTSET) || (bCN == NOTSET) || (bNH == NOTSET) || (aCNH == NOTSET);

        /* calculate */
        dH = bCN - bNH * std::cos(aCNH);
        rH = bNH * std::sin(aCNH);
        /* we assume the H's are symmetrically distributed */
        rHx = rH * std::cos(gmx::c_deg2Rad * 30);
        rHy = rH * std::sin(gmx::c_deg2Rad * 30);
        rM  = 0.5 * bMM;
        dM  = std::sqrt(gmx::square(bCM) - gmx::square(rM));
        a   = 0.5 * ((dH / dM) - (rHy / rM));
        b   = 0.5 * ((dH / dM) + (rHy / rM));
        c   = rHx / (2 * dM * rM);
    }
    else
    {
        /* this is the general construction */

        bij  = get_bond_length(bonds, vsite->ai(), vsite->aj());
        aijk = get_angle(angles, vsite->ai(), vsite->aj(), vsite->ak());
        aijl = get_angle(angles, vsite->ai(), vsite->aj(), vsite->al());
        akjl = get_angle(angles, vsite->ak(), vsite->aj(), vsite->al());
        bError = bError || (bij == NOTSET) || (aijk == NOTSET) || (aijl == NOTSET) || (akjl == NOTSET);

        pijk = std::cos(aijk) * bij;
        pijl = std::cos(aijl) * bij;
        a = (pijk + (pijk * std::cos(akjl) - pijl) * std::cos(akjl) / gmx::square(std::sin(akjl))) / bjk;
        b = (pijl + (pijl * std::cos(akjl) - pijk) * std::cos(akjl) / gmx::square(std::sin(akjl))) / bjl;
        c = -std::sqrt(gmx::square(bij)
                       - (gmx::square(pijk) - 2 * pijk * pijl * std::cos(akjl) + gmx::square(pijl))
                                 / gmx::square(std::sin(akjl)))
            / (bjk * bjl * std::sin(akjl));
    }

    vsite->setForceParameter(0, a);
    vsite->setForceParameter(1, b);
    if (bSwapParity)
    {
        vsite->setForceParameter(2, -c);
    }
    else
    {
        vsite->setForceParameter(2, c);
    }
    return bError;
}

static bool calc_vsite4fd_param(InteractionOfType*                          vsite,
                                gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                                gmx::ArrayRef<const VsiteBondedInteraction> angles,
                                const gmx::MDLogger&                        logger)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j-m
     * k,l,m = 2nd bonded atoms  |    `l
     */

    bool bError;
    real bij, bjk, bjl, bjm, aijk, aijl, aijm, akjm, akjl;
    real pk, pl, pm, cosakl, cosakm, sinakl, sinakm, cl, cm;

    bij  = get_bond_length(bonds, vsite->ai(), vsite->aj());
    bjk  = get_bond_length(bonds, vsite->aj(), vsite->ak());
    bjl  = get_bond_length(bonds, vsite->aj(), vsite->al());
    bjm  = get_bond_length(bonds, vsite->aj(), vsite->am());
    aijk = get_angle(angles, vsite->ai(), vsite->aj(), vsite->ak());
    aijl = get_angle(angles, vsite->ai(), vsite->aj(), vsite->al());
    aijm = get_angle(angles, vsite->ai(), vsite->aj(), vsite->am());
    akjm = get_angle(angles, vsite->ak(), vsite->aj(), vsite->am());
    akjl = get_angle(angles, vsite->ak(), vsite->aj(), vsite->al());
    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) || (bjm == NOTSET) || (aijk == NOTSET)
             || (aijl == NOTSET) || (aijm == NOTSET) || (akjm == NOTSET) || (akjl == NOTSET);

    if (!bError)
    {
        pk     = bjk * std::sin(aijk);
        pl     = bjl * std::sin(aijl);
        pm     = bjm * std::sin(aijm);
        cosakl = (std::cos(akjl) - std::cos(aijk) * std::cos(aijl)) / (std::sin(aijk) * std::sin(aijl));
        cosakm = (std::cos(akjm) - std::cos(aijk) * std::cos(aijm)) / (std::sin(aijk) * std::sin(aijm));
        if (cosakl < -1 || cosakl > 1 || cosakm < -1 || cosakm > 1)
        {
            GMX_LOG(logger.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "virtual site %d: angle ijk = %f, angle ijl = %f, angle ijm = %f",
                            vsite->ai() + 1,
                            gmx::c_rad2Deg * aijk,
                            gmx::c_rad2Deg * aijl,
                            gmx::c_rad2Deg * aijm);
            gmx_fatal(FARGS,
                      "invalid construction in calc_vsite4fd for atom %d: "
                      "cosakl=%f, cosakm=%f\n",
                      vsite->ai() + 1,
                      cosakl,
                      cosakm);
        }
        sinakl = std::sqrt(1 - gmx::square(cosakl));
        sinakm = std::sqrt(1 - gmx::square(cosakm));

        /* note: there is a '+' because of the way the sines are calculated */
        cl = -pk / (pl * cosakl - pk + pl * sinakl * (pm * cosakm - pk) / (pm * sinakm));
        cm = -pk / (pm * cosakm - pk + pm * sinakm * (pl * cosakl - pk) / (pl * sinakl));

        vsite->setForceParameter(0, cl);
        vsite->setForceParameter(1, cm);
        vsite->setForceParameter(2, -bij);
    }

    return bError;
}


static bool calc_vsite4fdn_param(InteractionOfType*                          vsite,
                                 gmx::ArrayRef<const VsiteBondedInteraction> bonds,
                                 gmx::ArrayRef<const VsiteBondedInteraction> angles,
                                 const gmx::MDLogger&                        logger)
{
    /* i = virtual site          |    ,k
     * j = 1st bonded heavy atom | i-j-m
     * k,l,m = 2nd bonded atoms  |    `l
     */

    bool bError;
    real bij, bjk, bjl, bjm, aijk, aijl, aijm;
    real pk, pl, pm, a, b;

    bij  = get_bond_length(bonds, vsite->ai(), vsite->aj());
    bjk  = get_bond_length(bonds, vsite->aj(), vsite->ak());
    bjl  = get_bond_length(bonds, vsite->aj(), vsite->al());
    bjm  = get_bond_length(bonds, vsite->aj(), vsite->am());
    aijk = get_angle(angles, vsite->ai(), vsite->aj(), vsite->ak());
    aijl = get_angle(angles, vsite->ai(), vsite->aj(), vsite->al());
    aijm = get_angle(angles, vsite->ai(), vsite->aj(), vsite->am());

    bError = (bij == NOTSET) || (bjk == NOTSET) || (bjl == NOTSET) || (bjm == NOTSET)
             || (aijk == NOTSET) || (aijl == NOTSET) || (aijm == NOTSET);

    if (!bError)
    {

        /* Calculate component of bond j-k along the direction i-j */
        pk = -bjk * std::cos(aijk);

        /* Calculate component of bond j-l along the direction i-j */
        pl = -bjl * std::cos(aijl);

        /* Calculate component of bond j-m along the direction i-j */
        pm = -bjm * std::cos(aijm);

        if (std::fabs(pl) < 1000 * GMX_REAL_MIN || std::fabs(pm) < 1000 * GMX_REAL_MIN)
        {
            GMX_LOG(logger.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "virtual site %d: angle ijk = %f, angle ijl = %f, angle ijm = %f",
                            vsite->ai() + 1,
                            gmx::c_rad2Deg * aijk,
                            gmx::c_rad2Deg * aijl,
                            gmx::c_rad2Deg * aijm);
            gmx_fatal(FARGS,
                      "invalid construction in calc_vsite4fdn for atom %d: "
                      "pl=%f, pm=%f\n",
                      vsite->ai() + 1,
                      pl,
                      pm);
        }

        a = pk / pl;
        b = pk / pm;

        vsite->setForceParameter(0, a);
        vsite->setForceParameter(1, b);
        vsite->setForceParameter(2, bij);
    }

    return bError;
}


int set_vsites(bool                              bVerbose,
               t_atoms*                          atoms,
               PreprocessingAtomTypes*           atypes,
               gmx::ArrayRef<InteractionsOfType> plist,
               const gmx::MDLogger&              logger)
{
    int  ftype;
    int  nvsite;
    bool bFirst, bERROR;

    bFirst = TRUE;
    nvsite = 0;

    /* Make a reverse list to avoid ninteractions^2 operations */
    std::vector<Atom2VsiteBond> at2vb = make_at2vsitebond(atoms->nr, plist);

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nvsite += plist[ftype].size();

            if (ftype == F_VSITEN)
            {
                /* We don't calculate parameters for VSITEN */
                continue;
            }

            int i = 0;
            for (auto& param : plist[ftype].interactionTypes)
            {
                /* check if all parameters are set */
                bool                      bSet       = true;
                gmx::ArrayRef<const real> forceParam = param.forceParam();
                for (int j = 0; (j < NRFP(ftype)) && bSet; j++)
                {
                    bSet = forceParam[j] != NOTSET;
                }

                if (debug)
                {
                    fprintf(debug, "bSet=%s ", gmx::boolToString(bSet));
                    printInteractionOfType(debug, ftype, i, plist[ftype].interactionTypes[i]);
                }
                if (!bSet)
                {
                    if (bVerbose && bFirst)
                    {
                        GMX_LOG(logger.info)
                                .asParagraph()
                                .appendTextFormatted("Calculating parameters for virtual sites");
                        bFirst = FALSE;
                    }

                    /* now set the vsite parameters: */
                    AllVsiteBondedInteractions allVsiteBondeds =
                            createVsiteBondedInformation(NRAL(ftype), param.atoms(), at2vb);
                    if (debug)
                    {
                        fprintf(debug,
                                "Found %zu bonds, %zu angles and %zu idihs "
                                "for virtual site %d (%s)\n",
                                allVsiteBondeds.bonds.size(),
                                allVsiteBondeds.angles.size(),
                                allVsiteBondeds.dihedrals.size(),
                                param.ai() + 1,
                                interaction_function[ftype].longname);
                        print_bad(debug,
                                  allVsiteBondeds.bonds,
                                  allVsiteBondeds.angles,
                                  allVsiteBondeds.dihedrals);
                    } /* debug */
                    switch (ftype)
                    {
                        case F_VSITE3:
                            bERROR = calc_vsite3_param(
                                    atypes, &param, atoms, allVsiteBondeds.bonds, allVsiteBondeds.angles);
                            break;
                        case F_VSITE3FD:
                            bERROR = calc_vsite3fd_param(
                                    &param, allVsiteBondeds.bonds, allVsiteBondeds.angles);
                            break;
                        case F_VSITE3FAD:
                            bERROR = calc_vsite3fad_param(
                                    &param, allVsiteBondeds.bonds, allVsiteBondeds.angles);
                            break;
                        case F_VSITE3OUT:
                            bERROR = calc_vsite3out_param(
                                    atypes, &param, atoms, allVsiteBondeds.bonds, allVsiteBondeds.angles);
                            break;
                        case F_VSITE4FD:
                            bERROR = calc_vsite4fd_param(
                                    &param, allVsiteBondeds.bonds, allVsiteBondeds.angles, logger);
                            break;
                        case F_VSITE4FDN:
                            bERROR = calc_vsite4fdn_param(
                                    &param, allVsiteBondeds.bonds, allVsiteBondeds.angles, logger);
                            break;
                        default:
                            gmx_fatal(FARGS,
                                      "Automatic parameter generation not supported "
                                      "for %s atom %d",
                                      interaction_function[ftype].longname,
                                      param.ai() + 1);
                            bERROR = TRUE;
                    } /* switch */
                    if (bERROR)
                    {
                        gmx_fatal(FARGS,
                                  "Automatic parameter generation not supported "
                                  "for %s atom %d for this bonding configuration",
                                  interaction_function[ftype].longname,
                                  param.ai() + 1);
                    }
                } /* if bSet */
                i++;
            }
        } /* if IF_VSITE */
    }
    return nvsite;
}

void set_vsites_ptype(bool bVerbose, gmx_moltype_t* molt, const gmx::MDLogger& logger)
{
    int ftype, i;

    if (bVerbose)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Setting particle type to V for virtual sites");
    }
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        InteractionList* il = &molt->ilist[ftype];
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            const int                nra = interaction_function[ftype].nratoms;
            const int                nrd = il->size();
            gmx::ArrayRef<const int> ia  = il->iatoms;

            if (debug && nrd)
            {
                GMX_LOG(logger.info)
                        .asParagraph()
                        .appendTextFormatted("doing %d %s virtual sites",
                                             (nrd / (nra + 1)),
                                             interaction_function[ftype].longname);
            }

            for (i = 0; (i < nrd);)
            {
                /* The virtual site */
                int avsite                     = ia[i + 1];
                molt->atoms.atom[avsite].ptype = ParticleType::VSite;

                i += nra + 1;
            }
        }
    }
}

/*! \brief
 *  Convenience typedef for linking function type to parameter numbers.
 *
 *  The entries in this datastructure are valid if the particle participates in
 *  a virtual site interaction and has a valid vsite function type other than VSITEN.
 *  \todo Change to remove empty constructor when gmx::compat::optional is available.
 */
class VsiteAtomMapping
{
public:
    //! Only construct with all information in place or nothing
    VsiteAtomMapping(int functionType, int interactionIndex) :
        functionType_(functionType), interactionIndex_(interactionIndex)
    {
    }
    VsiteAtomMapping() : functionType_(-1), interactionIndex_(-1) {}
    //! Get function type.
    const int& functionType() const { return functionType_; }
    //! Get parameter number.
    const int& interactionIndex() const { return interactionIndex_; }

private:
    //! Function type for the linked parameter.
    int functionType_;
    //! The linked parameter.
    int interactionIndex_;
};

static void check_vsite_constraints(gmx::ArrayRef<InteractionsOfType> plist,
                                    int                               cftype,
                                    const int                         vsite_type[],
                                    const gmx::MDLogger&              logger)
{
    int n = 0;
    for (const auto& param : plist[cftype].interactionTypes)
    {
        gmx::ArrayRef<const int> atoms = param.atoms();
        for (int k = 0; k < 2; k++)
        {
            int atom = atoms[k];
            if (vsite_type[atom] != NOTSET)
            {
                GMX_LOG(logger.info)
                        .asParagraph()
                        .appendTextFormatted(
                                "ERROR: Cannot have constraint (%d-%d) with virtual site (%d)",
                                param.ai() + 1,
                                param.aj() + 1,
                                atom + 1);
                n++;
            }
        }
    }
    if (n)
    {
        gmx_fatal(FARGS, "There were %d virtual sites involved in constraints", n);
    }
}

static void clean_vsite_bonds(gmx::ArrayRef<InteractionsOfType>     plist,
                              gmx::ArrayRef<const VsiteAtomMapping> pindex,
                              int                                   cftype,
                              const int                             vsite_type[],
                              const gmx::MDLogger&                  logger)
{
    int                 ftype, nOut;
    int                 nconverted, nremoved;
    int                 oatom, at1, at2;
    bool                bKeep, bRemove, bAllFD;
    InteractionsOfType* ps;

    if (cftype == F_CONNBONDS)
    {
        return;
    }

    ps         = &(plist[cftype]);
    nconverted = 0;
    nremoved   = 0;
    nOut       = 0;
    for (auto bond = ps->interactionTypes.begin(); bond != ps->interactionTypes.end();)
    {
        int        vsnral      = 0;
        const int* first_atoms = nullptr;

        bKeep   = false;
        bRemove = false;
        bAllFD  = true;
        /* check if all virtual sites are constructed from the same atoms */
        int                      nvsite = 0;
        gmx::ArrayRef<const int> atoms  = bond->atoms();
        for (int k = 0; (k < 2) && !bKeep && !bRemove; k++)
        {
            /* for all atoms in the bond */
            int atom = atoms[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                nvsite++;
                bool bThisFD  = ((pindex[atom].functionType() == F_VSITE3FD)
                                || (pindex[atom].functionType() == F_VSITE3FAD)
                                || (pindex[atom].functionType() == F_VSITE4FD)
                                || (pindex[atom].functionType() == F_VSITE4FDN));
                bool bThisOUT = ((pindex[atom].functionType() == F_VSITE3OUT)
                                 && ((interaction_function[cftype].flags & IF_CONSTRAINT) != 0U));
                bAllFD        = bAllFD && bThisFD;
                if (bThisFD || bThisOUT)
                {
                    oatom = atoms[1 - k]; /* the other atom */
                    if (vsite_type[oatom] == NOTSET
                        && oatom
                                   == plist[pindex[atom].functionType()]
                                              .interactionTypes[pindex[atom].interactionIndex()]
                                              .aj())
                    {
                        /* if the other atom isn't a vsite, and it is AI */
                        bRemove = true;
                        if (bThisOUT)
                        {
                            nOut++;
                        }
                    }
                }
                if (!bRemove)
                {
                    /* TODO This fragment, and corresponding logic in
                       clean_vsite_angles and clean_vsite_dihs should
                       be refactored into a common function */
                    if (nvsite == 1)
                    {
                        /* if this is the first vsite we encounter then
                           store construction atoms */
                        /* TODO This would be nicer to implement with
                           a C++ "vector view" class" with an
                           STL-container-like interface. */
                        vsnral      = NRAL(pindex[atom].functionType()) - 1;
                        first_atoms = plist[pindex[atom].functionType()]
                                              .interactionTypes[pindex[atom].interactionIndex()]
                                              .atoms()
                                              .data()
                                      + 1;
                    }
                    else
                    {
                        GMX_ASSERT(vsnral != 0, "nvsite > 1 must have vsnral != 0");
                        GMX_ASSERT(first_atoms != nullptr,
                                   "nvsite > 1 must have first_atoms != NULL");
                        /* if it is not the first then
                           check if this vsite is constructed from the same atoms */
                        if (vsnral == NRAL(pindex[atom].functionType()) - 1)
                        {
                            for (int m = 0; (m < vsnral) && !bKeep; m++)
                            {
                                const int* atoms;

                                bool bPresent = false;
                                atoms         = plist[pindex[atom].functionType()]
                                                .interactionTypes[pindex[atom].interactionIndex()]
                                                .atoms()
                                                .data()
                                        + 1;
                                for (int n = 0; (n < vsnral) && !bPresent; n++)
                                {
                                    if (atoms[m] == first_atoms[n])
                                    {
                                        bPresent = true;
                                    }
                                }
                                if (!bPresent)
                                {
                                    bKeep = true;
                                }
                            }
                        }
                        else
                        {
                            bKeep = true;
                        }
                    }
                }
            }
        }

        if (bRemove)
        {
            bKeep = false;
        }
        else
        {
            /* if we have no virtual sites in this bond, keep it */
            if (nvsite == 0)
            {
                bKeep = true;
            }

            /* TODO This loop and the corresponding loop in
               check_vsite_angles should be refactored into a common
               function */
            /* check if all non-vsite atoms are used in construction: */
            bool bFirstTwo = true;
            for (int k = 0; (k < 2) && !bKeep; k++) /* for all atoms in the bond */
            {
                int atom = atoms[k];
                if (vsite_type[atom] == NOTSET)
                {
                    bool bUsed = false;
                    for (int m = 0; (m < vsnral) && !bUsed; m++)
                    {
                        GMX_ASSERT(first_atoms != nullptr,
                                   "If we've seen a vsite before, we know what its first atom "
                                   "index was");

                        if (atom == first_atoms[m])
                        {
                            bUsed     = true;
                            bFirstTwo = bFirstTwo && m < 2;
                        }
                    }
                    if (!bUsed)
                    {
                        bKeep = true;
                    }
                }
            }

            if (!(bAllFD && bFirstTwo))
            {
                /* Two atom bonded interactions include constraints.
                 * We need to remove constraints between vsite pairs that have
                 * a fixed distance due to being constructed from the same
                 * atoms, since this can be numerically unstable.
                 */
                for (int m = 0; m < vsnral && !bKeep; m++) /* all constr. atoms */
                {
                    at1           = first_atoms[m];
                    at2           = first_atoms[(m + 1) % vsnral];
                    bool bPresent = false;
                    for (ftype = 0; ftype < F_NRE; ftype++)
                    {
                        if (interaction_function[ftype].flags & IF_CONSTRAINT)
                        {
                            for (auto entry = plist[ftype].interactionTypes.begin();
                                 (entry != plist[ftype].interactionTypes.end()) && !bPresent;
                                 entry++)
                            {
                                /* all constraints until one matches */
                                bPresent = (((entry->ai() == at1) && (entry->aj() == at2))
                                            || ((entry->ai() == at2) && (entry->aj() == at1)));
                            }
                        }
                    }
                    if (!bPresent)
                    {
                        bKeep = true;
                    }
                }
            }
        }

        if (bKeep)
        {
            ++bond;
        }
        else if (IS_CHEMBOND(cftype))
        {
            plist[F_CONNBONDS].interactionTypes.emplace_back(*bond);
            bond = ps->interactionTypes.erase(bond);
            nconverted++;
        }
        else
        {
            bond = ps->interactionTypes.erase(bond);
            nremoved++;
        }
    }

    if (nremoved)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Removed   %4d %15ss with virtual sites, %zu left",
                                     nremoved,
                                     interaction_function[cftype].longname,
                                     ps->size());
    }
    if (nconverted)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "Converted %4d %15ss with virtual sites to connections, %zu left",
                        nconverted,
                        interaction_function[cftype].longname,
                        ps->size());
    }
    if (nOut)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "Warning: removed %d %ss with vsite with %s construction\n"
                        "         This vsite construction does not guarantee constant "
                        "bond-length\n"
                        "         If the constructions were generated by pdb2gmx ignore "
                        "this warning",
                        nOut,
                        interaction_function[cftype].longname,
                        interaction_function[F_VSITE3OUT].longname);
    }
}

static void clean_vsite_angles(gmx::ArrayRef<InteractionsOfType>         plist,
                               gmx::ArrayRef<VsiteAtomMapping>           pindex,
                               int                                       cftype,
                               const int                                 vsite_type[],
                               gmx::ArrayRef<const Atom2VsiteConnection> at2vc,
                               const gmx::MDLogger&                      logger)
{
    int                 atom, at1, at2;
    InteractionsOfType* ps;

    ps          = &(plist[cftype]);
    int oldSize = ps->size();
    for (auto angle = ps->interactionTypes.begin(); angle != ps->interactionTypes.end();)
    {
        int        vsnral      = 0;
        const int* first_atoms = nullptr;

        bool bKeep    = false;
        bool bAll3FAD = true;
        /* check if all virtual sites are constructed from the same atoms */
        int                      nvsite = 0;
        gmx::ArrayRef<const int> atoms  = angle->atoms();
        for (int k = 0; (k < 3) && !bKeep; k++) /* for all atoms in the angle */
        {
            int atom = atoms[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                nvsite++;
                bAll3FAD = bAll3FAD && (pindex[atom].functionType() == F_VSITE3FAD);
                if (nvsite == 1)
                {
                    /* store construction atoms of first vsite */
                    vsnral      = NRAL(pindex[atom].functionType()) - 1;
                    first_atoms = plist[pindex[atom].functionType()]
                                          .interactionTypes[pindex[atom].interactionIndex()]
                                          .atoms()
                                          .data()
                                  + 1;
                }
                else
                {
                    GMX_ASSERT(vsnral != 0,
                               "If we've seen a vsite before, we know how many constructing atoms "
                               "it had");
                    GMX_ASSERT(
                            first_atoms != nullptr,
                            "If we've seen a vsite before, we know what its first atom index was");
                    /* check if this vsite is constructed from the same atoms */
                    if (vsnral == NRAL(pindex[atom].functionType()) - 1)
                    {
                        for (int m = 0; (m < vsnral) && !bKeep; m++)
                        {
                            const int* subAtoms;

                            bool bPresent = false;
                            subAtoms      = plist[pindex[atom].functionType()]
                                               .interactionTypes[pindex[atom].interactionIndex()]
                                               .atoms()
                                               .data()
                                       + 1;
                            for (int n = 0; (n < vsnral) && !bPresent; n++)
                            {
                                if (subAtoms[m] == first_atoms[n])
                                {
                                    bPresent = true;
                                }
                            }
                            if (!bPresent)
                            {
                                bKeep = true;
                            }
                        }
                    }
                    else
                    {
                        bKeep = true;
                    }
                }
            }
        }

        /* keep all angles with no virtual sites in them or
           with virtual sites with more than 3 constr. atoms */
        if (nvsite == 0 && vsnral > 3)
        {
            bKeep = true;
        }

        /* check if all non-vsite atoms are used in construction: */
        bool bFirstTwo = true;
        for (int k = 0; (k < 3) && !bKeep; k++) /* for all atoms in the angle */
        {
            atom = atoms[k];
            if (vsite_type[atom] == NOTSET)
            {
                bool bUsed = false;
                for (int m = 0; (m < vsnral) && !bUsed; m++)
                {
                    GMX_ASSERT(
                            first_atoms != nullptr,
                            "If we've seen a vsite before, we know what its first atom index was");

                    if (atom == first_atoms[m])
                    {
                        bUsed     = true;
                        bFirstTwo = bFirstTwo && m < 2;
                    }
                }
                if (!bUsed)
                {
                    bKeep = true;
                }
            }
        }

        if (!(bAll3FAD && bFirstTwo))
        {
            /* check if all constructing atoms are constrained together */
            for (int m = 0; m < vsnral && !bKeep; m++) /* all constr. atoms */
            {
                at1           = first_atoms[m];
                at2           = first_atoms[(m + 1) % vsnral];
                bool bPresent = false;
                auto found    = std::find(at2vc[at1].begin(), at2vc[at1].end(), at2);
                if (found != at2vc[at1].end())
                {
                    bPresent = true;
                }
                if (!bPresent)
                {
                    bKeep = true;
                }
            }
        }

        if (bKeep)
        {
            ++angle;
        }
        else
        {
            angle = ps->interactionTypes.erase(angle);
        }
    }

    if (oldSize != gmx::ssize(*ps))
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Removed   %4zu %15ss with virtual sites, %zu left",
                                     oldSize - ps->size(),
                                     interaction_function[cftype].longname,
                                     ps->size());
    }
}

static void clean_vsite_dihs(gmx::ArrayRef<InteractionsOfType>     plist,
                             gmx::ArrayRef<const VsiteAtomMapping> pindex,
                             int                                   cftype,
                             const int                             vsite_type[],
                             const gmx::MDLogger&                  logger)
{
    InteractionsOfType* ps;

    ps = &(plist[cftype]);

    int oldSize = ps->size();
    for (auto dih = ps->interactionTypes.begin(); dih != ps->interactionTypes.end();)
    {
        int        vsnral      = 0;
        const int* first_atoms = nullptr;
        int        atom;

        gmx::ArrayRef<const int> atoms = dih->atoms();
        bool                     bKeep = false;
        /* check if all virtual sites are constructed from the same atoms */
        int nvsite = 0;
        for (int k = 0; (k < 4) && !bKeep; k++) /* for all atoms in the dihedral */
        {
            atom = atoms[k];
            if (vsite_type[atom] != NOTSET && vsite_type[atom] != F_VSITEN)
            {
                if (nvsite == 0)
                {
                    /* store construction atoms of first vsite */
                    vsnral      = NRAL(pindex[atom].functionType()) - 1;
                    first_atoms = plist[pindex[atom].functionType()]
                                          .interactionTypes[pindex[atom].interactionIndex()]
                                          .atoms()
                                          .data()
                                  + 1;
                }
                else
                {
                    GMX_ASSERT(vsnral != 0,
                               "If we've seen a vsite before, we know how many constructing atoms "
                               "it had");
                    GMX_ASSERT(
                            first_atoms != nullptr,
                            "If we've seen a vsite before, we know what its first atom index was");
                    /* check if this vsite is constructed from the same atoms */
                    if (vsnral == NRAL(pindex[atom].functionType()) - 1)
                    {
                        for (int m = 0; (m < vsnral) && !bKeep; m++)
                        {
                            const int* subAtoms;

                            bool bPresent = false;
                            subAtoms      = plist[pindex[atom].functionType()]
                                               .interactionTypes[pindex[atom].interactionIndex()]
                                               .atoms()
                                               .data()
                                       + 1;
                            for (int n = 0; (n < vsnral) && !bPresent; n++)
                            {
                                if (subAtoms[m] == first_atoms[n])
                                {
                                    bPresent = true;
                                }
                            }
                            if (!bPresent)
                            {
                                bKeep = true;
                            }
                        }
                    }
                }
                /* TODO clean_site_bonds and _angles do this increment
                   at the top of the loop. Refactor this for
                   consistency */
                nvsite++;
            }
        }

        /* keep all dihedrals with no virtual sites in them */
        if (nvsite == 0)
        {
            bKeep = true;
        }

        /* check if all atoms in dihedral are either virtual sites, or used in
           construction of virtual sites. If so, keep it, if not throw away: */
        for (int k = 0; (k < 4) && !bKeep; k++) /* for all atoms in the dihedral */
        {
            GMX_ASSERT(vsnral != 0,
                       "If we've seen a vsite before, we know how many constructing atoms it had");
            GMX_ASSERT(first_atoms != nullptr,
                       "If we've seen a vsite before, we know what its first atom index was");
            atom = atoms[k];
            if (vsite_type[atom] == NOTSET)
            {
                /* vsnral will be set here, we don't get here with nvsite==0 */
                bool bUsed = false;
                for (int m = 0; (m < vsnral) && !bUsed; m++)
                {
                    if (atom == first_atoms[m])
                    {
                        bUsed = true;
                    }
                }
                if (!bUsed)
                {
                    bKeep = true;
                }
            }
        }

        if (bKeep)
        {
            ++dih;
        }
        else
        {
            dih = ps->interactionTypes.erase(dih);
        }
    }

    if (oldSize != gmx::ssize(*ps))
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Removed   %4zu %15ss with virtual sites, %zu left",
                                     oldSize - ps->size(),
                                     interaction_function[cftype].longname,
                                     ps->size());
    }
}

// TODO use gmx::compat::optional for pindex.
void clean_vsite_bondeds(gmx::ArrayRef<InteractionsOfType> plist,
                         int                               natoms,
                         bool                              bRmVSiteBds,
                         const gmx::MDLogger&              logger)
{
    int                               nvsite, vsite;
    int*                              vsite_type;
    std::vector<VsiteAtomMapping>     pindex;
    std::vector<Atom2VsiteConnection> at2vc;

    /* make vsite_type array */
    snew(vsite_type, natoms);
    for (int i = 0; i < natoms; i++)
    {
        vsite_type[i] = NOTSET;
    }
    nvsite = 0;
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nvsite += plist[ftype].size();
            int i = 0;
            while (i < gmx::ssize(plist[ftype]))
            {
                vsite = plist[ftype].interactionTypes[i].ai();
                if (vsite_type[vsite] == NOTSET)
                {
                    vsite_type[vsite] = ftype;
                }
                else
                {
                    gmx_fatal(FARGS, "multiple vsite constructions for atom %d", vsite + 1);
                }
                if (ftype == F_VSITEN)
                {
                    while (i < gmx::ssize(plist[ftype]) && plist[ftype].interactionTypes[i].ai() == vsite)
                    {
                        i++;
                    }
                }
                else
                {
                    i++;
                }
            }
        }
    }

    /* the rest only if we have virtual sites: */
    if (nvsite)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Cleaning up constraints %swith virtual sites",
                                     bRmVSiteBds ? "and constant bonded interactions " : "");

        /* Make a reverse list to avoid ninteractions^2 operations */
        at2vc = make_at2vsitecon(natoms, plist);

        pindex.resize(natoms);
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            /* Here we skip VSITEN. In neary all practical use cases this
             * is not an issue, since VSITEN is intended for constructing
             * additional interaction sites, not for replacing normal atoms
             * with bonded interactions. Thus we do not expect constant
             * bonded interactions. If VSITEN does get used with constant
             * bonded interactions, these are not removed which only leads
             * to very minor extra computation and constant energy.
             * The only problematic case is onstraints between VSITEN
             * constructions with fixed distance (which is anyhow useless).
             * This will generate a fatal error in check_vsite_constraints.
             */
            if ((interaction_function[ftype].flags & IF_VSITE) && ftype != F_VSITEN)
            {
                for (gmx::Index interactionIndex = 0; interactionIndex < gmx::ssize(plist[ftype]);
                     interactionIndex++)
                {
                    int k     = plist[ftype].interactionTypes[interactionIndex].ai();
                    pindex[k] = VsiteAtomMapping(ftype, interactionIndex);
                }
            }
        }

        /* remove interactions that include virtual sites */
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (((interaction_function[ftype].flags & IF_BOND) && bRmVSiteBds)
                || (interaction_function[ftype].flags & IF_CONSTRAINT))
            {
                if (interaction_function[ftype].flags & (IF_BTYPE | IF_CONSTRAINT))
                {
                    clean_vsite_bonds(plist, pindex, ftype, vsite_type, logger);
                }
                else if (interaction_function[ftype].flags & IF_ATYPE)
                {
                    clean_vsite_angles(plist, pindex, ftype, vsite_type, at2vc, logger);
                }
                else if ((ftype == F_PDIHS) || (ftype == F_IDIHS))
                {
                    clean_vsite_dihs(plist, pindex, ftype, vsite_type, logger);
                }
            }
        }
        /* check that no remaining constraints include virtual sites */
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_CONSTRAINT)
            {
                check_vsite_constraints(plist, ftype, vsite_type, logger);
            }
        }
    }
    sfree(vsite_type);
}
