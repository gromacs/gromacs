/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gen_vsite.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#define MAXNAME 32
#define OPENDIR     '[' /* starting sign for directive		*/
#define CLOSEDIR    ']' /* ending sign for directive		*/

typedef struct {
    char       atomtype[MAXNAME];  /* Type for the XH3/XH2 atom */
    gmx_bool   isplanar;           /* If true, the atomtype above and the three connected
                                    * ones are in a planar geometry. The two next entries
                                    * are undefined in that case
                                    */
    int    nhydrogens;             /* number of connected hydrogens */
    char   nextheavytype[MAXNAME]; /* Type for the heavy atom bonded to XH2/XH3 */
    char   dummymass[MAXNAME];     /* The type of MNH* or MCH3* dummy mass to use */
} t_vsiteconf;


/* Structure to represent average bond and angles values in vsite aromatic
 * residues. Note that these are NOT necessarily the bonds and angles from the
 * forcefield; many forcefields (like Amber, OPLS) have some inherent strain in
 * 5-rings (i.e. the sum of angles is !=540, but impropers keep it planar)
 */
typedef struct {
    char resname[MAXNAME];
    int  nbonds;
    int  nangles;
    struct vsitetop_bond {
        char   atom1[MAXNAME];
        char   atom2[MAXNAME];
        float  value;
    } *bond; /* list of bonds */
    struct vsitetop_angle {
        char   atom1[MAXNAME];
        char   atom2[MAXNAME];
        char   atom3[MAXNAME];
        float  value;
    } *angle; /* list of angles */
} t_vsitetop;


enum {
    DDB_CH3, DDB_NH3, DDB_NH2, DDB_PHE, DDB_TYR,
    DDB_TRP, DDB_HISA, DDB_HISB, DDB_HISH, DDB_DIR_NR
};

typedef char t_dirname[STRLEN];

static const t_dirname ddb_dirnames[DDB_DIR_NR] = {
    "CH3",
    "NH3",
    "NH2",
    "PHE",
    "TYR",
    "TRP",
    "HISA",
    "HISB",
    "HISH"
};

static int ddb_name2dir(char *name)
{
    /* Translate a directive name to the number of the directive.
     * HID/HIE/HIP names are translated to the ones we use in Gromacs.
     */

    int i, index;

    index = -1;

    for (i = 0; i < DDB_DIR_NR && index < 0; i++)
    {
        if (!gmx_strcasecmp(name, ddb_dirnames[i]))
        {
            index = i;
        }
    }

    return index;
}


static void read_vsite_database(const char *ddbname,
                                t_vsiteconf **pvsiteconflist, int *nvsiteconf,
                                t_vsitetop **pvsitetoplist, int *nvsitetop)
{
    /* This routine is a quick hack to fix the problem with hardcoded atomtypes
     * and aromatic vsite parameters by reading them from a ff???.vsd file.
     *
     * The file can contain sections [ NH3 ], [ CH3 ], [ NH2 ], and ring residue names.
     * For the NH3 and CH3 section each line has three fields. The first is the atomtype
     * (nb: not bonded type) of the N/C atom to be replaced, the second field is
     * the type of the next heavy atom it is bonded to, and the third field the type
     * of dummy mass that will be used for this group.
     *
     * If the NH2 group planar (sp2 N) a different vsite construct is used, so in this
     * case the second field should just be the word planar.
     */

    FILE        *ddb;
    char         dirstr[STRLEN];
    char         pline[STRLEN];
    int          i, n, k, nvsite, ntop, curdir;
    t_vsiteconf *vsiteconflist;
    t_vsitetop  *vsitetoplist;
    char        *ch;
    char         s1[MAXNAME], s2[MAXNAME], s3[MAXNAME], s4[MAXNAME];

    ddb = libopen(ddbname);

    nvsite        = *nvsiteconf;
    vsiteconflist = *pvsiteconflist;
    ntop          = *nvsitetop;
    vsitetoplist  = *pvsitetoplist;

    curdir = -1;

    snew(vsiteconflist, 1);
    snew(vsitetoplist, 1);

    while (fgets2(pline, STRLEN-2, ddb) != nullptr)
    {
        strip_comment(pline);
        trim(pline);
        if (strlen(pline) > 0)
        {
            if (pline[0] == OPENDIR)
            {
                strncpy(dirstr, pline+1, STRLEN-2);
                if ((ch = strchr (dirstr, CLOSEDIR)) != nullptr)
                {
                    (*ch) = 0;
                }
                trim (dirstr);

                if (!gmx_strcasecmp(dirstr, "HID") ||
                    !gmx_strcasecmp(dirstr, "HISD"))
                {
                    sprintf(dirstr, "HISA");
                }
                else if (!gmx_strcasecmp(dirstr, "HIE") ||
                         !gmx_strcasecmp(dirstr, "HISE"))
                {
                    sprintf(dirstr, "HISB");
                }
                else if (!gmx_strcasecmp(dirstr, "HIP"))
                {
                    sprintf(dirstr, "HISH");
                }

                curdir = ddb_name2dir(dirstr);
                if (curdir < 0)
                {
                    gmx_fatal(FARGS, "Invalid directive %s in vsite database %s",
                              dirstr, ddbname);
                }
            }
            else
            {
                switch (curdir)
                {
                    case -1:
                        gmx_fatal(FARGS, "First entry in vsite database must be a directive.\n");
                        break;
                    case DDB_CH3:
                    case DDB_NH3:
                    case DDB_NH2:
                        n = sscanf(pline, "%s%s%s", s1, s2, s3);
                        if (n < 3 && !gmx_strcasecmp(s2, "planar"))
                        {
                            srenew(vsiteconflist, nvsite+1);
                            strncpy(vsiteconflist[nvsite].atomtype, s1, MAXNAME-1);
                            vsiteconflist[nvsite].isplanar         = TRUE;
                            vsiteconflist[nvsite].nextheavytype[0] = 0;
                            vsiteconflist[nvsite].dummymass[0]     = 0;
                            vsiteconflist[nvsite].nhydrogens       = 2;
                            nvsite++;
                        }
                        else if (n == 3)
                        {
                            srenew(vsiteconflist, (nvsite+1));
                            strncpy(vsiteconflist[nvsite].atomtype, s1, MAXNAME-1);
                            vsiteconflist[nvsite].isplanar = FALSE;
                            strncpy(vsiteconflist[nvsite].nextheavytype, s2, MAXNAME-1);
                            strncpy(vsiteconflist[nvsite].dummymass, s3, MAXNAME-1);
                            if (curdir == DDB_NH2)
                            {
                                vsiteconflist[nvsite].nhydrogens = 2;
                            }
                            else
                            {
                                vsiteconflist[nvsite].nhydrogens = 3;
                            }
                            nvsite++;
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Not enough directives in vsite database line: %s\n", pline);
                        }
                        break;
                    case DDB_PHE:
                    case DDB_TYR:
                    case DDB_TRP:
                    case DDB_HISA:
                    case DDB_HISB:
                    case DDB_HISH:
                        i = 0;
                        while ((i < ntop) && gmx_strcasecmp(dirstr, vsitetoplist[i].resname))
                        {
                            i++;
                        }
                        /* Allocate a new topology entry if this is a new residue */
                        if (i == ntop)
                        {
                            srenew(vsitetoplist, ntop+1);
                            ntop++; /* i still points to current vsite topology entry */
                            strncpy(vsitetoplist[i].resname, dirstr, MAXNAME-1);
                            vsitetoplist[i].nbonds = vsitetoplist[i].nangles = 0;
                            snew(vsitetoplist[i].bond, 1);
                            snew(vsitetoplist[i].angle, 1);
                        }
                        n = sscanf(pline, "%s%s%s%s", s1, s2, s3, s4);
                        if (n == 3)
                        {
                            /* bond */
                            k = vsitetoplist[i].nbonds++;
                            srenew(vsitetoplist[i].bond, k+1);
                            strncpy(vsitetoplist[i].bond[k].atom1, s1, MAXNAME-1);
                            strncpy(vsitetoplist[i].bond[k].atom2, s2, MAXNAME-1);
                            vsitetoplist[i].bond[k].value = strtod(s3, nullptr);
                        }
                        else if (n == 4)
                        {
                            /* angle */
                            k = vsitetoplist[i].nangles++;
                            srenew(vsitetoplist[i].angle, k+1);
                            strncpy(vsitetoplist[i].angle[k].atom1, s1, MAXNAME-1);
                            strncpy(vsitetoplist[i].angle[k].atom2, s2, MAXNAME-1);
                            strncpy(vsitetoplist[i].angle[k].atom3, s3, MAXNAME-1);
                            vsitetoplist[i].angle[k].value = strtod(s4, nullptr);
                        }
                        else
                        {
                            gmx_fatal(FARGS, "Need 3 or 4 values to specify bond/angle values in %s: %s\n", ddbname, pline);
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "Didnt find a case for directive %s in read_vsite_database\n", dirstr);
                }
            }
        }
    }

    *pvsiteconflist = vsiteconflist;
    *pvsitetoplist  = vsitetoplist;
    *nvsiteconf     = nvsite;
    *nvsitetop      = ntop;

    gmx_ffclose(ddb);
}

static int nitrogen_is_planar(t_vsiteconf vsiteconflist[], int nvsiteconf, char atomtype[])
{
    /* Return 1 if atomtype exists in database list and is planar, 0 if not,
     * and -1 if not found.
     */
    int      i, res;
    gmx_bool found = FALSE;
    for (i = 0; i < nvsiteconf && !found; i++)
    {
        found = (!gmx_strcasecmp(vsiteconflist[i].atomtype, atomtype) && (vsiteconflist[i].nhydrogens == 2));
    }
    if (found)
    {
        res = (vsiteconflist[i-1].isplanar == TRUE);
    }
    else
    {
        res = -1;
    }

    return res;
}

static char *get_dummymass_name(t_vsiteconf vsiteconflist[], int nvsiteconf, char atom[], char nextheavy[])
{
    /* Return the dummy mass name if found, or NULL if not set in ddb database */
    int      i;
    gmx_bool found = FALSE;
    for (i = 0; i < nvsiteconf && !found; i++)
    {
        found = (!gmx_strcasecmp(vsiteconflist[i].atomtype, atom) &&
                 !gmx_strcasecmp(vsiteconflist[i].nextheavytype, nextheavy));
    }
    if (found)
    {
        return vsiteconflist[i-1].dummymass;
    }
    else
    {
        return nullptr;
    }
}



static real get_ddb_bond(t_vsitetop *vsitetop, int nvsitetop,
                         const char res[],
                         const char atom1[], const char atom2[])
{
    int i, j;

    i = 0;
    while (i < nvsitetop && gmx_strcasecmp(res, vsitetop[i].resname))
    {
        i++;
    }
    if (i == nvsitetop)
    {
        gmx_fatal(FARGS, "No vsite information for residue %s found in vsite database.\n", res);
    }
    j = 0;
    while (j < vsitetop[i].nbonds &&
           ( strcmp(atom1, vsitetop[i].bond[j].atom1) || strcmp(atom2, vsitetop[i].bond[j].atom2)) &&
           ( strcmp(atom2, vsitetop[i].bond[j].atom1) || strcmp(atom1, vsitetop[i].bond[j].atom2)))
    {
        j++;
    }
    if (j == vsitetop[i].nbonds)
    {
        gmx_fatal(FARGS, "Couldnt find bond %s-%s for residue %s in vsite database.\n", atom1, atom2, res);
    }

    return vsitetop[i].bond[j].value;
}


static real get_ddb_angle(t_vsitetop *vsitetop, int nvsitetop,
                          const char res[], const char atom1[],
                          const char atom2[], const char atom3[])
{
    int i, j;

    i = 0;
    while (i < nvsitetop && gmx_strcasecmp(res, vsitetop[i].resname))
    {
        i++;
    }
    if (i == nvsitetop)
    {
        gmx_fatal(FARGS, "No vsite information for residue %s found in vsite database.\n", res);
    }
    j = 0;
    while (j < vsitetop[i].nangles &&
           ( strcmp(atom1, vsitetop[i].angle[j].atom1) ||
             strcmp(atom2, vsitetop[i].angle[j].atom2) ||
             strcmp(atom3, vsitetop[i].angle[j].atom3)) &&
           ( strcmp(atom3, vsitetop[i].angle[j].atom1) ||
             strcmp(atom2, vsitetop[i].angle[j].atom2) ||
             strcmp(atom1, vsitetop[i].angle[j].atom3)))
    {
        j++;
    }
    if (j == vsitetop[i].nangles)
    {
        gmx_fatal(FARGS, "Couldnt find angle %s-%s-%s for residue %s in vsite database.\n", atom1, atom2, atom3, res);
    }

    return vsitetop[i].angle[j].value;
}


static void count_bonds(int atom, t_params *psb, char ***atomname,
                        int *nrbonds, int *nrHatoms, int Hatoms[], int *Heavy,
                        int *nrheavies, int heavies[])
{
    int i, heavy, other, nrb, nrH, nrhv;

    /* find heavy atom bound to this hydrogen */
    heavy = NOTSET;
    for (i = 0; (i < psb->nr) && (heavy == NOTSET); i++)
    {
        if (psb->param[i].ai() == atom)
        {
            heavy = psb->param[i].aj();
        }
        else if (psb->param[i].aj() == atom)
        {
            heavy = psb->param[i].ai();
        }
    }
    if (heavy == NOTSET)
    {
        gmx_fatal(FARGS, "unbound hydrogen atom %d", atom+1);
    }
    /* find all atoms bound to heavy atom */
    other = NOTSET;
    nrb   = 0;
    nrH   = 0;
    nrhv  = 0;
    for (i = 0; i < psb->nr; i++)
    {
        if (psb->param[i].ai() == heavy)
        {
            other = psb->param[i].aj();
        }
        else if (psb->param[i].aj() == heavy)
        {
            other = psb->param[i].ai();
        }
        if (other != NOTSET)
        {
            nrb++;
            if (is_hydrogen(*(atomname[other])))
            {
                Hatoms[nrH] = other;
                nrH++;
            }
            else
            {
                heavies[nrhv] = other;
                nrhv++;
            }
            other = NOTSET;
        }
    }
    *Heavy     = heavy;
    *nrbonds   = nrb;
    *nrHatoms  = nrH;
    *nrheavies = nrhv;
}

static void print_bonds(FILE *fp, int o2n[],
                        int nrHatoms, int Hatoms[], int Heavy,
                        int nrheavies, int heavies[])
{
    int i;

    fprintf(fp, "Found: %d Hatoms: ", nrHatoms);
    for (i = 0; i < nrHatoms; i++)
    {
        fprintf(fp, " %d", o2n[Hatoms[i]]+1);
    }
    fprintf(fp, "; %d Heavy atoms: %d", nrheavies+1, o2n[Heavy]+1);
    for (i = 0; i < nrheavies; i++)
    {
        fprintf(fp, " %d", o2n[heavies[i]]+1);
    }
    fprintf(fp, "\n");
}

static int get_atype(int atom, t_atoms *at, int nrtp, t_restp rtp[],
                     gmx_residuetype_t *rt)
{
    int      type;
    gmx_bool bNterm;
    int      j;
    t_restp *rtpp;

    if (at->atom[atom].m)
    {
        type = at->atom[atom].type;
    }
    else
    {
        /* get type from rtp */
        rtpp   = get_restp(*(at->resinfo[at->atom[atom].resind].name), nrtp, rtp);
        bNterm = gmx_residuetype_is_protein(rt, *(at->resinfo[at->atom[atom].resind].name)) &&
            (at->atom[atom].resind == 0);
        j    = search_jtype(rtpp, *(at->atomname[atom]), bNterm);
        type = rtpp->atom[j].type;
    }
    return type;
}

static int vsite_nm2type(const char *name, gpp_atomtype_t atype)
{
    int tp;

    tp = get_atomtype_type(name, atype);
    if (tp == NOTSET)
    {
        gmx_fatal(FARGS, "Dummy mass type (%s) not found in atom type database",
                  name);
    }

    return tp;
}

static real get_amass(int atom, t_atoms *at, int nrtp, t_restp rtp[],
                      gmx_residuetype_t *rt)
{
    real     mass;
    gmx_bool bNterm;
    int      j;
    t_restp *rtpp;

    if (at->atom[atom].m)
    {
        mass = at->atom[atom].m;
    }
    else
    {
        /* get mass from rtp */
        rtpp   = get_restp(*(at->resinfo[at->atom[atom].resind].name), nrtp, rtp);
        bNterm = gmx_residuetype_is_protein(rt, *(at->resinfo[at->atom[atom].resind].name)) &&
            (at->atom[atom].resind == 0);
        j    = search_jtype(rtpp, *(at->atomname[atom]), bNterm);
        mass = rtpp->atom[j].m;
    }
    return mass;
}

static void my_add_param(t_params *plist, int ai, int aj, real b)
{
    static real c[MAXFORCEPARAM] =
    { NOTSET, NOTSET, NOTSET, NOTSET, NOTSET, NOTSET };

    c[0] = b;
    add_param(plist, ai, aj, c, nullptr);
}

static void add_vsites(t_params plist[], int vsite_type[],
                       int Heavy, int nrHatoms, int Hatoms[],
                       int nrheavies, int heavies[])
{
    int      i, j, ftype, other, moreheavy;
    gmx_bool bSwapParity;

    for (i = 0; i < nrHatoms; i++)
    {
        ftype = vsite_type[Hatoms[i]];
        /* Errors in setting the vsite_type should really be caugth earlier,
         * because here it's not possible to print any useful error message.
         * But it's still better to print a message than to segfault.
         */
        if (ftype == NOTSET)
        {
            gmx_incons("Undetected error in setting up virtual sites");
        }
        bSwapParity           = (ftype < 0);
        vsite_type[Hatoms[i]] = ftype = abs(ftype);
        if (ftype == F_BONDS)
        {
            if ( (nrheavies != 1) && (nrHatoms != 1) )
            {
                gmx_fatal(FARGS, "cannot make constraint in add_vsites for %d heavy "
                          "atoms and %d hydrogen atoms", nrheavies, nrHatoms);
            }
            my_add_param(&(plist[F_CONSTRNC]), Hatoms[i], heavies[0], NOTSET);
        }
        else
        {
            switch (ftype)
            {
                case F_VSITE3:
                case F_VSITE3FD:
                case F_VSITE3OUT:
                    if (nrheavies < 2)
                    {
                        gmx_fatal(FARGS, "Not enough heavy atoms (%d) for %s (min 3)",
                                  nrheavies+1,
                                  interaction_function[vsite_type[Hatoms[i]]].name);
                    }
                    add_vsite3_atoms(&plist[ftype], Hatoms[i], Heavy, heavies[0], heavies[1],
                                     bSwapParity);
                    break;
                case F_VSITE3FAD:
                {
                    if (nrheavies > 1)
                    {
                        moreheavy = heavies[1];
                    }
                    else
                    {
                        /* find more heavy atoms */
                        other = moreheavy = NOTSET;
                        for (j = 0; (j < plist[F_BONDS].nr) && (moreheavy == NOTSET); j++)
                        {
                            if (plist[F_BONDS].param[j].ai() == heavies[0])
                            {
                                other = plist[F_BONDS].param[j].aj();
                            }
                            else if (plist[F_BONDS].param[j].aj() == heavies[0])
                            {
                                other = plist[F_BONDS].param[j].ai();
                            }
                            if ( (other != NOTSET) && (other != Heavy) )
                            {
                                moreheavy = other;
                            }
                        }
                        if (moreheavy == NOTSET)
                        {
                            gmx_fatal(FARGS, "Unbound molecule part %d-%d", Heavy+1, Hatoms[0]+1);
                        }
                    }
                    add_vsite3_atoms(&plist[ftype], Hatoms[i], Heavy, heavies[0], moreheavy,
                                     bSwapParity);
                    break;
                }
                case F_VSITE4FD:
                case F_VSITE4FDN:
                    if (nrheavies < 3)
                    {
                        gmx_fatal(FARGS, "Not enough heavy atoms (%d) for %s (min 4)",
                                  nrheavies+1,
                                  interaction_function[vsite_type[Hatoms[i]]].name);
                    }
                    add_vsite4_atoms(&plist[ftype],
                                     Hatoms[0], Heavy, heavies[0], heavies[1], heavies[2]);
                    break;

                default:
                    gmx_fatal(FARGS, "can't use add_vsites for interaction function %s",
                              interaction_function[vsite_type[Hatoms[i]]].name);
            } /* switch ftype */
        }     /* else */
    }         /* for i */
}

#define ANGLE_6RING (DEG2RAD*120)

/* cosine rule: a^2 = b^2 + c^2 - 2 b c cos(alpha) */
/* get a^2 when a, b and alpha are given: */
#define cosrule(b, c, alpha) ( gmx::square(b) + gmx::square(c) - 2*b*c*std::cos(alpha) )
/* get cos(alpha) when a, b and c are given: */
#define acosrule(a, b, c) ( (gmx::square(b)+gmx::square(c)-gmx::square(a))/(2*b*c) )

static int gen_vsites_6ring(t_atoms *at, int *vsite_type[], t_params plist[],
                            int nrfound, int *ats, real bond_cc, real bond_ch,
                            real xcom, gmx_bool bDoZ)
{
    /* these MUST correspond to the atnms array in do_vsite_aromatics! */
    enum {
        atCG, atCD1, atHD1, atCD2, atHD2, atCE1, atHE1, atCE2, atHE2,
        atCZ, atHZ, atNR
    };

    int  i, nvsite;
    real a, b, dCGCE, tmp1, tmp2, mtot, mG, mrest;
    real xCG;
    /* CG, CE1 and CE2 stay and each get a part of the total mass,
     * so the c-o-m stays the same.
     */

    if (bDoZ)
    {
        if (atNR != nrfound)
        {
            gmx_incons("Generating vsites on 6-rings");
        }
    }

    /* constraints between CG, CE1 and CE2: */
    dCGCE = std::sqrt( cosrule(bond_cc, bond_cc, ANGLE_6RING) );
    my_add_param(&(plist[F_CONSTRNC]), ats[atCG], ats[atCE1], dCGCE);
    my_add_param(&(plist[F_CONSTRNC]), ats[atCG], ats[atCE2], dCGCE);
    my_add_param(&(plist[F_CONSTRNC]), ats[atCE1], ats[atCE2], dCGCE);

    /* rest will be vsite3 */
    mtot   = 0;
    nvsite = 0;
    for (i = 0; i <  (bDoZ ? atNR : atHZ); i++)
    {
        mtot += at->atom[ats[i]].m;
        if (i != atCG && i != atCE1 && i != atCE2 && (bDoZ || (i != atHZ && i != atCZ) ) )
        {
            at->atom[ats[i]].m    = at->atom[ats[i]].mB = 0;
            (*vsite_type)[ats[i]] = F_VSITE3;
            nvsite++;
        }
    }
    /* Distribute mass so center-of-mass stays the same.
     * The center-of-mass in the call is defined with x=0 at
     * the CE1-CE2 bond and y=0 at the line from CG to the middle of CE1-CE2 bond.
     */
    xCG  = -bond_cc+bond_cc*std::cos(ANGLE_6RING);

    mG                             = at->atom[ats[atCG]].m = at->atom[ats[atCG]].mB = xcom*mtot/xCG;
    mrest                          = mtot-mG;
    at->atom[ats[atCE1]].m         = at->atom[ats[atCE1]].mB =
            at->atom[ats[atCE2]].m = at->atom[ats[atCE2]].mB = mrest / 2;

    /* vsite3 construction: r_d = r_i + a r_ij + b r_ik */
    tmp1  = dCGCE*std::sin(ANGLE_6RING*0.5);
    tmp2  = bond_cc*std::cos(0.5*ANGLE_6RING) + tmp1;
    tmp1 *= 2;
    a     = b = -bond_ch / tmp1;
    /* HE1 and HE2: */
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atHE1], ats[atCE1], ats[atCE2], ats[atCG], a, b);
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atHE2], ats[atCE2], ats[atCE1], ats[atCG], a, b);
    /* CD1, CD2 and CZ: */
    a = b = tmp2 / tmp1;
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atCD1], ats[atCE2], ats[atCE1], ats[atCG], a, b);
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atCD2], ats[atCE1], ats[atCE2], ats[atCG], a, b);
    if (bDoZ)
    {
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atCZ], ats[atCG], ats[atCE1], ats[atCE2], a, b);
    }
    /* HD1, HD2 and HZ: */
    a = b = ( bond_ch + tmp2 ) / tmp1;
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atHD1], ats[atCE2], ats[atCE1], ats[atCG], a, b);
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atHD2], ats[atCE1], ats[atCE2], ats[atCG], a, b);
    if (bDoZ)
    {
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atHZ], ats[atCG], ats[atCE1], ats[atCE2], a, b);
    }

    return nvsite;
}

static int gen_vsites_phe(t_atoms *at, int *vsite_type[], t_params plist[],
                          int nrfound, int *ats, t_vsitetop *vsitetop, int nvsitetop)
{
    real bond_cc, bond_ch;
    real xcom, mtot;
    int  i;
    /* these MUST correspond to the atnms array in do_vsite_aromatics! */
    enum {
        atCG, atCD1, atHD1, atCD2, atHD2, atCE1, atHE1, atCE2, atHE2,
        atCZ, atHZ, atNR
    };
    real x[atNR];
    /* Aromatic rings have 6-fold symmetry, so we only need one bond length.
     * (angle is always 120 degrees).
     */
    bond_cc = get_ddb_bond(vsitetop, nvsitetop, "PHE", "CD1", "CE1");
    bond_ch = get_ddb_bond(vsitetop, nvsitetop, "PHE", "CD1", "HD1");

    x[atCG]  = -bond_cc+bond_cc*std::cos(ANGLE_6RING);
    x[atCD1] = -bond_cc;
    x[atHD1] = x[atCD1]+bond_ch*std::cos(ANGLE_6RING);
    x[atCE1] = 0;
    x[atHE1] = x[atCE1]-bond_ch*std::cos(ANGLE_6RING);
    x[atCD2] = x[atCD1];
    x[atHD2] = x[atHD1];
    x[atCE2] = x[atCE1];
    x[atHE2] = x[atHE1];
    x[atCZ]  = bond_cc*std::cos(0.5*ANGLE_6RING);
    x[atHZ]  = x[atCZ]+bond_ch;

    xcom = mtot = 0;
    for (i = 0; i < atNR; i++)
    {
        xcom += x[i]*at->atom[ats[i]].m;
        mtot += at->atom[ats[i]].m;
    }
    xcom /= mtot;

    return gen_vsites_6ring(at, vsite_type, plist, nrfound, ats, bond_cc, bond_ch, xcom, TRUE);
}

static void calc_vsite3_param(real xd, real yd, real xi, real yi, real xj, real yj,
                              real xk, real yk, real *a, real *b)
{
    /* determine parameters by solving the equation system, since we know the
     * virtual site coordinates here.
     */
    real dx_ij, dx_ik, dy_ij, dy_ik;

    dx_ij = xj-xi;
    dy_ij = yj-yi;
    dx_ik = xk-xi;
    dy_ik = yk-yi;

    *a = ( (xd-xi)*dy_ik - dx_ik*(yd-yi) ) / (dx_ij*dy_ik - dx_ik*dy_ij);
    *b = ( yd - yi - (*a)*dy_ij ) / dy_ik;
}


static int gen_vsites_trp(gpp_atomtype_t atype, rvec *newx[],
                          t_atom *newatom[], char ***newatomname[],
                          int *o2n[], int *newvsite_type[], int *newcgnr[],
                          t_symtab *symtab, int *nadd, rvec x[], int *cgnr[],
                          t_atoms *at, int *vsite_type[], t_params plist[],
                          int nrfound, int *ats, int add_shift,
                          t_vsitetop *vsitetop, int nvsitetop)
{
#define NMASS 2
    /* these MUST correspond to the atnms array in do_vsite_aromatics! */
    enum {
        atCB,  atCG,  atCD1, atHD1, atCD2, atNE1, atHE1, atCE2, atCE3, atHE3,
        atCZ2, atHZ2, atCZ3, atHZ3, atCH2, atHH2, atNR
    };
    /* weights for determining the COM's of both rings (M1 and M2): */
    real mw[NMASS][atNR] = {
        {   0,     1,     1,     1,   0.5,     1,     1,   0.5,     0,     0,
            0,     0,     0,     0,     0,     0 },
        {   0,     0,     0,     0,   0.5,     0,     0,   0.5,     1,     1,
            1,     1,     1,     1,     1,     1 }
    };

    real xi[atNR], yi[atNR];
    real xcom[NMASS], ycom[NMASS], alpha;
    real b_CD2_CE2, b_NE1_CE2, b_CG_CD2, b_CH2_HH2, b_CE2_CZ2;
    real b_NE1_HE1, b_CD2_CE3, b_CE3_CZ3, b_CB_CG;
    real b_CZ2_CH2, b_CZ2_HZ2, b_CD1_HD1, b_CE3_HE3;
    real b_CG_CD1, b_CZ3_HZ3;
    real a_NE1_CE2_CD2, a_CE2_CD2_CG, a_CB_CG_CD2, a_CE2_CD2_CE3;
    real a_CD2_CG_CD1, a_CE2_CZ2_HZ2, a_CZ2_CH2_HH2;
    real a_CD2_CE2_CZ2, a_CD2_CE3_CZ3, a_CE3_CZ3_HZ3, a_CG_CD1_HD1;
    real a_CE2_CZ2_CH2, a_HE1_NE1_CE2, a_CD2_CE3_HE3;
    int  atM[NMASS], tpM, i, i0, j, nvsite;
    real mM[NMASS], dCBM1, dCBM2, dM1M2;
    real a, b;
    rvec r_ij, r_ik, t1, t2;
    char name[10];

    if (atNR != nrfound)
    {
        gmx_incons("atom types in gen_vsites_trp");
    }
    /* Get geometry from database */
    b_CD2_CE2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CD2", "CE2");
    b_NE1_CE2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "NE1", "CE2");
    b_CG_CD1  = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CG", "CD1");
    b_CG_CD2  = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CG", "CD2");
    b_CB_CG   = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CB", "CG");
    b_CE2_CZ2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CE2", "CZ2");
    b_CD2_CE3 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CD2", "CE3");
    b_CE3_CZ3 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CE3", "CZ3");
    b_CZ2_CH2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CZ2", "CH2");

    b_CD1_HD1 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CD1", "HD1");
    b_CZ2_HZ2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CZ2", "HZ2");
    b_NE1_HE1 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "NE1", "HE1");
    b_CH2_HH2 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CH2", "HH2");
    b_CE3_HE3 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CE3", "HE3");
    b_CZ3_HZ3 = get_ddb_bond(vsitetop, nvsitetop, "TRP", "CZ3", "HZ3");

    a_NE1_CE2_CD2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "NE1", "CE2", "CD2");
    a_CE2_CD2_CG  = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CE2", "CD2", "CG");
    a_CB_CG_CD2   = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CB", "CG", "CD2");
    a_CD2_CG_CD1  = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CD2", "CG", "CD1");
    /*a_CB_CG_CD1   = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CB", "CG", "CD1"); unused */

    a_CE2_CD2_CE3 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CE2", "CD2", "CE3");
    a_CD2_CE2_CZ2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CD2", "CE2", "CZ2");
    a_CD2_CE3_CZ3 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CD2", "CE3", "CZ3");
    a_CE3_CZ3_HZ3 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CE3", "CZ3", "HZ3");
    a_CZ2_CH2_HH2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CZ2", "CH2", "HH2");
    a_CE2_CZ2_HZ2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CE2", "CZ2", "HZ2");
    a_CE2_CZ2_CH2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CE2", "CZ2", "CH2");
    a_CG_CD1_HD1  = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CG", "CD1", "HD1");
    a_HE1_NE1_CE2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "HE1", "NE1", "CE2");
    a_CD2_CE3_HE3 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TRP", "CD2", "CE3", "HE3");

    /* Calculate local coordinates.
     * y-axis (x=0) is the bond CD2-CE2.
     * x-axis (y=0) is perpendicular to the bond CD2-CE2 and
     * intersects the middle of the bond.
     */
    xi[atCD2] = 0;
    yi[atCD2] = -0.5*b_CD2_CE2;

    xi[atCE2] = 0;
    yi[atCE2] = 0.5*b_CD2_CE2;

    xi[atNE1] = -b_NE1_CE2*std::sin(a_NE1_CE2_CD2);
    yi[atNE1] = yi[atCE2]-b_NE1_CE2*std::cos(a_NE1_CE2_CD2);

    xi[atCG] = -b_CG_CD2*std::sin(a_CE2_CD2_CG);
    yi[atCG] = yi[atCD2]+b_CG_CD2*std::cos(a_CE2_CD2_CG);

    alpha    = a_CE2_CD2_CG + M_PI - a_CB_CG_CD2;
    xi[atCB] = xi[atCG]-b_CB_CG*std::sin(alpha);
    yi[atCB] = yi[atCG]+b_CB_CG*std::cos(alpha);

    alpha     = a_CE2_CD2_CG + a_CD2_CG_CD1 - M_PI;
    xi[atCD1] = xi[atCG]-b_CG_CD1*std::sin(alpha);
    yi[atCD1] = yi[atCG]+b_CG_CD1*std::cos(alpha);

    xi[atCE3] = b_CD2_CE3*std::sin(a_CE2_CD2_CE3);
    yi[atCE3] = yi[atCD2]+b_CD2_CE3*std::cos(a_CE2_CD2_CE3);

    xi[atCZ2] = b_CE2_CZ2*std::sin(a_CD2_CE2_CZ2);
    yi[atCZ2] = yi[atCE2]-b_CE2_CZ2*std::cos(a_CD2_CE2_CZ2);

    alpha     = a_CE2_CD2_CE3 + a_CD2_CE3_CZ3 - M_PI;
    xi[atCZ3] = xi[atCE3]+b_CE3_CZ3*std::sin(alpha);
    yi[atCZ3] = yi[atCE3]+b_CE3_CZ3*std::cos(alpha);

    alpha     = a_CD2_CE2_CZ2 + a_CE2_CZ2_CH2 - M_PI;
    xi[atCH2] = xi[atCZ2]+b_CZ2_CH2*std::sin(alpha);
    yi[atCH2] = yi[atCZ2]-b_CZ2_CH2*std::cos(alpha);

    /* hydrogens */
    alpha     = a_CE2_CD2_CG + a_CD2_CG_CD1 - a_CG_CD1_HD1;
    xi[atHD1] = xi[atCD1]-b_CD1_HD1*std::sin(alpha);
    yi[atHD1] = yi[atCD1]+b_CD1_HD1*std::cos(alpha);

    alpha     = a_NE1_CE2_CD2 + M_PI - a_HE1_NE1_CE2;
    xi[atHE1] = xi[atNE1]-b_NE1_HE1*std::sin(alpha);
    yi[atHE1] = yi[atNE1]-b_NE1_HE1*std::cos(alpha);

    alpha     = a_CE2_CD2_CE3 + M_PI - a_CD2_CE3_HE3;
    xi[atHE3] = xi[atCE3]+b_CE3_HE3*std::sin(alpha);
    yi[atHE3] = yi[atCE3]+b_CE3_HE3*std::cos(alpha);

    alpha     = a_CD2_CE2_CZ2 + M_PI - a_CE2_CZ2_HZ2;
    xi[atHZ2] = xi[atCZ2]+b_CZ2_HZ2*std::sin(alpha);
    yi[atHZ2] = yi[atCZ2]-b_CZ2_HZ2*std::cos(alpha);

    alpha     = a_CD2_CE2_CZ2 + a_CE2_CZ2_CH2 - a_CZ2_CH2_HH2;
    xi[atHZ3] = xi[atCZ3]+b_CZ3_HZ3*std::sin(alpha);
    yi[atHZ3] = yi[atCZ3]+b_CZ3_HZ3*std::cos(alpha);

    alpha     = a_CE2_CD2_CE3 + a_CD2_CE3_CZ3 - a_CE3_CZ3_HZ3;
    xi[atHH2] = xi[atCH2]+b_CH2_HH2*std::sin(alpha);
    yi[atHH2] = yi[atCH2]-b_CH2_HH2*std::cos(alpha);

    /* Calculate masses for each ring and put it on the dummy masses */
    for (j = 0; j < NMASS; j++)
    {
        mM[j] = xcom[j] = ycom[j] = 0;
    }
    for (i = 0; i < atNR; i++)
    {
        if (i != atCB)
        {
            for (j = 0; j < NMASS; j++)
            {
                mM[j]   += mw[j][i] * at->atom[ats[i]].m;
                xcom[j] += xi[i] * mw[j][i] * at->atom[ats[i]].m;
                ycom[j] += yi[i] * mw[j][i] * at->atom[ats[i]].m;
            }
        }
    }
    for (j = 0; j < NMASS; j++)
    {
        xcom[j] /= mM[j];
        ycom[j] /= mM[j];
    }

    /* get dummy mass type */
    tpM = vsite_nm2type("MW", atype);
    /* make space for 2 masses: shift all atoms starting with CB */
    i0 = ats[atCB];
    for (j = 0; j < NMASS; j++)
    {
        atM[j] = i0+*nadd+j;
    }
    if (debug)
    {
        fprintf(stderr, "Inserting %d dummy masses at %d\n", NMASS, (*o2n)[i0]+1);
    }
    *nadd += NMASS;
    for (j = i0; j < at->nr; j++)
    {
        (*o2n)[j] = j+*nadd;
    }
    srenew(*newx, at->nr+*nadd);
    srenew(*newatom, at->nr+*nadd);
    srenew(*newatomname, at->nr+*nadd);
    srenew(*newvsite_type, at->nr+*nadd);
    srenew(*newcgnr, at->nr+*nadd);
    for (j = 0; j < NMASS; j++)
    {
        (*newatomname)[at->nr+*nadd-1-j] = nullptr;
    }

    /* Dummy masses will be placed at the center-of-mass in each ring. */

    /* calc initial position for dummy masses in real (non-local) coordinates.
     * Cheat by using the routine to calculate virtual site parameters. It is
     * much easier when we have the coordinates expressed in terms of
     * CB, CG, CD2.
     */
    rvec_sub(x[ats[atCB]], x[ats[atCG]], r_ij);
    rvec_sub(x[ats[atCD2]], x[ats[atCG]], r_ik);
    calc_vsite3_param(xcom[0], ycom[0], xi[atCG], yi[atCG], xi[atCB], yi[atCB],
                      xi[atCD2], yi[atCD2], &a, &b);
    svmul(a, r_ij, t1);
    svmul(b, r_ik, t2);
    rvec_add(t1, t2, t1);
    rvec_add(t1, x[ats[atCG]], (*newx)[atM[0]]);

    calc_vsite3_param(xcom[1], ycom[1], xi[atCG], yi[atCG], xi[atCB], yi[atCB],
                      xi[atCD2], yi[atCD2], &a, &b);
    svmul(a, r_ij, t1);
    svmul(b, r_ik, t2);
    rvec_add(t1, t2, t1);
    rvec_add(t1, x[ats[atCG]], (*newx)[atM[1]]);

    /* set parameters for the masses */
    for (j = 0; j < NMASS; j++)
    {
        sprintf(name, "MW%d", j+1);
        (*newatomname)  [atM[j]]         = put_symtab(symtab, name);
        (*newatom)      [atM[j]].m       = (*newatom)[atM[j]].mB    = mM[j];
        (*newatom)      [atM[j]].q       = (*newatom)[atM[j]].qB    = 0.0;
        (*newatom)      [atM[j]].type    = (*newatom)[atM[j]].typeB = tpM;
        (*newatom)      [atM[j]].ptype   = eptAtom;
        (*newatom)      [atM[j]].resind  = at->atom[i0].resind;
        (*newatom)      [atM[j]].elem[0] = 'M';
        (*newatom)      [atM[j]].elem[1] = '\0';
        (*newvsite_type)[atM[j]]         = NOTSET;
        (*newcgnr)      [atM[j]]         = (*cgnr)[i0];
    }
    /* renumber cgnr: */
    for (i = i0; i < at->nr; i++)
    {
        (*cgnr)[i]++;
    }

    /* constraints between CB, M1 and M2 */
    /* 'add_shift' says which atoms won't be renumbered afterwards */
    dCBM1 = std::hypot( xcom[0]-xi[atCB], ycom[0]-yi[atCB] );
    dM1M2 = std::hypot( xcom[0]-xcom[1], ycom[0]-ycom[1] );
    dCBM2 = std::hypot( xcom[1]-xi[atCB], ycom[1]-yi[atCB] );
    my_add_param(&(plist[F_CONSTRNC]), ats[atCB],       add_shift+atM[0], dCBM1);
    my_add_param(&(plist[F_CONSTRNC]), ats[atCB],       add_shift+atM[1], dCBM2);
    my_add_param(&(plist[F_CONSTRNC]), add_shift+atM[0], add_shift+atM[1], dM1M2);

    /* rest will be vsite3 */
    nvsite = 0;
    for (i = 0; i < atNR; i++)
    {
        if (i != atCB)
        {
            at->atom[ats[i]].m    = at->atom[ats[i]].mB = 0;
            (*vsite_type)[ats[i]] = F_VSITE3;
            nvsite++;
        }
    }

    /* now define all vsites from M1, M2, CB, ie:
       r_d = r_M1 + a r_M1_M2 + b r_M1_CB */
    for (i = 0; i < atNR; i++)
    {
        if ( (*vsite_type)[ats[i]] == F_VSITE3)
        {
            calc_vsite3_param(xi[i], yi[i], xcom[0], ycom[0], xcom[1], ycom[1], xi[atCB], yi[atCB], &a, &b);
            add_vsite3_param(&plist[F_VSITE3],
                             ats[i], add_shift+atM[0], add_shift+atM[1], ats[atCB], a, b);
        }
    }
    return nvsite;
#undef NMASS
}


static int gen_vsites_tyr(gpp_atomtype_t atype, rvec *newx[],
                          t_atom *newatom[], char ***newatomname[],
                          int *o2n[], int *newvsite_type[], int *newcgnr[],
                          t_symtab *symtab, int *nadd, rvec x[], int *cgnr[],
                          t_atoms *at, int *vsite_type[], t_params plist[],
                          int nrfound, int *ats, int add_shift,
                          t_vsitetop *vsitetop, int nvsitetop)
{
    int  nvsite, i, i0, j, atM, tpM;
    real dCGCE, dCEOH, dCGM, tmp1, a, b;
    real bond_cc, bond_ch, bond_co, bond_oh, angle_coh;
    real xcom, mtot;
    real vdist, mM;
    rvec r1;
    char name[10];

    /* these MUST correspond to the atnms array in do_vsite_aromatics! */
    enum {
        atCG, atCD1, atHD1, atCD2, atHD2, atCE1, atHE1, atCE2, atHE2,
        atCZ, atOH, atHH, atNR
    };
    real xi[atNR];
    /* CG, CE1, CE2 (as in general 6-ring) and OH and HH stay,
       rest gets virtualized.
       Now we have two linked triangles with one improper keeping them flat */
    if (atNR != nrfound)
    {
        gmx_incons("Number of atom types in gen_vsites_tyr");
    }

    /* Aromatic rings have 6-fold symmetry, so we only need one bond length
     * for the ring part (angle is always 120 degrees).
     */
    bond_cc   = get_ddb_bond(vsitetop, nvsitetop, "TYR", "CD1", "CE1");
    bond_ch   = get_ddb_bond(vsitetop, nvsitetop, "TYR", "CD1", "HD1");
    bond_co   = get_ddb_bond(vsitetop, nvsitetop, "TYR", "CZ", "OH");
    bond_oh   = get_ddb_bond(vsitetop, nvsitetop, "TYR", "OH", "HH");
    angle_coh = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, "TYR", "CZ", "OH", "HH");

    xi[atCG]  = -bond_cc+bond_cc*std::cos(ANGLE_6RING);
    xi[atCD1] = -bond_cc;
    xi[atHD1] = xi[atCD1]+bond_ch*std::cos(ANGLE_6RING);
    xi[atCE1] = 0;
    xi[atHE1] = xi[atCE1]-bond_ch*std::cos(ANGLE_6RING);
    xi[atCD2] = xi[atCD1];
    xi[atHD2] = xi[atHD1];
    xi[atCE2] = xi[atCE1];
    xi[atHE2] = xi[atHE1];
    xi[atCZ]  = bond_cc*std::cos(0.5*ANGLE_6RING);
    xi[atOH]  = xi[atCZ]+bond_co;

    xcom = mtot = 0;
    for (i = 0; i < atOH; i++)
    {
        xcom += xi[i]*at->atom[ats[i]].m;
        mtot += at->atom[ats[i]].m;
    }
    xcom /= mtot;

    /* first do 6 ring as default,
       except CZ (we'll do that different) and HZ (we don't have that): */
    nvsite = gen_vsites_6ring(at, vsite_type, plist, nrfound, ats, bond_cc, bond_ch, xcom, FALSE);

    /* then construct CZ from the 2nd triangle */
    /* vsite3 construction: r_d = r_i + a r_ij + b r_ik */
    a = b = 0.5 * bond_co / ( bond_co - bond_cc*std::cos(ANGLE_6RING) );
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atCZ], ats[atOH], ats[atCE1], ats[atCE2], a, b);
    at->atom[ats[atCZ]].m = at->atom[ats[atCZ]].mB = 0;

    /* constraints between CE1, CE2 and OH */
    dCGCE = std::sqrt( cosrule(bond_cc, bond_cc, ANGLE_6RING) );
    dCEOH = std::sqrt( cosrule(bond_cc, bond_co, ANGLE_6RING) );
    my_add_param(&(plist[F_CONSTRNC]), ats[atCE1], ats[atOH], dCEOH);
    my_add_param(&(plist[F_CONSTRNC]), ats[atCE2], ats[atOH], dCEOH);

    /* We also want to constrain the angle C-O-H, but since CZ is constructed
     * we need to introduce a constraint to CG.
     * CG is much further away, so that will lead to instabilities in LINCS
     * when we constrain both CG-HH and OH-HH distances. Instead of requiring
     * the use of lincs_order=8 we introduce a dummy mass three times further
     * away from OH than HH. The mass is accordingly a third, with the remaining
     * 2/3 moved to OH. This shouldn't cause any problems since the forces will
     * apply to the HH constructed atom and not directly on the virtual mass.
     */

    vdist                   = 2.0*bond_oh;
    mM                      = at->atom[ats[atHH]].m/2.0;
    at->atom[ats[atOH]].m  += mM; /* add 1/2 of original H mass */
    at->atom[ats[atOH]].mB += mM; /* add 1/2 of original H mass */
    at->atom[ats[atHH]].m   = at->atom[ats[atHH]].mB = 0;

    /* get dummy mass type */
    tpM = vsite_nm2type("MW", atype);
    /* make space for 1 mass: shift HH only */
    i0  = ats[atHH];
    atM = i0+*nadd;
    if (debug)
    {
        fprintf(stderr, "Inserting 1 dummy mass at %d\n", (*o2n)[i0]+1);
    }
    (*nadd)++;
    for (j = i0; j < at->nr; j++)
    {
        (*o2n)[j] = j+*nadd;
    }
    srenew(*newx, at->nr+*nadd);
    srenew(*newatom, at->nr+*nadd);
    srenew(*newatomname, at->nr+*nadd);
    srenew(*newvsite_type, at->nr+*nadd);
    srenew(*newcgnr, at->nr+*nadd);
    (*newatomname)[at->nr+*nadd-1] = nullptr;

    /* Calc the dummy mass initial position */
    rvec_sub(x[ats[atHH]], x[ats[atOH]], r1);
    svmul(2.0, r1, r1);
    rvec_add(r1, x[ats[atHH]], (*newx)[atM]);

    strcpy(name, "MW1");
    (*newatomname)  [atM]         = put_symtab(symtab, name);
    (*newatom)      [atM].m       = (*newatom)[atM].mB    = mM;
    (*newatom)      [atM].q       = (*newatom)[atM].qB    = 0.0;
    (*newatom)      [atM].type    = (*newatom)[atM].typeB = tpM;
    (*newatom)      [atM].ptype   = eptAtom;
    (*newatom)      [atM].resind  = at->atom[i0].resind;
    (*newatom)      [atM].elem[0] = 'M';
    (*newatom)      [atM].elem[1] = '\0';
    (*newvsite_type)[atM]         = NOTSET;
    (*newcgnr)      [atM]         = (*cgnr)[i0];
    /* renumber cgnr: */
    for (i = i0; i < at->nr; i++)
    {
        (*cgnr)[i]++;
    }

    (*vsite_type)[ats[atHH]] = F_VSITE2;
    nvsite++;
    /* assume we also want the COH angle constrained: */
    tmp1 = bond_cc*std::cos(0.5*ANGLE_6RING) + dCGCE*std::sin(ANGLE_6RING*0.5) + bond_co;
    dCGM = std::sqrt( cosrule(tmp1, vdist, angle_coh) );
    my_add_param(&(plist[F_CONSTRNC]), ats[atCG], add_shift+atM, dCGM);
    my_add_param(&(plist[F_CONSTRNC]), ats[atOH], add_shift+atM, vdist);

    add_vsite2_param(&plist[F_VSITE2],
                     ats[atHH], ats[atOH], add_shift+atM, 1.0/2.0);
    return nvsite;
}

static int gen_vsites_his(t_atoms *at, int *vsite_type[], t_params plist[],
                          int nrfound, int *ats, t_vsitetop *vsitetop, int nvsitetop)
{
    int  nvsite, i;
    real a, b, alpha, dCGCE1, dCGNE2;
    real sinalpha, cosalpha;
    real xcom, ycom, mtot;
    real mG, mrest, mCE1, mNE2;
    real b_CG_ND1, b_ND1_CE1, b_CE1_NE2, b_CG_CD2, b_CD2_NE2;
    real b_ND1_HD1, b_NE2_HE2, b_CE1_HE1, b_CD2_HD2;
    real a_CG_ND1_CE1, a_CG_CD2_NE2, a_ND1_CE1_NE2, a_CE1_NE2_CD2;
    real a_NE2_CE1_HE1, a_NE2_CD2_HD2, a_CE1_ND1_HD1, a_CE1_NE2_HE2;
    char resname[10];

    /* these MUST correspond to the atnms array in do_vsite_aromatics! */
    enum {
        atCG, atND1, atHD1, atCD2, atHD2, atCE1, atHE1, atNE2, atHE2, atNR
    };
    real x[atNR], y[atNR];

    /* CG, CE1 and NE2 stay, each gets part of the total mass,
       rest gets virtualized */
    /* check number of atoms, 3 hydrogens may be missing: */
    /* assert( nrfound >= atNR-3 || nrfound <= atNR );
     * Don't understand the above logic. Shouldn't it be && rather than || ???
     */
    if ((nrfound < atNR-3) || (nrfound > atNR))
    {
        gmx_incons("Generating vsites for HIS");
    }

    /* avoid warnings about uninitialized variables */
    b_ND1_HD1 = b_NE2_HE2 = b_CE1_HE1 = b_CD2_HD2 = a_NE2_CE1_HE1 =
                        a_NE2_CD2_HD2 = a_CE1_ND1_HD1 = a_CE1_NE2_HE2 = 0;

    if (ats[atHD1] != NOTSET)
    {
        if (ats[atHE2] != NOTSET)
        {
            sprintf(resname, "HISH");
        }
        else
        {
            sprintf(resname, "HISA");
        }
    }
    else
    {
        sprintf(resname, "HISB");
    }

    /* Get geometry from database */
    b_CG_ND1      = get_ddb_bond(vsitetop, nvsitetop, resname, "CG", "ND1");
    b_ND1_CE1     = get_ddb_bond(vsitetop, nvsitetop, resname, "ND1", "CE1");
    b_CE1_NE2     = get_ddb_bond(vsitetop, nvsitetop, resname, "CE1", "NE2");
    b_CG_CD2      = get_ddb_bond(vsitetop, nvsitetop, resname, "CG", "CD2");
    b_CD2_NE2     = get_ddb_bond(vsitetop, nvsitetop, resname, "CD2", "NE2");
    a_CG_ND1_CE1  = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "CG", "ND1", "CE1");
    a_CG_CD2_NE2  = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "CG", "CD2", "NE2");
    a_ND1_CE1_NE2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "ND1", "CE1", "NE2");
    a_CE1_NE2_CD2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "CE1", "NE2", "CD2");

    if (ats[atHD1] != NOTSET)
    {
        b_ND1_HD1     = get_ddb_bond(vsitetop, nvsitetop, resname, "ND1", "HD1");
        a_CE1_ND1_HD1 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "CE1", "ND1", "HD1");
    }
    if (ats[atHE2] != NOTSET)
    {
        b_NE2_HE2     = get_ddb_bond(vsitetop, nvsitetop, resname, "NE2", "HE2");
        a_CE1_NE2_HE2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "CE1", "NE2", "HE2");
    }
    if (ats[atHD2] != NOTSET)
    {
        b_CD2_HD2     = get_ddb_bond(vsitetop, nvsitetop, resname, "CD2", "HD2");
        a_NE2_CD2_HD2 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "NE2", "CD2", "HD2");
    }
    if (ats[atHE1] != NOTSET)
    {
        b_CE1_HE1     = get_ddb_bond(vsitetop, nvsitetop, resname, "CE1", "HE1");
        a_NE2_CE1_HE1 = DEG2RAD*get_ddb_angle(vsitetop, nvsitetop, resname, "NE2", "CE1", "HE1");
    }

    /* constraints between CG, CE1 and NE1 */
    dCGCE1   = std::sqrt( cosrule(b_CG_ND1, b_ND1_CE1, a_CG_ND1_CE1) );
    dCGNE2   = std::sqrt( cosrule(b_CG_CD2, b_CD2_NE2, a_CG_CD2_NE2) );

    my_add_param(&(plist[F_CONSTRNC]), ats[atCG], ats[atCE1], dCGCE1);
    my_add_param(&(plist[F_CONSTRNC]), ats[atCG], ats[atNE2], dCGNE2);
    /* we already have a constraint CE1-NE2, so we don't add it again */

    /* calculate the positions in a local frame of reference.
     * The x-axis is the line from CG that makes a right angle
     * with the bond CE1-NE2, and the y-axis the bond CE1-NE2.
     */
    /* First calculate the x-axis intersection with y-axis (=yCE1).
     * Get cos(angle CG-CE1-NE2) :
     */
    cosalpha = acosrule(dCGNE2, dCGCE1, b_CE1_NE2);
    x[atCE1] = 0;
    y[atCE1] = cosalpha*dCGCE1;
    x[atNE2] = 0;
    y[atNE2] = y[atCE1]-b_CE1_NE2;
    sinalpha = std::sqrt(1-cosalpha*cosalpha);
    x[atCG]  = -sinalpha*dCGCE1;
    y[atCG]  = 0;
    x[atHE1] = x[atHE2] = x[atHD1] = x[atHD2] = 0;
    y[atHE1] = y[atHE2] = y[atHD1] = y[atHD2] = 0;

    /* calculate ND1 and CD2 positions from CE1 and NE2 */

    x[atND1] = -b_ND1_CE1*std::sin(a_ND1_CE1_NE2);
    y[atND1] = y[atCE1]-b_ND1_CE1*std::cos(a_ND1_CE1_NE2);

    x[atCD2] = -b_CD2_NE2*std::sin(a_CE1_NE2_CD2);
    y[atCD2] = y[atNE2]+b_CD2_NE2*std::cos(a_CE1_NE2_CD2);

    /* And finally the hydrogen positions */
    if (ats[atHE1] != NOTSET)
    {
        x[atHE1] = x[atCE1] + b_CE1_HE1*std::sin(a_NE2_CE1_HE1);
        y[atHE1] = y[atCE1] - b_CE1_HE1*std::cos(a_NE2_CE1_HE1);
    }
    /* HD2 - first get (ccw) angle from (positive) y-axis */
    if (ats[atHD2] != NOTSET)
    {
        alpha    = a_CE1_NE2_CD2 + M_PI - a_NE2_CD2_HD2;
        x[atHD2] = x[atCD2] - b_CD2_HD2*std::sin(alpha);
        y[atHD2] = y[atCD2] + b_CD2_HD2*std::cos(alpha);
    }
    if (ats[atHD1] != NOTSET)
    {
        /* HD1 - first get (cw) angle from (positive) y-axis */
        alpha    = a_ND1_CE1_NE2 + M_PI - a_CE1_ND1_HD1;
        x[atHD1] = x[atND1] - b_ND1_HD1*std::sin(alpha);
        y[atHD1] = y[atND1] - b_ND1_HD1*std::cos(alpha);
    }
    if (ats[atHE2] != NOTSET)
    {
        x[atHE2] = x[atNE2] + b_NE2_HE2*std::sin(a_CE1_NE2_HE2);
        y[atHE2] = y[atNE2] + b_NE2_HE2*std::cos(a_CE1_NE2_HE2);
    }
    /* Have all coordinates now */

    /* calc center-of-mass; keep atoms CG, CE1, NE2 and
     * set the rest to vsite3
     */
    mtot   = xcom = ycom = 0;
    nvsite = 0;
    for (i = 0; i < atNR; i++)
    {
        if (ats[i] != NOTSET)
        {
            mtot += at->atom[ats[i]].m;
            xcom += x[i]*at->atom[ats[i]].m;
            ycom += y[i]*at->atom[ats[i]].m;
            if (i != atCG && i != atCE1 && i != atNE2)
            {
                at->atom[ats[i]].m    = at->atom[ats[i]].mB = 0;
                (*vsite_type)[ats[i]] = F_VSITE3;
                nvsite++;
            }
        }
    }
    if (nvsite+3 != nrfound)
    {
        gmx_incons("Generating vsites for HIS");
    }

    xcom /= mtot;
    ycom /= mtot;

    /* distribute mass so that com stays the same */
    mG    = xcom*mtot/x[atCG];
    mrest = mtot-mG;
    mCE1  = (ycom-y[atNE2])*mrest/(y[atCE1]-y[atNE2]);
    mNE2  = mrest-mCE1;

    at->atom[ats[atCG]].m  = at->atom[ats[atCG]].mB = mG;
    at->atom[ats[atCE1]].m = at->atom[ats[atCE1]].mB = mCE1;
    at->atom[ats[atNE2]].m = at->atom[ats[atNE2]].mB = mNE2;

    /* HE1 */
    if (ats[atHE1] != NOTSET)
    {
        calc_vsite3_param(x[atHE1], y[atHE1], x[atCE1], y[atCE1], x[atNE2], y[atNE2],
                          x[atCG], y[atCG], &a, &b);
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atHE1], ats[atCE1], ats[atNE2], ats[atCG], a, b);
    }
    /* HE2 */
    if (ats[atHE2] != NOTSET)
    {
        calc_vsite3_param(x[atHE2], y[atHE2], x[atNE2], y[atNE2], x[atCE1], y[atCE1],
                          x[atCG], y[atCG], &a, &b);
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atHE2], ats[atNE2], ats[atCE1], ats[atCG], a, b);
    }

    /* ND1 */
    calc_vsite3_param(x[atND1], y[atND1], x[atNE2], y[atNE2], x[atCE1], y[atCE1],
                      x[atCG], y[atCG], &a, &b);
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atND1], ats[atNE2], ats[atCE1], ats[atCG], a, b);

    /* CD2 */
    calc_vsite3_param(x[atCD2], y[atCD2], x[atCE1], y[atCE1], x[atNE2], y[atNE2],
                      x[atCG], y[atCG], &a, &b);
    add_vsite3_param(&plist[F_VSITE3],
                     ats[atCD2], ats[atCE1], ats[atNE2], ats[atCG], a, b);

    /* HD1 */
    if (ats[atHD1] != NOTSET)
    {
        calc_vsite3_param(x[atHD1], y[atHD1], x[atNE2], y[atNE2], x[atCE1], y[atCE1],
                          x[atCG], y[atCG], &a, &b);
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atHD1], ats[atNE2], ats[atCE1], ats[atCG], a, b);
    }
    /* HD2 */
    if (ats[atHD2] != NOTSET)
    {
        calc_vsite3_param(x[atHD2], y[atHD2], x[atCE1], y[atCE1], x[atNE2], y[atNE2],
                          x[atCG], y[atCG], &a, &b);
        add_vsite3_param(&plist[F_VSITE3],
                         ats[atHD2], ats[atCE1], ats[atNE2], ats[atCG], a, b);
    }
    return nvsite;
}

static gmx_bool is_vsite(int vsite_type)
{
    if (vsite_type == NOTSET)
    {
        return FALSE;
    }
    switch (abs(vsite_type) )
    {
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3OUT:
        case F_VSITE3FAD:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            return TRUE;
        default:
            return FALSE;
    }
}

static char atomnamesuffix[] = "1234";

void do_vsites(int nrtp, t_restp rtp[], gpp_atomtype_t atype,
               t_atoms *at, t_symtab *symtab, rvec *x[],
               t_params plist[], int *vsite_type[], int *cgnr[],
               real mHmult, gmx_bool bVsiteAromatics,
               const char *ffdir)
{
#define MAXATOMSPERRESIDUE 16
    int               i, j, k, m, i0, ni0, whatres, resind, add_shift, ftype, nvsite, nadd;
    int               ai, aj, ak, al;
    int               nrfound = 0, needed, nrbonds, nrHatoms, Heavy, nrheavies, tpM, tpHeavy;
    int               Hatoms[4], heavies[4];
    gmx_bool          bWARNING, bAddVsiteParam, bFirstWater;
    matrix            tmpmat;
    gmx_bool         *bResProcessed;
    real              mHtot, mtot, fact, fact2;
    rvec              rpar, rperp, temp;
    char              name[10], tpname[32], nexttpname[32], *ch;
    rvec             *newx;
    int              *o2n, *newvsite_type, *newcgnr, ats[MAXATOMSPERRESIDUE];
    t_atom           *newatom;
    t_params         *params;
    char           ***newatomname;
    char             *resnm = nullptr;
    int               ndb, f;
    char            **db;
    int               nvsiteconf, nvsitetop, cmplength;
    gmx_bool          isN, planarN, bFound;
    gmx_residuetype_t*rt;

    t_vsiteconf      *vsiteconflist;
    /* pointer to a list of CH3/NH3/NH2 configuration entries.
     * See comments in read_vsite_database. It isnt beautiful,
     * but it had to be fixed, and I dont even want to try to
     * maintain this part of the code...
     */
    t_vsitetop *vsitetop;
    /* Pointer to a list of geometry (bond/angle) entries for
     * residues like PHE, TRP, TYR, HIS, etc., where we need
     * to know the geometry to construct vsite aromatics.
     * Note that equilibrium geometry isnt necessarily the same
     * as the individual bond and angle values given in the
     * force field (rings can be strained).
     */

    /* if bVsiteAromatics=TRUE do_vsites will specifically convert atoms in
       PHE, TRP, TYR and HIS to a construction of virtual sites */
    enum                    {
        resPHE, resTRP, resTYR, resHIS, resNR
    };
    const char *resnms[resNR]   = {   "PHE",  "TRP",  "TYR",  "HIS" };
    /* Amber03 alternative names for termini */
    const char *resnmsN[resNR]  = {  "NPHE", "NTRP", "NTYR", "NHIS" };
    const char *resnmsC[resNR]  = {  "CPHE", "CTRP", "CTYR", "CHIS" };
    /* HIS can be known as HISH, HIS1, HISA, HID, HIE, HIP, etc. too */
    gmx_bool    bPartial[resNR]  = {  FALSE,  FALSE,  FALSE,   TRUE  };
    /* the atnms for every residue MUST correspond to the enums in the
       gen_vsites_* (one for each residue) routines! */
    /* also the atom names in atnms MUST be in the same order as in the .rtp! */
    const char *atnms[resNR][MAXATOMSPERRESIDUE+1] = {
        { "CG", /* PHE */
          "CD1", "HD1", "CD2", "HD2",
          "CE1", "HE1", "CE2", "HE2",
          "CZ", "HZ", nullptr },
        { "CB", /* TRP */
          "CG",
          "CD1", "HD1", "CD2",
          "NE1", "HE1", "CE2", "CE3", "HE3",
          "CZ2", "HZ2", "CZ3", "HZ3",
          "CH2", "HH2", nullptr },
        { "CG", /* TYR */
          "CD1", "HD1", "CD2", "HD2",
          "CE1", "HE1", "CE2", "HE2",
          "CZ", "OH", "HH", nullptr },
        { "CG", /* HIS */
          "ND1", "HD1", "CD2", "HD2",
          "CE1", "HE1", "NE2", "HE2", nullptr }
    };

    if (debug)
    {
        printf("Searching for atoms to make virtual sites ...\n");
        fprintf(debug, "# # # VSITES # # #\n");
    }

    ndb           = fflib_search_file_end(ffdir, ".vsd", FALSE, &db);
    nvsiteconf    = 0;
    vsiteconflist = nullptr;
    nvsitetop     = 0;
    vsitetop      = nullptr;
    for (f = 0; f < ndb; f++)
    {
        read_vsite_database(db[f], &vsiteconflist, &nvsiteconf, &vsitetop, &nvsitetop);
        sfree(db[f]);
    }
    sfree(db);

    bFirstWater = TRUE;
    nvsite      = 0;
    nadd        = 0;
    /* we need a marker for which atoms should *not* be renumbered afterwards */
    add_shift = 10*at->nr;
    /* make arrays where masses can be inserted into */
    snew(newx, at->nr);
    snew(newatom, at->nr);
    snew(newatomname, at->nr);
    snew(newvsite_type, at->nr);
    snew(newcgnr, at->nr);
    /* make index array to tell where the atoms go to when masses are inserted */
    snew(o2n, at->nr);
    for (i = 0; i < at->nr; i++)
    {
        o2n[i] = i;
    }
    /* make index to tell which residues were already processed */
    snew(bResProcessed, at->nres);

    gmx_residuetype_init(&rt);

    /* generate vsite constructions */
    /* loop over all atoms */
    resind = -1;
    for (i = 0; (i < at->nr); i++)
    {
        if (at->atom[i].resind != resind)
        {
            resind = at->atom[i].resind;
            resnm  = *(at->resinfo[resind].name);
        }
        /* first check for aromatics to virtualize */
        /* don't waste our effort on DNA, water etc. */
        /* Only do the vsite aromatic stuff when we reach the
         * CA atom, since there might be an X2/X3 group on the
         * N-terminus that must be treated first.
         */
        if (bVsiteAromatics &&
            !strcmp(*(at->atomname[i]), "CA") &&
            !bResProcessed[resind] &&
            gmx_residuetype_is_protein(rt, *(at->resinfo[resind].name)) )
        {
            /* mark this residue */
            bResProcessed[resind] = TRUE;
            /* find out if this residue needs converting */
            whatres = NOTSET;
            for (j = 0; j < resNR && whatres == NOTSET; j++)
            {

                cmplength = bPartial[j] ? strlen(resnm)-1 : strlen(resnm);

                bFound = ((gmx_strncasecmp(resnm, resnms[j], cmplength) == 0) ||
                          (gmx_strncasecmp(resnm, resnmsN[j], cmplength) == 0) ||
                          (gmx_strncasecmp(resnm, resnmsC[j], cmplength) == 0));

                if (bFound)
                {
                    whatres = j;
                    /* get atoms we will be needing for the conversion */
                    nrfound = 0;
                    for (k = 0; atnms[j][k]; k++)
                    {
                        ats[k] = NOTSET;
                        for (m = i; m < at->nr && at->atom[m].resind == resind && ats[k] == NOTSET; m++)
                        {
                            if (gmx_strcasecmp(*(at->atomname[m]), atnms[j][k]) == 0)
                            {
                                ats[k] = m;
                                nrfound++;
                            }
                        }
                    }

                    /* now k is number of atom names in atnms[j] */
                    if (j == resHIS)
                    {
                        needed = k-3;
                    }
                    else
                    {
                        needed = k;
                    }
                    if (nrfound < needed)
                    {
                        gmx_fatal(FARGS, "not enough atoms found (%d, need %d) in "
                                  "residue %s %d while\n             "
                                  "generating aromatics virtual site construction",
                                  nrfound, needed, resnm, at->resinfo[resind].nr);
                    }
                    /* Advance overall atom counter */
                    i++;
                }
            }
            /* the enums for every residue MUST correspond to atnms[residue] */
            switch (whatres)
            {
                case resPHE:
                    if (debug)
                    {
                        fprintf(stderr, "PHE at %d\n", o2n[ats[0]]+1);
                    }
                    nvsite += gen_vsites_phe(at, vsite_type, plist, nrfound, ats, vsitetop, nvsitetop);
                    break;
                case resTRP:
                    if (debug)
                    {
                        fprintf(stderr, "TRP at %d\n", o2n[ats[0]]+1);
                    }
                    nvsite += gen_vsites_trp(atype, &newx, &newatom, &newatomname, &o2n,
                                             &newvsite_type, &newcgnr, symtab, &nadd, *x, cgnr,
                                             at, vsite_type, plist, nrfound, ats, add_shift, vsitetop, nvsitetop);
                    break;
                case resTYR:
                    if (debug)
                    {
                        fprintf(stderr, "TYR at %d\n", o2n[ats[0]]+1);
                    }
                    nvsite += gen_vsites_tyr(atype, &newx, &newatom, &newatomname, &o2n,
                                             &newvsite_type, &newcgnr, symtab, &nadd, *x, cgnr,
                                             at, vsite_type, plist, nrfound, ats, add_shift, vsitetop, nvsitetop);
                    break;
                case resHIS:
                    if (debug)
                    {
                        fprintf(stderr, "HIS at %d\n", o2n[ats[0]]+1);
                    }
                    nvsite += gen_vsites_his(at, vsite_type, plist, nrfound, ats, vsitetop, nvsitetop);
                    break;
                case NOTSET:
                    /* this means this residue won't be processed */
                    break;
                default:
                    gmx_fatal(FARGS, "DEATH HORROR in do_vsites (%s:%d)",
                              __FILE__, __LINE__);
            } /* switch whatres */
              /* skip back to beginning of residue */
            while (i > 0 && at->atom[i-1].resind == resind)
            {
                i--;
            }
        } /* if bVsiteAromatics & is protein */

        /* now process the rest of the hydrogens */
        /* only process hydrogen atoms which are not already set */
        if ( ((*vsite_type)[i] == NOTSET) && is_hydrogen(*(at->atomname[i])))
        {
            /* find heavy atom, count #bonds from it and #H atoms bound to it
               and return H atom numbers (Hatoms) and heavy atom numbers (heavies) */
            count_bonds(i, &plist[F_BONDS], at->atomname,
                        &nrbonds, &nrHatoms, Hatoms, &Heavy, &nrheavies, heavies);
            /* get Heavy atom type */
            tpHeavy = get_atype(Heavy, at, nrtp, rtp, rt);
            strcpy(tpname, get_atomtype_name(tpHeavy, atype));

            bWARNING       = FALSE;
            bAddVsiteParam = TRUE;
            /* nested if's which check nrHatoms, nrbonds and atomname */
            if (nrHatoms == 1)
            {
                switch (nrbonds)
                {
                    case 2: /* -O-H */
                        (*vsite_type)[i] = F_BONDS;
                        break;
                    case 3: /* =CH-, -NH- or =NH+- */
                        (*vsite_type)[i] = F_VSITE3FD;
                        break;
                    case 4: /* --CH- (tert) */
                        /* The old type 4FD had stability issues, so
                         * all new constructs should use 4FDN
                         */
                        (*vsite_type)[i] = F_VSITE4FDN;

                        /* Check parity of heavy atoms from coordinates */
                        ai = Heavy;
                        aj = heavies[0];
                        ak = heavies[1];
                        al = heavies[2];
                        rvec_sub((*x)[aj], (*x)[ai], tmpmat[0]);
                        rvec_sub((*x)[ak], (*x)[ai], tmpmat[1]);
                        rvec_sub((*x)[al], (*x)[ai], tmpmat[2]);

                        if (det(tmpmat) > 0)
                        {
                            /* swap parity */
                            heavies[1] = aj;
                            heavies[0] = ak;
                        }

                        break;
                    default: /* nrbonds != 2, 3 or 4 */
                        bWARNING = TRUE;
                }

            }
            else if ( (nrHatoms == 2) && (nrbonds == 2) &&
                      (at->atom[Heavy].atomnumber == 8) )
            {
                bAddVsiteParam = FALSE; /* this is water: skip these hydrogens */
                if (bFirstWater)
                {
                    bFirstWater = FALSE;
                    if (debug)
                    {
                        fprintf(debug,
                                "Not converting hydrogens in water to virtual sites\n");
                    }
                }
            }
            else if ( (nrHatoms == 2) && (nrbonds == 4) )
            {
                /* -CH2- , -NH2+- */
                (*vsite_type)[Hatoms[0]] = F_VSITE3OUT;
                (*vsite_type)[Hatoms[1]] = -F_VSITE3OUT;
            }
            else
            {
                /* 2 or 3 hydrogen atom, with 3 or 4 bonds in total to the heavy atom.
                 * If it is a nitrogen, first check if it is planar.
                 */
                isN = planarN = FALSE;
                if ((nrHatoms == 2) && ((*at->atomname[Heavy])[0] == 'N'))
                {
                    isN = TRUE;
                    j   = nitrogen_is_planar(vsiteconflist, nvsiteconf, tpname);
                    if (j < 0)
                    {
                        gmx_fatal(FARGS, "No vsite database NH2 entry for type %s\n", tpname);
                    }
                    planarN = (j == 1);
                }
                if ( (nrHatoms == 2) && (nrbonds == 3) && ( !isN || planarN ) )
                {
                    /* =CH2 or, if it is a nitrogen NH2, it is a planar one */
                    (*vsite_type)[Hatoms[0]] = F_VSITE3FAD;
                    (*vsite_type)[Hatoms[1]] = -F_VSITE3FAD;
                }
                else if ( ( (nrHatoms == 2) && (nrbonds == 3) &&
                            ( isN && !planarN ) ) ||
                          ( (nrHatoms == 3) && (nrbonds == 4) ) )
                {
                    /* CH3, NH3 or non-planar NH2 group */
                    int      Hat_vsite_type[3] = { F_VSITE3, F_VSITE3OUT, F_VSITE3OUT };
                    gmx_bool Hat_SwapParity[3] = { FALSE,    TRUE,        FALSE };

                    if (debug)
                    {
                        fprintf(stderr, "-XH3 or nonplanar NH2 group at %d\n", i+1);
                    }
                    bAddVsiteParam = FALSE; /* we'll do this ourselves! */
                    /* -NH2 (umbrella), -NH3+ or -CH3 */
                    (*vsite_type)[Heavy]       = F_VSITE3;
                    for (j = 0; j < nrHatoms; j++)
                    {
                        (*vsite_type)[Hatoms[j]] = Hat_vsite_type[j];
                    }
                    /* get dummy mass type from first char of heavy atom type (N or C) */

                    strcpy(nexttpname, get_atomtype_name(get_atype(heavies[0], at, nrtp, rtp, rt), atype));
                    ch = get_dummymass_name(vsiteconflist, nvsiteconf, tpname, nexttpname);

                    if (ch == nullptr)
                    {
                        if (ndb > 0)
                        {
                            gmx_fatal(FARGS, "Can't find dummy mass for type %s bonded to type %s in the virtual site database (.vsd files). Add it to the database!\n", tpname, nexttpname);
                        }
                        else
                        {
                            gmx_fatal(FARGS, "A dummy mass for type %s bonded to type %s is required, but no virtual site database (.vsd) files where found.\n", tpname, nexttpname);
                        }
                    }
                    else
                    {
                        strcpy(name, ch);
                    }

                    tpM = vsite_nm2type(name, atype);
                    /* make space for 2 masses: shift all atoms starting with 'Heavy' */
#define NMASS 2
                    i0  = Heavy;
                    ni0 = i0+nadd;
                    if (debug)
                    {
                        fprintf(stderr, "Inserting %d dummy masses at %d\n", NMASS, o2n[i0]+1);
                    }
                    nadd += NMASS;
                    for (j = i0; j < at->nr; j++)
                    {
                        o2n[j] = j+nadd;
                    }

                    srenew(newx, at->nr+nadd);
                    srenew(newatom, at->nr+nadd);
                    srenew(newatomname, at->nr+nadd);
                    srenew(newvsite_type, at->nr+nadd);
                    srenew(newcgnr, at->nr+nadd);

                    for (j = 0; j < NMASS; j++)
                    {
                        newatomname[at->nr+nadd-1-j] = nullptr;
                    }

                    /* calculate starting position for the masses */
                    mHtot = 0;
                    /* get atom masses, and set Heavy and Hatoms mass to zero */
                    for (j = 0; j < nrHatoms; j++)
                    {
                        mHtot                += get_amass(Hatoms[j], at, nrtp, rtp, rt);
                        at->atom[Hatoms[j]].m = at->atom[Hatoms[j]].mB = 0;
                    }
                    mtot              = mHtot + get_amass(Heavy, at, nrtp, rtp, rt);
                    at->atom[Heavy].m = at->atom[Heavy].mB = 0;
                    if (mHmult != 1.0)
                    {
                        mHtot *= mHmult;
                    }
                    fact2 = mHtot/mtot;
                    fact  = std::sqrt(fact2);
                    /* generate vectors parallel and perpendicular to rotational axis:
                     * rpar  = Heavy -> Hcom
                     * rperp = Hcom  -> H1   */
                    clear_rvec(rpar);
                    for (j = 0; j < nrHatoms; j++)
                    {
                        rvec_inc(rpar, (*x)[Hatoms[j]]);
                    }
                    svmul(1.0/nrHatoms, rpar, rpar); /* rpar = ( H1+H2+H3 ) / 3 */
                    rvec_dec(rpar, (*x)[Heavy]);     /*        - Heavy          */
                    rvec_sub((*x)[Hatoms[0]], (*x)[Heavy], rperp);
                    rvec_dec(rperp, rpar);           /* rperp = H1 - Heavy - rpar */
                    /* calc mass positions */
                    svmul(fact2, rpar, temp);
                    for (j = 0; (j < NMASS); j++) /* xM = xN + fact2 * rpar +/- fact * rperp */
                    {
                        rvec_add((*x)[Heavy], temp, newx[ni0+j]);
                    }
                    svmul(fact, rperp, temp);
                    rvec_inc(newx[ni0  ], temp);
                    rvec_dec(newx[ni0+1], temp);
                    /* set atom parameters for the masses */
                    for (j = 0; (j < NMASS); j++)
                    {
                        /* make name: "M??#" or "M?#" (? is atomname, # is number) */
                        name[0] = 'M';
                        for (k = 0; (*at->atomname[Heavy])[k] && ( k < NMASS ); k++)
                        {
                            name[k+1] = (*at->atomname[Heavy])[k];
                        }
                        name[k+1]              = atomnamesuffix[j];
                        name[k+2]              = '\0';
                        newatomname[ni0+j]     = put_symtab(symtab, name);
                        newatom[ni0+j].m       = newatom[ni0+j].mB    = mtot/NMASS;
                        newatom[ni0+j].q       = newatom[ni0+j].qB    = 0.0;
                        newatom[ni0+j].type    = newatom[ni0+j].typeB = tpM;
                        newatom[ni0+j].ptype   = eptAtom;
                        newatom[ni0+j].resind  = at->atom[i0].resind;
                        newatom[ni0+j].elem[0] = 'M';
                        newatom[ni0+j].elem[1] = '\0';
                        newvsite_type[ni0+j]   = NOTSET;
                        newcgnr[ni0+j]         = (*cgnr)[i0];
                    }
                    /* add constraints between dummy masses and to heavies[0] */
                    /* 'add_shift' says which atoms won't be renumbered afterwards */
                    my_add_param(&(plist[F_CONSTRNC]), heavies[0],  add_shift+ni0,  NOTSET);
                    my_add_param(&(plist[F_CONSTRNC]), heavies[0],  add_shift+ni0+1, NOTSET);
                    my_add_param(&(plist[F_CONSTRNC]), add_shift+ni0, add_shift+ni0+1, NOTSET);

                    /* generate Heavy, H1, H2 and H3 from M1, M2 and heavies[0] */
                    /* note that vsite_type cannot be NOTSET, because we just set it */
                    add_vsite3_atoms  (&plist[(*vsite_type)[Heavy]],
                                       Heavy,     heavies[0], add_shift+ni0, add_shift+ni0+1,
                                       FALSE);
                    for (j = 0; j < nrHatoms; j++)
                    {
                        add_vsite3_atoms(&plist[(*vsite_type)[Hatoms[j]]],
                                         Hatoms[j], heavies[0], add_shift+ni0, add_shift+ni0+1,
                                         Hat_SwapParity[j]);
                    }
#undef NMASS
                }
                else
                {
                    bWARNING = TRUE;
                }

            }
            if (bWARNING)
            {
                fprintf(stderr,
                        "Warning: cannot convert atom %d %s (bound to a heavy atom "
                        "%s with \n"
                        "         %d bonds and %d bound hydrogens atoms) to virtual site\n",
                        i+1, *(at->atomname[i]), tpname, nrbonds, nrHatoms);
            }
            if (bAddVsiteParam)
            {
                /* add vsite parameters to topology,
                   also get rid of negative vsite_types */
                add_vsites(plist, (*vsite_type), Heavy, nrHatoms, Hatoms,
                           nrheavies, heavies);
                /* transfer mass of virtual site to Heavy atom */
                for (j = 0; j < nrHatoms; j++)
                {
                    if (is_vsite((*vsite_type)[Hatoms[j]]))
                    {
                        at->atom[Heavy].m    += at->atom[Hatoms[j]].m;
                        at->atom[Heavy].mB    = at->atom[Heavy].m;
                        at->atom[Hatoms[j]].m = at->atom[Hatoms[j]].mB = 0;
                    }
                }
            }
            nvsite += nrHatoms;
            if (debug)
            {
                fprintf(debug, "atom %d: ", o2n[i]+1);
                print_bonds(debug, o2n, nrHatoms, Hatoms, Heavy, nrheavies, heavies);
            }
        } /* if vsite NOTSET & is hydrogen */

    }     /* for i < at->nr */

    gmx_residuetype_destroy(rt);

    if (debug)
    {
        fprintf(debug, "Before inserting new atoms:\n");
        for (i = 0; i < at->nr; i++)
        {
            fprintf(debug, "%4d %4d %4s %4d %4s %6d %-10s\n", i+1, o2n[i]+1,
                    at->atomname[i] ? *(at->atomname[i]) : "(NULL)",
                    at->resinfo[at->atom[i].resind].nr,
                    at->resinfo[at->atom[i].resind].name ?
                    *(at->resinfo[at->atom[i].resind].name) : "(NULL)",
                    (*cgnr)[i],
                    ((*vsite_type)[i] == NOTSET) ?
                    "NOTSET" : interaction_function[(*vsite_type)[i]].name);
        }
        fprintf(debug, "new atoms to be inserted:\n");
        for (i = 0; i < at->nr+nadd; i++)
        {
            if (newatomname[i])
            {
                fprintf(debug, "%4d %4s %4d %6d %-10s\n", i+1,
                        newatomname[i] ? *(newatomname[i]) : "(NULL)",
                        newatom[i].resind, newcgnr[i],
                        (newvsite_type[i] == NOTSET) ?
                        "NOTSET" : interaction_function[newvsite_type[i]].name);
            }
        }
    }

    /* add all original atoms to the new arrays, using o2n index array */
    for (i = 0; i < at->nr; i++)
    {
        newatomname  [o2n[i]] = at->atomname [i];
        newatom      [o2n[i]] = at->atom     [i];
        newvsite_type[o2n[i]] = (*vsite_type)[i];
        newcgnr      [o2n[i]] = (*cgnr)      [i];
        copy_rvec((*x)[i], newx[o2n[i]]);
    }
    /* throw away old atoms */
    sfree(at->atom);
    sfree(at->atomname);
    sfree(*vsite_type);
    sfree(*cgnr);
    sfree(*x);
    /* put in the new ones */
    at->nr      += nadd;
    at->atom     = newatom;
    at->atomname = newatomname;
    *vsite_type  = newvsite_type;
    *cgnr        = newcgnr;
    *x           = newx;
    if (at->nr > add_shift)
    {
        gmx_fatal(FARGS, "Added impossible amount of dummy masses "
                  "(%d on a total of %d atoms)\n", nadd, at->nr-nadd);
    }

    if (debug)
    {
        fprintf(debug, "After inserting new atoms:\n");
        for (i = 0; i < at->nr; i++)
        {
            fprintf(debug, "%4d %4s %4d %4s %6d %-10s\n", i+1,
                    at->atomname[i] ? *(at->atomname[i]) : "(NULL)",
                    at->resinfo[at->atom[i].resind].nr,
                    at->resinfo[at->atom[i].resind].name ?
                    *(at->resinfo[at->atom[i].resind].name) : "(NULL)",
                    (*cgnr)[i],
                    ((*vsite_type)[i] == NOTSET) ?
                    "NOTSET" : interaction_function[(*vsite_type)[i]].name);
        }
    }

    /* now renumber all the interactions because of the added atoms */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        params = &(plist[ftype]);
        if (debug)
        {
            fprintf(debug, "Renumbering %d %s\n", params->nr,
                    interaction_function[ftype].longname);
        }
        for (i = 0; i < params->nr; i++)
        {
            for (j = 0; j < NRAL(ftype); j++)
            {
                if (params->param[i].a[j] >= add_shift)
                {
                    if (debug)
                    {
                        fprintf(debug, " [%d -> %d]", params->param[i].a[j],
                                params->param[i].a[j]-add_shift);
                    }
                    params->param[i].a[j] = params->param[i].a[j]-add_shift;
                }
                else
                {
                    if (debug)
                    {
                        fprintf(debug, " [%d -> %d]", params->param[i].a[j],
                                o2n[params->param[i].a[j]]);
                    }
                    params->param[i].a[j] = o2n[params->param[i].a[j]];
                }
            }
            if (debug)
            {
                fprintf(debug, "\n");
            }
        }
    }
    /* now check if atoms in the added constraints are in increasing order */
    params = &(plist[F_CONSTRNC]);
    for (i = 0; i < params->nr; i++)
    {
        if (params->param[i].ai() > params->param[i].aj())
        {
            j                     = params->param[i].aj();
            params->param[i].aj() = params->param[i].ai();
            params->param[i].ai() = j;
        }
    }

    /* clean up */
    sfree(o2n);

    /* tell the user what we did */
    fprintf(stderr, "Marked %d virtual sites\n", nvsite);
    fprintf(stderr, "Added %d dummy masses\n", nadd);
    fprintf(stderr, "Added %d new constraints\n", plist[F_CONSTRNC].nr);
}

void do_h_mass(t_params *psb, int vsite_type[], t_atoms *at, real mHmult,
               gmx_bool bDeuterate)
{
    int i, j, a;

    /* loop over all atoms */
    for (i = 0; i < at->nr; i++)
    {
        /* adjust masses if i is hydrogen and not a virtual site */
        if (!is_vsite(vsite_type[i]) && is_hydrogen(*(at->atomname[i])) )
        {
            /* find bonded heavy atom */
            a = NOTSET;
            for (j = 0; (j < psb->nr) && (a == NOTSET); j++)
            {
                /* if other atom is not a virtual site, it is the one we want */
                if ( (psb->param[j].ai() == i) &&
                     !is_vsite(vsite_type[psb->param[j].aj()]) )
                {
                    a = psb->param[j].aj();
                }
                else if ( (psb->param[j].aj() == i) &&
                          !is_vsite(vsite_type[psb->param[j].ai()]) )
                {
                    a = psb->param[j].ai();
                }
            }
            if (a == NOTSET)
            {
                gmx_fatal(FARGS, "Unbound hydrogen atom (%d) found while adjusting mass",
                          i+1);
            }

            /* adjust mass of i (hydrogen) with mHmult
               and correct mass of a (bonded atom) with same amount */
            if (!bDeuterate)
            {
                at->atom[a].m  -= (mHmult-1.0)*at->atom[i].m;
                at->atom[a].mB -= (mHmult-1.0)*at->atom[i].m;
            }
            at->atom[i].m  *= mHmult;
            at->atom[i].mB *= mHmult;
        }
    }
}
