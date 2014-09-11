/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "pdbio.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int ai, aj;
} gmx_conection_t;

typedef struct gmx_conect_t {
    int              nconect;
    gmx_bool         bSorted;
    gmx_conection_t *conect;
} gmx_conect_t;

static const char *pdbtp[epdbNR] = {
    "ATOM  ", "HETATM", "ANISOU", "CRYST1",
    "COMPND", "MODEL", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK",
    "CONECT"
};


#define REMARK_SIM_BOX "REMARK    THIS IS A SIMULATION BOX"

static void xlate_atomname_pdb2gmx(char *name)
{
    int  i, length;
    char temp;

    length = strlen(name);
    if (length > 3 && isdigit(name[0]))
    {
        temp = name[0];
        for (i = 1; i < length; i++)
        {
            name[i-1] = name[i];
        }
        name[length-1] = temp;
    }
}

static void xlate_atomname_gmx2pdb(char *name)
{
    int  i, length;
    char temp;

    length = strlen(name);
    if (length > 3 && isdigit(name[length-1]))
    {
        temp = name[length-1];
        for (i = length-1; i > 0; --i)
        {
            name[i] = name[i-1];
        }
        name[0] = temp;
    }
}


void gmx_write_pdb_box(FILE *out, int ePBC, matrix box)
{
    real alpha, beta, gamma;

    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }

    if (ePBC == epbcNONE)
    {
        return;
    }

    if (norm2(box[YY])*norm2(box[ZZ]) != 0)
    {
        alpha = RAD2DEG*acos(cos_angle_no_table(box[YY], box[ZZ]));
    }
    else
    {
        alpha = 90;
    }
    if (norm2(box[XX])*norm2(box[ZZ]) != 0)
    {
        beta  = RAD2DEG*acos(cos_angle_no_table(box[XX], box[ZZ]));
    }
    else
    {
        beta  = 90;
    }
    if (norm2(box[XX])*norm2(box[YY]) != 0)
    {
        gamma = RAD2DEG*acos(cos_angle_no_table(box[XX], box[YY]));
    }
    else
    {
        gamma = 90;
    }
    fprintf(out, "REMARK    THIS IS A SIMULATION BOX\n");
    if (ePBC != epbcSCREW)
    {
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                10*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 1", 1);
    }
    else
    {
        /* Double the a-vector length and write the correct space group */
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                20*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 21 1 1", 1);

    }
}

static void read_cryst1(char *line, int *ePBC, matrix box)
{
#define SG_SIZE 11
    char   sa[12], sb[12], sc[12], sg[SG_SIZE+1], ident;
    double fa, fb, fc, alpha, beta, gamma, cosa, cosb, cosg, sing;
    int    syma, symb, symc;
    int    ePBC_file;

    sscanf(line, "%*s%s%s%s%lf%lf%lf", sa, sb, sc, &alpha, &beta, &gamma);

    ePBC_file = -1;
    if (strlen(line) >= 55)
    {
        strncpy(sg, line+55, SG_SIZE);
        sg[SG_SIZE] = '\0';
        ident       = ' ';
        syma        = 0;
        symb        = 0;
        symc        = 0;
        sscanf(sg, "%c %d %d %d", &ident, &syma, &symb, &symc);
        if (ident == 'P' && syma ==  1 && symb <= 1 && symc <= 1)
        {
            fc        = strtod(sc, NULL)*0.1;
            ePBC_file = (fc > 0 ? epbcXYZ : epbcXY);
        }
        if (ident == 'P' && syma == 21 && symb == 1 && symc == 1)
        {
            ePBC_file = epbcSCREW;
        }
    }
    if (ePBC)
    {
        *ePBC = ePBC_file;
    }

    if (box)
    {
        fa = strtod(sa, NULL)*0.1;
        fb = strtod(sb, NULL)*0.1;
        fc = strtod(sc, NULL)*0.1;
        if (ePBC_file == epbcSCREW)
        {
            fa *= 0.5;
        }
        clear_mat(box);
        box[XX][XX] = fa;
        if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0))
        {
            if (alpha != 90.0)
            {
                cosa = cos(alpha*DEG2RAD);
            }
            else
            {
                cosa = 0;
            }
            if (beta != 90.0)
            {
                cosb = cos(beta*DEG2RAD);
            }
            else
            {
                cosb = 0;
            }
            if (gamma != 90.0)
            {
                cosg = cos(gamma*DEG2RAD);
                sing = sin(gamma*DEG2RAD);
            }
            else
            {
                cosg = 0;
                sing = 1;
            }
            box[YY][XX] = fb*cosg;
            box[YY][YY] = fb*sing;
            box[ZZ][XX] = fc*cosb;
            box[ZZ][YY] = fc*(cosa - cosb*cosg)/sing;
            box[ZZ][ZZ] = sqrt(fc*fc
                               - box[ZZ][XX]*box[ZZ][XX] - box[ZZ][YY]*box[ZZ][YY]);
        }
        else
        {
            box[YY][YY] = fb;
            box[ZZ][ZZ] = fc;
        }
    }
}

void write_pdbfile_indexed(FILE *out, const char *title,
                           t_atoms *atoms, rvec x[],
                           int ePBC, matrix box, char chainid,
                           int model_nr, atom_id nindex, const atom_id index[],
                           gmx_conect conect, gmx_bool bTerSepChains)
{
    gmx_conect_t     *gc = (gmx_conect_t *)conect;
    char              resnm[6], nm[6], pukestring[100];
    atom_id           i, ii;
    int               resind, resnr;
    enum PDB_record   type;
    unsigned char     resic, ch;
    char              altloc;
    real              occup, bfac;
    gmx_bool          bOccup;
    int               nlongname = 0;
    int               chainnum, lastchainnum;
    int               lastresind, lastchainresind;
    gmx_residuetype_t*rt;
    const char       *p_restype;
    const char       *p_lastrestype;

    gmx_residuetype_init(&rt);

    bromacs(pukestring, 99);
    fprintf(out, "TITLE     %s\n", (title && title[0]) ? title : pukestring);
    if (box && ( norm2(box[XX]) || norm2(box[YY]) || norm2(box[ZZ]) ) )
    {
        gmx_write_pdb_box(out, ePBC, box);
    }
    if (atoms->pdbinfo)
    {
        /* Check whether any occupancies are set, in that case leave it as is,
         * otherwise set them all to one
         */
        bOccup = TRUE;
        for (ii = 0; (ii < nindex) && bOccup; ii++)
        {
            i      = index[ii];
            bOccup = bOccup && (atoms->pdbinfo[i].occup == 0.0);
        }
    }
    else
    {
        bOccup = FALSE;
    }

    fprintf(out, "MODEL %8d\n", model_nr > 0 ? model_nr : 1);

    lastchainresind   = -1;
    lastchainnum      = -1;
    resind            = -1;
    p_restype         = NULL;

    for (ii = 0; ii < nindex; ii++)
    {
        i             = index[ii];
        lastresind    = resind;
        resind        = atoms->atom[i].resind;
        chainnum      = atoms->resinfo[resind].chainnum;
        p_lastrestype = p_restype;
        gmx_residuetype_get_type(rt, *atoms->resinfo[resind].name, &p_restype);

        /* Add a TER record if we changed chain, and if either the previous or this chain is protein/DNA/RNA. */
        if (bTerSepChains && ii > 0 && chainnum != lastchainnum)
        {
            /* Only add TER if the previous chain contained protein/DNA/RNA. */
            if (gmx_residuetype_is_protein(rt, p_lastrestype) || gmx_residuetype_is_dna(rt, p_lastrestype) || gmx_residuetype_is_rna(rt, p_lastrestype))
            {
                fprintf(out, "TER\n");
            }
            lastchainnum    = chainnum;
            lastchainresind = lastresind;
        }

        strncpy(resnm, *atoms->resinfo[resind].name, sizeof(resnm)-1);
        resnm[sizeof(resnm)-1] = 0;
        strncpy(nm, *atoms->atomname[i], sizeof(nm)-1);
        nm[sizeof(nm)-1] = 0;

        /* rename HG12 to 2HG1, etc. */
        xlate_atomname_gmx2pdb(nm);
        resnr = atoms->resinfo[resind].nr;
        resic = atoms->resinfo[resind].ic;
        if (chainid != ' ')
        {
            ch = chainid;
        }
        else
        {
            ch = atoms->resinfo[resind].chainid;

            if (ch == 0)
            {
                ch = ' ';
            }
        }
        if (resnr >= 10000)
        {
            resnr = resnr % 10000;
        }
        if (atoms->pdbinfo)
        {
            type   = (enum PDB_record)(atoms->pdbinfo[i].type);
            altloc = atoms->pdbinfo[i].altloc;
            if (!isalnum(altloc))
            {
                altloc = ' ';
            }
            occup = bOccup ? 1.0 : atoms->pdbinfo[i].occup;
            bfac  = atoms->pdbinfo[i].bfac;
        }
        else
        {
            type   = epdbATOM;
            occup  = 1.0;
            bfac   = 0.0;
            altloc = ' ';
        }

        gmx_fprintf_pdb_atomline(out,
                                 type,
                                 i+1,
                                 nm,
                                 altloc,
                                 resnm,
                                 ch,
                                 resnr,
                                 resic,
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ],
                                 occup,
                                 bfac,
                                 atoms->atom[i].elem);

        if (atoms->pdbinfo && atoms->pdbinfo[i].bAnisotropic)
        {
            fprintf(out, "ANISOU%5d  %-4.4s%4.4s%c%4d%c %7d%7d%7d%7d%7d%7d\n",
                    (i+1)%100000, nm, resnm, ch, resnr,
                    (resic == '\0') ? ' ' : resic,
                    atoms->pdbinfo[i].uij[0], atoms->pdbinfo[i].uij[1],
                    atoms->pdbinfo[i].uij[2], atoms->pdbinfo[i].uij[3],
                    atoms->pdbinfo[i].uij[4], atoms->pdbinfo[i].uij[5]);
        }
    }

    fprintf(out, "TER\n");
    fprintf(out, "ENDMDL\n");

    if (NULL != gc)
    {
        /* Write conect records */
        for (i = 0; (i < gc->nconect); i++)
        {
            fprintf(out, "CONECT%5d%5d\n", gc->conect[i].ai+1, gc->conect[i].aj+1);
        }
    }

    gmx_residuetype_destroy(rt);
}

void write_pdbfile(FILE *out, const char *title, t_atoms *atoms, rvec x[],
                   int ePBC, matrix box, char chainid, int model_nr, gmx_conect conect, gmx_bool bTerSepChains)
{
    atom_id i, *index;

    snew(index, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        index[i] = i;
    }
    write_pdbfile_indexed(out, title, atoms, x, ePBC, box, chainid, model_nr,
                          atoms->nr, index, conect, bTerSepChains);
    sfree(index);
}

int line2type(char *line)
{
    int  k;
    char type[8];

    for (k = 0; (k < 6); k++)
    {
        type[k] = line[k];
    }
    type[k] = '\0';

    for (k = 0; (k < epdbNR); k++)
    {
        if (strncmp(type, pdbtp[k], strlen(pdbtp[k])) == 0)
        {
            break;
        }
    }

    return k;
}

static void read_anisou(char line[], int natom, t_atoms *atoms)
{
    int  i, j, k, atomnr;
    char nc = '\0';
    char anr[12], anm[12];

    /* Skip over type */
    j = 6;
    for (k = 0; (k < 5); k++, j++)
    {
        anr[k] = line[j];
    }
    anr[k] = nc;
    j++;
    for (k = 0; (k < 4); k++, j++)
    {
        anm[k] = line[j];
    }
    anm[k] = nc;
    j++;

    /* Strip off spaces */
    trim(anm);

    /* Search backwards for number and name only */
    atomnr = strtol(anr, NULL, 10);
    for (i = natom-1; (i >= 0); i--)
    {
        if ((strcmp(anm, *(atoms->atomname[i])) == 0) &&
            (atomnr == atoms->pdbinfo[i].atomnr))
        {
            break;
        }
    }
    if (i < 0)
    {
        fprintf(stderr, "Skipping ANISOU record (atom %s %d not found)\n",
                anm, atomnr);
    }
    else
    {
        if (sscanf(line+29, "%d%d%d%d%d%d",
                   &atoms->pdbinfo[i].uij[U11], &atoms->pdbinfo[i].uij[U22],
                   &atoms->pdbinfo[i].uij[U33], &atoms->pdbinfo[i].uij[U12],
                   &atoms->pdbinfo[i].uij[U13], &atoms->pdbinfo[i].uij[U23])
            == 6)
        {
            atoms->pdbinfo[i].bAnisotropic = TRUE;
        }
        else
        {
            fprintf(stderr, "Invalid ANISOU record for atom %d\n", i);
            atoms->pdbinfo[i].bAnisotropic = FALSE;
        }
    }
}

void get_pdb_atomnumber(t_atoms *atoms, gmx_atomprop_t aps)
{
    int    i, atomnumber, len;
    size_t k;
    char   anm[6], anm_copy[6], *ptr;
    char   nc = '\0';
    real   eval;

    if (!atoms->pdbinfo)
    {
        gmx_incons("Trying to deduce atomnumbers when no pdb information is present");
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        strcpy(anm, atoms->pdbinfo[i].atomnm);
        strcpy(anm_copy, atoms->pdbinfo[i].atomnm);
        len        = strlen(anm);
        atomnumber = NOTSET;
        if ((anm[0] != ' ') && ((len <= 2) || ((len > 2) && !isdigit(anm[2]))))
        {
            anm_copy[2] = nc;
            if (gmx_atomprop_query(aps, epropElement, "???", anm_copy, &eval))
            {
                atomnumber = gmx_nint(eval);
            }
            else
            {
                anm_copy[1] = nc;
                if (gmx_atomprop_query(aps, epropElement, "???", anm_copy, &eval))
                {
                    atomnumber = gmx_nint(eval);
                }
            }
        }
        if (atomnumber == NOTSET)
        {
            k = 0;
            while ((k < strlen(anm)) && (isspace(anm[k]) || isdigit(anm[k])))
            {
                k++;
            }
            anm_copy[0] = anm[k];
            anm_copy[1] = nc;
            if (gmx_atomprop_query(aps, epropElement, "???", anm_copy, &eval))
            {
                atomnumber = gmx_nint(eval);
            }
        }
        atoms->atom[i].atomnumber = atomnumber;
        ptr = gmx_atomprop_element(aps, atomnumber);
        strncpy(atoms->atom[i].elem, ptr == NULL ? "" : ptr, 4);
        if (debug)
        {
            fprintf(debug, "Atomnumber for atom '%s' is %d\n", anm, atomnumber);
        }
    }
}

static int read_atom(t_symtab *symtab,
                     char line[], int type, int natom,
                     t_atoms *atoms, rvec x[], int chainnum, gmx_bool bChange)
{
    t_atom       *atomn;
    int           j, k;
    char          nc = '\0';
    char          anr[12], anm[12], anm_copy[12], altloc, resnm[12], rnr[12], elem[3];
    char          xc[12], yc[12], zc[12], occup[12], bfac[12];
    unsigned char resic;
    char          chainid;
    int           resnr, atomnumber;

    if (natom >= atoms->nr)
    {
        gmx_fatal(FARGS, "\nFound more atoms (%d) in pdb file than expected (%d)",
                  natom+1, atoms->nr);
    }

    /* Skip over type */
    j = 6;
    for (k = 0; (k < 5); k++, j++)
    {
        anr[k] = line[j];
    }
    anr[k] = nc;
    trim(anr);
    j++;
    for (k = 0; (k < 4); k++, j++)
    {
        anm[k] = line[j];
    }
    anm[k] = nc;
    strcpy(anm_copy, anm);
    rtrim(anm_copy);
    atomnumber = NOTSET;
    trim(anm);
    altloc = line[j];
    j++;
    for (k = 0; (k < 4); k++, j++)
    {
        resnm[k] = line[j];
    }
    resnm[k] = nc;
    trim(resnm);

    chainid = line[j];
    j++;

    for (k = 0; (k < 4); k++, j++)
    {
        rnr[k] = line[j];
    }
    rnr[k] = nc;
    trim(rnr);
    resnr = strtol(rnr, NULL, 10);
    resic = line[j];
    j    += 4;

    /* X,Y,Z Coordinate */
    for (k = 0; (k < 8); k++, j++)
    {
        xc[k] = line[j];
    }
    xc[k] = nc;
    for (k = 0; (k < 8); k++, j++)
    {
        yc[k] = line[j];
    }
    yc[k] = nc;
    for (k = 0; (k < 8); k++, j++)
    {
        zc[k] = line[j];
    }
    zc[k] = nc;

    /* Weight */
    for (k = 0; (k < 6); k++, j++)
    {
        occup[k] = line[j];
    }
    occup[k] = nc;

    /* B-Factor */
    for (k = 0; (k < 7); k++, j++)
    {
        bfac[k] = line[j];
    }
    bfac[k] = nc;

    /* 10 blanks */
    j += 10;

    /* Element name */
    for (k = 0; (k < 2); k++, j++)
    {
        elem[k] = line[j];
    }
    elem[k] = nc;
    trim(elem);

    if (atoms->atom)
    {
        atomn = &(atoms->atom[natom]);
        if ((natom == 0) ||
            atoms->resinfo[atoms->atom[natom-1].resind].nr != resnr ||
            atoms->resinfo[atoms->atom[natom-1].resind].ic != resic ||
            (strcmp(*atoms->resinfo[atoms->atom[natom-1].resind].name, resnm) != 0))
        {
            if (natom == 0)
            {
                atomn->resind = 0;
            }
            else
            {
                atomn->resind = atoms->atom[natom-1].resind + 1;
            }
            atoms->nres = atomn->resind + 1;
            t_atoms_set_resinfo(atoms, natom, symtab, resnm, resnr, resic, chainnum, chainid);
        }
        else
        {
            atomn->resind = atoms->atom[natom-1].resind;
        }
        if (bChange)
        {
            xlate_atomname_pdb2gmx(anm);
        }
        atoms->atomname[natom] = put_symtab(symtab, anm);
        atomn->m               = 0.0;
        atomn->q               = 0.0;
        atomn->atomnumber      = atomnumber;
        strncpy(atomn->elem, elem, 4);
    }
    x[natom][XX] = strtod(xc, NULL)*0.1;
    x[natom][YY] = strtod(yc, NULL)*0.1;
    x[natom][ZZ] = strtod(zc, NULL)*0.1;
    if (atoms->pdbinfo)
    {
        atoms->pdbinfo[natom].type   = type;
        atoms->pdbinfo[natom].atomnr = strtol(anr, NULL, 10);
        atoms->pdbinfo[natom].altloc = altloc;
        strcpy(atoms->pdbinfo[natom].atomnm, anm_copy);
        atoms->pdbinfo[natom].bfac  = strtod(bfac, NULL);
        atoms->pdbinfo[natom].occup = strtod(occup, NULL);
    }
    natom++;

    return natom;
}

gmx_bool is_hydrogen(const char *nm)
{
    char buf[30];

    strcpy(buf, nm);
    trim(buf);

    if (buf[0] == 'H')
    {
        return TRUE;
    }
    else if ((isdigit(buf[0])) && (buf[1] == 'H'))
    {
        return TRUE;
    }
    return FALSE;
}

gmx_bool is_dummymass(const char *nm)
{
    char buf[30];

    strcpy(buf, nm);
    trim(buf);

    if ((buf[0] == 'M') && isdigit(buf[strlen(buf)-1]))
    {
        return TRUE;
    }

    return FALSE;
}

static void gmx_conect_addline(gmx_conect_t *con, char *line)
{
    int  n, ai, aj;
    char format[32], form2[32];

    sprintf(form2, "%%*s");
    sprintf(format, "%s%%d", form2);
    if (sscanf(line, format, &ai) == 1)
    {
        do
        {
            strcat(form2, "%*s");
            sprintf(format, "%s%%d", form2);
            n = sscanf(line, format, &aj);
            if (n == 1)
            {
                srenew(con->conect, ++con->nconect);
                con->conect[con->nconect-1].ai = ai-1;
                con->conect[con->nconect-1].aj = aj-1;
            }
        }
        while (n == 1);
    }
}

void gmx_conect_dump(FILE *fp, gmx_conect conect)
{
    gmx_conect_t *gc = (gmx_conect_t *)conect;
    int           i;

    for (i = 0; (i < gc->nconect); i++)
    {
        fprintf(fp, "%6s%5d%5d\n", "CONECT",
                gc->conect[i].ai+1, gc->conect[i].aj+1);
    }
}

gmx_conect gmx_conect_init()
{
    gmx_conect_t *gc;

    snew(gc, 1);

    return (gmx_conect) gc;
}

void gmx_conect_done(gmx_conect conect)
{
    gmx_conect_t *gc = (gmx_conect_t *)conect;

    sfree(gc->conect);
}

gmx_bool gmx_conect_exist(gmx_conect conect, int ai, int aj)
{
    gmx_conect_t *gc = (gmx_conect_t *)conect;
    int           i;

    /* if (!gc->bSorted)
       sort_conect(gc);*/

    for (i = 0; (i < gc->nconect); i++)
    {
        if (((gc->conect[i].ai == ai) &&
             (gc->conect[i].aj == aj)) ||
            ((gc->conect[i].aj == ai) &&
             (gc->conect[i].ai == aj)))
        {
            return TRUE;
        }
    }
    return FALSE;
}

void gmx_conect_add(gmx_conect conect, int ai, int aj)
{
    gmx_conect_t *gc = (gmx_conect_t *)conect;
    int           i;

    /* if (!gc->bSorted)
       sort_conect(gc);*/

    if (!gmx_conect_exist(conect, ai, aj))
    {
        srenew(gc->conect, ++gc->nconect);
        gc->conect[gc->nconect-1].ai = ai;
        gc->conect[gc->nconect-1].aj = aj;
    }
}

int read_pdbfile(FILE *in, char *title, int *model_nr,
                 t_atoms *atoms, rvec x[], int *ePBC, matrix box, gmx_bool bChange,
                 gmx_conect conect)
{
    gmx_conect_t *gc = (gmx_conect_t *)conect;
    t_symtab      symtab;
    gmx_bool      bCOMPND;
    gmx_bool      bConnWarn = FALSE;
    char          line[STRLEN+1];
    int           line_type;
    char         *c, *d;
    int           natom, chainnum, nres_ter_prev = 0;
    char          chidmax = ' ';
    gmx_bool      bStop   = FALSE;

    if (ePBC)
    {
        /* Only assume pbc when there is a CRYST1 entry */
        *ePBC = epbcNONE;
    }
    if (box != NULL)
    {
        clear_mat(box);
    }

    open_symtab(&symtab);

    bCOMPND  = FALSE;
    title[0] = '\0';
    natom    = 0;
    chainnum = 0;
    while (!bStop && (fgets2(line, STRLEN, in) != NULL))
    {
        line_type = line2type(line);

        switch (line_type)
        {
            case epdbATOM:
            case epdbHETATM:
                natom = read_atom(&symtab, line, line_type, natom, atoms, x, chainnum, bChange);
                break;

            case epdbANISOU:
                if (atoms->pdbinfo)
                {
                    read_anisou(line, natom, atoms);
                }
                break;

            case epdbCRYST1:
                read_cryst1(line, ePBC, box);
                break;

            case epdbTITLE:
            case epdbHEADER:
                if (strlen(line) > 6)
                {
                    c = line+6;
                    /* skip HEADER or TITLE and spaces */
                    while (c[0] != ' ')
                    {
                        c++;
                    }
                    while (c[0] == ' ')
                    {
                        c++;
                    }
                    /* truncate after title */
                    d = strstr(c, "      ");
                    if (d)
                    {
                        d[0] = '\0';
                    }
                    if (strlen(c) > 0)
                    {
                        strcpy(title, c);
                    }
                }
                break;

            case epdbCOMPND:
                if ((!strstr(line, ": ")) || (strstr(line+6, "MOLECULE:")))
                {
                    if (!(c = strstr(line+6, "MOLECULE:")) )
                    {
                        c = line;
                    }
                    /* skip 'MOLECULE:' and spaces */
                    while (c[0] != ' ')
                    {
                        c++;
                    }
                    while (c[0] == ' ')
                    {
                        c++;
                    }
                    /* truncate after title */
                    d = strstr(c, "   ");
                    if (d)
                    {
                        while ( (d[-1] == ';') && d > c)
                        {
                            d--;
                        }
                        d[0] = '\0';
                    }
                    if (strlen(c) > 0)
                    {
                        if (bCOMPND)
                        {
                            strcat(title, "; ");
                            strcat(title, c);
                        }
                        else
                        {
                            strcpy(title, c);
                        }
                    }
                    bCOMPND = TRUE;
                }
                break;

            case epdbTER:
                chainnum++;
                break;

            case epdbMODEL:
                if (model_nr)
                {
                    sscanf(line, "%*s%d", model_nr);
                }
                break;

            case epdbENDMDL:
                bStop = TRUE;
                break;
            case epdbCONECT:
                if (gc)
                {
                    gmx_conect_addline(gc, line);
                }
                else if (!bConnWarn)
                {
                    fprintf(stderr, "WARNING: all CONECT records are ignored\n");
                    bConnWarn = TRUE;
                }
                break;

            default:
                break;
        }
    }

    free_symtab(&symtab);
    return natom;
}

void get_pdb_coordnum(FILE *in, int *natoms)
{
    char line[STRLEN];

    *natoms = 0;
    while (fgets2(line, STRLEN, in))
    {
        if (strncmp(line, "ENDMDL", 6) == 0)
        {
            break;
        }
        if ((strncmp(line, "ATOM  ", 6) == 0) || (strncmp(line, "HETATM", 6) == 0))
        {
            (*natoms)++;
        }
    }
}

void read_pdb_conf(const char *infile, char *title,
                   t_atoms *atoms, rvec x[], int *ePBC, matrix box, gmx_bool bChange,
                   gmx_conect conect)
{
    FILE *in;

    in = gmx_fio_fopen(infile, "r");
    read_pdbfile(in, title, NULL, atoms, x, ePBC, box, bChange, conect);
    gmx_fio_fclose(in);
}

gmx_conect gmx_conect_generate(t_topology *top)
{
    int        f, i;
    gmx_conect gc;

    /* Fill the conect records */
    gc  = gmx_conect_init();

    for (f = 0; (f < F_NRE); f++)
    {
        if (IS_CHEMBOND(f))
        {
            for (i = 0; (i < top->idef.il[f].nr); i += interaction_function[f].nratoms+1)
            {
                gmx_conect_add(gc, top->idef.il[f].iatoms[i+1],
                               top->idef.il[f].iatoms[i+2]);
            }
        }
    }
    return gc;
}

int
gmx_fprintf_pdb_atomline(FILE *            fp,
                         enum PDB_record   record,
                         int               atom_seq_number,
                         const char *      atom_name,
                         char              alternate_location,
                         const char *      res_name,
                         char              chain_id,
                         int               res_seq_number,
                         char              res_insertion_code,
                         real              x,
                         real              y,
                         real              z,
                         real              occupancy,
                         real              b_factor,
                         const char *      element)
{
    char     tmp_atomname[6], tmp_resname[6];
    gmx_bool start_name_in_col13;
    int      n;

    if (record != epdbATOM && record != epdbHETATM)
    {
        gmx_fatal(FARGS, "Can only print PDB atom lines as ATOM or HETATM records");
    }

    /* Format atom name */
    if (atom_name != NULL)
    {
        /* If the atom name is an element name with two chars, it should start already in column 13.
         * Otherwise it should start in column 14, unless the name length is 4 chars.
         */
        if ( (element != NULL) && (strlen(element) >= 2) && (gmx_strncasecmp(atom_name, element, 2) == 0) )
        {
            start_name_in_col13 = TRUE;
        }
        else
        {
            start_name_in_col13 = (strlen(atom_name) >= 4);
        }
        sprintf(tmp_atomname, start_name_in_col13 ? "" : " ");
        strncat(tmp_atomname, atom_name, 4);
        tmp_atomname[5] = '\0';
    }
    else
    {
        tmp_atomname[0] = '\0';
    }

    /* Format residue name */
    strncpy(tmp_resname, (res_name != NULL) ? res_name : "", 4);
    /* Make sure the string is terminated if strlen was > 4 */
    tmp_resname[4] = '\0';
    /* String is properly terminated, so now we can use strcat. By adding a
     * space we can write it right-justified, and if the original name was
     * three characters or less there will be a space added on the right side.
     */
    strcat(tmp_resname, " ");

    /* Truncate integers so they fit */
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    n = fprintf(fp,
                "%-6s%5d %-4.4s%c%4.4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                pdbtp[record],
                atom_seq_number,
                tmp_atomname,
                alternate_location,
                tmp_resname,
                chain_id,
                res_seq_number,
                res_insertion_code,
                x, y, z,
                occupancy,
                b_factor,
                (element != NULL) ? element : "");

    return n;
}
