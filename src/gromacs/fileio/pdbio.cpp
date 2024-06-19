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

#include "gromacs/fileio/pdbio.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringtoenumvalueconverter.h"

typedef struct
{
    int ai, aj;
} gmx_conection_t;

typedef struct gmx_conect_t
{
    int              nconect;
    gmx_conection_t* conect;
} gmx_conect_t;

const char* enumValueToString(PdbRecordType enumValue)
{
    static constexpr gmx::EnumerationArray<PdbRecordType, const char*> pdbRecordTypeName = {
        "ATOM  ", "HETATM", "ANISOU", "CRYST1", "COMPND", "MODEL",
        "ENDMDL", "TER",    "HEADER", "TITLE",  "REMARK", "CONECT"
    };
    return pdbRecordTypeName[enumValue];
}

void gmx_write_pdb_box(FILE* out, PbcType pbcType, const matrix box)
{
    real alpha, beta, gamma;

    if (pbcType == PbcType::Unset)
    {
        pbcType = guessPbcType(box);
    }

    if (pbcType == PbcType::No)
    {
        return;
    }

    if (norm2(box[YY]) * norm2(box[ZZ]) != 0)
    {
        alpha = gmx::c_rad2Deg * gmx_angle(box[YY], box[ZZ]);
    }
    else
    {
        alpha = 90;
    }
    if (norm2(box[XX]) * norm2(box[ZZ]) != 0)
    {
        beta = gmx::c_rad2Deg * gmx_angle(box[XX], box[ZZ]);
    }
    else
    {
        beta = 90;
    }
    if (norm2(box[XX]) * norm2(box[YY]) != 0)
    {
        gamma = gmx::c_rad2Deg * gmx_angle(box[XX], box[YY]);
    }
    else
    {
        gamma = 90;
    }
    fprintf(out, "REMARK    THIS IS A SIMULATION BOX\n");
    if (pbcType != PbcType::Screw)
    {
        fprintf(out,
                "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                10 * norm(box[XX]),
                10 * norm(box[YY]),
                10 * norm(box[ZZ]),
                alpha,
                beta,
                gamma,
                "P 1",
                1);
    }
    else
    {
        /* Double the a-vector length and write the correct space group */
        fprintf(out,
                "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                20 * norm(box[XX]),
                10 * norm(box[YY]),
                10 * norm(box[ZZ]),
                alpha,
                beta,
                gamma,
                "P 21 1 1",
                1);
    }
}

static void read_cryst1(char* line, PbcType* pbcType, matrix box)
{
#define SG_SIZE 11
    char    sa[12], sb[12], sc[12], sg[SG_SIZE + 1], ident;
    double  fa, fb, fc, alpha, beta, gamma, cosa, cosb, cosg, sing;
    int     syma, symb, symc;
    PbcType pbcTypeFile;

    sscanf(line, "%*s%s%s%s%lf%lf%lf", sa, sb, sc, &alpha, &beta, &gamma);

    pbcTypeFile = PbcType::Unset;
    if (strlen(line) >= 55)
    {
        strncpy(sg, line + 55, SG_SIZE);
        sg[SG_SIZE] = '\0';
        ident       = ' ';
        syma        = 0;
        symb        = 0;
        symc        = 0;
        sscanf(sg, "%c %d %d %d", &ident, &syma, &symb, &symc);
        if (ident == 'P' && syma == 1 && symb <= 1 && symc <= 1)
        {
            fc          = strtod(sc, nullptr) * 0.1;
            pbcTypeFile = (fc > 0 ? PbcType::Xyz : PbcType::XY);
        }
        if (ident == 'P' && syma == 21 && symb == 1 && symc == 1)
        {
            pbcTypeFile = PbcType::Screw;
        }
    }
    if (pbcType)
    {
        *pbcType = pbcTypeFile;
    }

    if (box)
    {
        fa = strtod(sa, nullptr) * 0.1;
        fb = strtod(sb, nullptr) * 0.1;
        fc = strtod(sc, nullptr) * 0.1;
        if (pbcTypeFile == PbcType::Screw)
        {
            fa *= 0.5;
        }
        clear_mat(box);
        box[XX][XX] = fa;
        if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0))
        {
            if (alpha != 90.0)
            {
                cosa = std::cos(alpha * gmx::c_deg2Rad);
            }
            else
            {
                cosa = 0;
            }
            if (beta != 90.0)
            {
                cosb = std::cos(beta * gmx::c_deg2Rad);
            }
            else
            {
                cosb = 0;
            }
            if (gamma != 90.0)
            {
                cosg = std::cos(gamma * gmx::c_deg2Rad);
                sing = std::sin(gamma * gmx::c_deg2Rad);
            }
            else
            {
                cosg = 0;
                sing = 1;
            }
            box[YY][XX] = fb * cosg;
            box[YY][YY] = fb * sing;
            box[ZZ][XX] = fc * cosb;
            box[ZZ][YY] = fc * (cosa - cosb * cosg) / sing;
            box[ZZ][ZZ] = std::sqrt(fc * fc - box[ZZ][XX] * box[ZZ][XX] - box[ZZ][YY] * box[ZZ][YY]);
        }
        else
        {
            box[YY][YY] = fb;
            box[ZZ][ZZ] = fc;
        }
    }
}

static int gmx_fprintf_pqr_atomline(FILE*         fp,
                                    PdbRecordType record,
                                    int           atom_seq_number,
                                    const char*   atom_name,
                                    const char*   res_name,
                                    char          chain_id,
                                    int           res_seq_number,
                                    real          x,
                                    real          y,
                                    real          z,
                                    real          occupancy,
                                    real          b_factor)
{
    GMX_RELEASE_ASSERT(record == PdbRecordType::Atom || record == PdbRecordType::Hetatm,
                       "Can only print PQR atom lines as ATOM or HETATM records");

    /* Check atom name */
    GMX_RELEASE_ASSERT(atom_name != nullptr, "Need atom information to print pqr");

    /* Check residue name */
    GMX_RELEASE_ASSERT(res_name != nullptr, "Need residue information to print pqr");

    /* Truncate integers so they fit */
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    int n = fprintf(fp,
                    "%-6s%5d %-4.4s%4.4s%c%4d %8.3f %8.3f %8.3f %6.2f %6.2f\n",
                    enumValueToString(record),
                    atom_seq_number,
                    atom_name,
                    res_name,
                    chain_id,
                    res_seq_number,
                    x,
                    y,
                    z,
                    occupancy,
                    b_factor);

    return n;
}

void write_pdbfile_indexed(FILE*          out,
                           const char*    title,
                           const t_atoms* atoms,
                           const rvec     x[],
                           PbcType        pbcType,
                           const matrix   box,
                           char           chainid,
                           int            model_nr,
                           int            nindex,
                           const int      index[],
                           gmx_conect     conect,
                           bool           usePqrFormat,
                           bool           standardCompliantMode)
{
    gmx_conect_t* gc = static_cast<gmx_conect_t*>(conect);
    PdbRecordType type;
    char          altloc;
    real          occup, bfac;
    gmx_bool      bOccup;

    if (standardCompliantMode)
    {
        fprintf(out, "HEADER    GROMACS SIMULATION BOX                  01-JAN-00   0000\n");
        fprintf(out,
                "TITLE     Gromacs simulation box\n"
                "COMPND    MOL_ID:  1;                                                           \n"
                "COMPND   2 MOLECULE:  GROMACS SIMULATION BOX;                                   \n"
                "COMPND   3 CHAIN: A;  \n"
                "SOURCE    MOL_ID: 1;\n"
                "SOURCE   2 SYNTHETIC\n"
                "KEYWDS    GROMACS\n"
                "EXPDTA    PURE PRODUCT OF COMPUTER SIMULATION\n"
                "AUTHOR    GROMACS\n"
                "REVDAT   1   01-JAN-00 0000    0\n");
    }
    else
    {
        fprintf(out, "TITLE     %s\n", (title && title[0]) ? title : gmx::bromacs().c_str());
    }


    if (box && ((norm2(box[XX]) != 0.0F) || (norm2(box[YY]) != 0.0F) || (norm2(box[ZZ]) != 0.0F)))
    {
        gmx_write_pdb_box(out, pbcType, box);
    }
    if (atoms->havePdbInfo)
    {
        /* Check whether any occupancies are set, in that case leave it as is,
         * otherwise set them all to one
         */
        bOccup = TRUE;
        for (int ii = 0; (ii < nindex) && bOccup; ii++)
        {
            int i  = index[ii];
            bOccup = bOccup && (atoms->pdbinfo[i].occup == 0.0);
        }
    }
    else
    {
        bOccup = FALSE;
    }

    fprintf(out, "MODEL %8d\n", model_nr > 0 ? model_nr : 1);

    // Collect last printed values for TER record
    int         lastAtomNumber = 0;
    std::string lastResName;
    int         lastResNr = 0;

    for (int ii = 0; ii < nindex; ii++)
    {
        int         i      = index[ii];
        int         resind = atoms->atom[i].resind;
        std::string resnm  = *atoms->resinfo[resind].name;
        std::string nm     = *atoms->atomname[i];

        int           resnr = atoms->resinfo[resind].nr;
        unsigned char resic = atoms->resinfo[resind].ic;
        unsigned char ch;
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
        t_pdbinfo pdbinfo;
        if (atoms->pdbinfo != nullptr)
        {
            pdbinfo = atoms->pdbinfo[i];
        }
        else
        {
            gmx_pdbinfo_init_default(&pdbinfo);
        }
        type   = pdbinfo.type;
        altloc = pdbinfo.altloc;
        if (!isalnum(altloc))
        {
            altloc = ' ';
        }
        occup = bOccup ? 1.0 : pdbinfo.occup;
        bfac  = pdbinfo.bfac;
        if (!usePqrFormat)
        {
            gmx_fprintf_pdb_atomline(out,
                                     type,
                                     i + 1,
                                     nm.c_str(),
                                     altloc,
                                     resnm.c_str(),
                                     ch,
                                     resnr,
                                     resic,
                                     10 * x[i][XX],
                                     10 * x[i][YY],
                                     10 * x[i][ZZ],
                                     occup,
                                     bfac,
                                     atoms->atom[i].elem);

            lastAtomNumber = i + 1;
            lastResName    = resnm;
            lastResNr      = resnr;

            if (atoms->pdbinfo && atoms->pdbinfo[i].bAnisotropic)
            {
                fprintf(out,
                        "ANISOU%5d  %-4.4s%4.4s%c%4d%c %7d%7d%7d%7d%7d%7d\n",
                        (i + 1) % 100000,
                        nm.c_str(),
                        resnm.c_str(),
                        ch,
                        resnr,
                        (resic == '\0') ? ' ' : resic,
                        atoms->pdbinfo[i].uij[0],
                        atoms->pdbinfo[i].uij[1],
                        atoms->pdbinfo[i].uij[2],
                        atoms->pdbinfo[i].uij[3],
                        atoms->pdbinfo[i].uij[4],
                        atoms->pdbinfo[i].uij[5]);
            }
        }
        else
        {
            gmx_fprintf_pqr_atomline(out,
                                     type,
                                     i + 1,
                                     nm.c_str(),
                                     resnm.c_str(),
                                     ch,
                                     resnr,
                                     10 * x[i][XX],
                                     10 * x[i][YY],
                                     10 * x[i][ZZ],
                                     occup,
                                     bfac);
        }
    }

    if (standardCompliantMode)
    {
        fprintf(out, "TER   %5d      %4.4s%c%4d\n", lastAtomNumber, lastResName.c_str(), chainid, lastResNr);
    }
    else
    {
        fprintf(out, "TER\n");
    }

    fprintf(out, "ENDMDL\n");

    if (nullptr != gc)
    {
        /* Write conect records */
        for (int i = 0; (i < gc->nconect); i++)
        {
            fprintf(out, "CONECT%5d%5d\n", gc->conect[i].ai + 1, gc->conect[i].aj + 1);
        }
    }
}

void write_pdbfile(FILE*          out,
                   const char*    title,
                   const t_atoms* atoms,
                   const rvec     x[],
                   PbcType        pbcType,
                   const matrix   box,
                   char           chainid,
                   int            model_nr,
                   gmx_conect     conect)
{
    int i, *index;

    snew(index, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        index[i] = i;
    }
    write_pdbfile_indexed(
            out, title, atoms, x, pbcType, box, chainid, model_nr, atoms->nr, index, conect, false);
    sfree(index);
}

static void read_anisou(char line[], int natom, t_atoms* atoms)
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
    atomnr = std::strtol(anr, nullptr, 10);
    for (i = natom - 1; (i >= 0); i--)
    {
        if ((std::strcmp(anm, *(atoms->atomname[i])) == 0) && (atomnr == atoms->pdbinfo[i].atomnr))
        {
            break;
        }
    }
    if (i < 0)
    {
        fprintf(stderr, "Skipping ANISOU record (atom %s %d not found)\n", anm, atomnr);
    }
    else
    {
        if (sscanf(line + 29,
                   "%d%d%d%d%d%d",
                   &atoms->pdbinfo[i].uij[U11],
                   &atoms->pdbinfo[i].uij[U22],
                   &atoms->pdbinfo[i].uij[U33],
                   &atoms->pdbinfo[i].uij[U12],
                   &atoms->pdbinfo[i].uij[U13],
                   &atoms->pdbinfo[i].uij[U23])
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

void get_pdb_atomnumber(const t_atoms* atoms, AtomProperties* aps)
{
    int    i, atomnumber, len;
    size_t k;
    char   anm[6], anm_copy[6];
    char   nc = '\0';
    real   eval;

    if (!atoms->pdbinfo)
    {
        gmx_incons("Trying to deduce atomnumbers when no pdb information is present");
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        std::strcpy(anm, atoms->pdbinfo[i].atomnm);
        std::strcpy(anm_copy, atoms->pdbinfo[i].atomnm);
        bool atomNumberSet = false;
        len                = strlen(anm);
        if ((anm[0] != ' ') && ((len <= 2) || !std::isdigit(anm[2])))
        {
            anm_copy[2] = nc;
            if (aps->setAtomProperty(epropElement, "???", anm_copy, &eval))
            {
                atomnumber    = gmx::roundToInt(eval);
                atomNumberSet = true;
            }
            else
            {
                anm_copy[1] = nc;
                if (aps->setAtomProperty(epropElement, "???", anm_copy, &eval))
                {
                    atomnumber    = gmx::roundToInt(eval);
                    atomNumberSet = true;
                }
            }
        }
        if (!atomNumberSet)
        {
            k = 0;
            while ((k < std::strlen(anm)) && (std::isspace(anm[k]) || std::isdigit(anm[k])))
            {
                k++;
            }
            anm_copy[0] = anm[k];
            anm_copy[1] = nc;
            if (aps->setAtomProperty(epropElement, "???", anm_copy, &eval))
            {
                atomnumber    = gmx::roundToInt(eval);
                atomNumberSet = true;
            }
        }
        static constexpr size_t sc_maxElementNameLength = 3;
        static_assert(sizeof(atoms->atom[i].elem) >= sc_maxElementNameLength + 1);
        std::string element;
        if (atomNumberSet)
        {
            atoms->atom[i].atomnumber = atomnumber;
            element                   = aps->elementFromAtomNumber(atomnumber);
            if (debug)
            {
                fprintf(debug, "Atomnumber for atom '%s' is %d\n", anm, atomnumber);
            }
        }
        element.resize(sc_maxElementNameLength);
        std::strcpy(atoms->atom[i].elem, element.c_str());
    }
}

static int
read_atom(t_symtab* symtab, const char line[], PdbRecordType type, int natom, t_atoms* atoms, rvec x[], int chainnum)
{
    t_atom*       atomn;
    int           j, k;
    char          nc = '\0';
    char          anr[12], anm[12], anm_copy[12], altloc, resnm[12], rnr[12], elem[3];
    char          xc[12], yc[12], zc[12], occup[12], bfac[12];
    unsigned char resic;
    char          chainid;
    int           resnr, atomnumber;

    if (natom >= atoms->nr)
    {
        gmx_fatal(FARGS, "\nFound more atoms (%d) in pdb file than expected (%d)", natom + 1, atoms->nr);
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
    std::strcpy(anm_copy, anm);
    rtrim(anm_copy);
    atomnumber = 0;
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
    resnr = std::strtol(rnr, nullptr, 10);
    resic = line[j];
    j += 4;

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
        if ((natom == 0) || atoms->resinfo[atoms->atom[natom - 1].resind].nr != resnr
            || atoms->resinfo[atoms->atom[natom - 1].resind].ic != resic
            || (strcmp(*atoms->resinfo[atoms->atom[natom - 1].resind].name, resnm) != 0))
        {
            if (natom == 0)
            {
                atomn->resind = 0;
            }
            else
            {
                atomn->resind = atoms->atom[natom - 1].resind + 1;
            }
            atoms->nres = atomn->resind + 1;
            t_atoms_set_resinfo(atoms, natom, symtab, resnm, resnr, resic, chainnum, chainid);
        }
        else
        {
            atomn->resind = atoms->atom[natom - 1].resind;
        }
        atoms->atomname[natom] = put_symtab(symtab, anm);
        atomn->m               = 0.0;
        atomn->q               = 0.0;
        atomn->atomnumber      = atomnumber;
        strncpy(atomn->elem, elem, 4);
    }
    x[natom][XX] = strtod(xc, nullptr) * 0.1;
    x[natom][YY] = strtod(yc, nullptr) * 0.1;
    x[natom][ZZ] = strtod(zc, nullptr) * 0.1;
    if (atoms->pdbinfo)
    {
        atoms->pdbinfo[natom].type   = type;
        atoms->pdbinfo[natom].atomnr = strtol(anr, nullptr, 10);
        atoms->pdbinfo[natom].altloc = altloc;
        strcpy(atoms->pdbinfo[natom].atomnm, anm_copy);
        atoms->pdbinfo[natom].bfac  = strtod(bfac, nullptr);
        atoms->pdbinfo[natom].occup = strtod(occup, nullptr);
    }
    natom++;

    return natom;
}

gmx_bool is_hydrogen(const char* nm)
{
    char buf[30];

    std::strcpy(buf, nm);
    trim(buf);

    return ((buf[0] == 'H') || ((std::isdigit(buf[0]) != 0) && (buf[1] == 'H')));
}

gmx_bool is_dummymass(const char* nm)
{
    char buf[30];

    std::strcpy(buf, nm);
    trim(buf);

    return (buf[0] == 'M') && (std::isdigit(buf[strlen(buf) - 1]) != 0);
}

static void gmx_conect_addline(gmx_conect_t* con, char* line)
{
    int n, ai, aj;

    std::string form2  = "%*s";
    std::string format = form2 + "%d";
    if (sscanf(line, format.c_str(), &ai) == 1)
    {
        do
        {
            form2 += "%*s";
            format = form2 + "%d";
            n      = sscanf(line, format.c_str(), &aj);
            if (n == 1)
            {
                gmx_conect_add(con, ai - 1, aj - 1); /* to prevent duplicated records */
            }
        } while (n == 1);
    }
}

void gmx_conect_dump(FILE* fp, gmx_conect conect)
{
    gmx_conect_t* gc = static_cast<gmx_conect_t*>(conect);
    int           i;

    for (i = 0; (i < gc->nconect); i++)
    {
        fprintf(fp, "%6s%5d%5d\n", "CONECT", gc->conect[i].ai + 1, gc->conect[i].aj + 1);
    }
}

gmx_conect gmx_conect_init()
{
    gmx_conect_t* gc;

    snew(gc, 1);

    return gc;
}

void gmx_conect_done(gmx_conect conect)
{
    gmx_conect_t* gc = conect;

    sfree(gc->conect);
}

gmx_bool gmx_conect_exist(gmx_conect conect, int ai, int aj)
{
    gmx_conect_t* gc = conect;
    int           i;

    for (i = 0; (i < gc->nconect); i++)
    {
        if (((gc->conect[i].ai == ai) && (gc->conect[i].aj == aj))
            || ((gc->conect[i].aj == ai) && (gc->conect[i].ai == aj)))
        {
            return TRUE;
        }
    }
    return FALSE;
}

void gmx_conect_add(gmx_conect conect, int ai, int aj)
{
    gmx_conect_t* gc = static_cast<gmx_conect_t*>(conect);

    if (!gmx_conect_exist(conect, ai, aj))
    {
        srenew(gc->conect, ++gc->nconect);
        gc->conect[gc->nconect - 1].ai = ai;
        gc->conect[gc->nconect - 1].aj = aj;
    }
}

int read_pdbfile(FILE*      in,
                 char*      title,
                 int*       model_nr,
                 t_atoms*   atoms,
                 t_symtab*  symtab,
                 rvec       x[],
                 PbcType*   pbcType,
                 matrix     box,
                 gmx_conect conect)
{
    gmx_conect_t* gc = conect;
    gmx_bool      bCOMPND;
    gmx_bool      bConnWarn = FALSE;
    char          line[STRLEN + 1];
    char *        c, *d;
    int           natom, chainnum;
    gmx_bool      bStop = FALSE;

    if (pbcType)
    {
        /* Only assume pbc when there is a CRYST1 entry */
        *pbcType = PbcType::No;
    }
    if (box != nullptr)
    {
        clear_mat(box);
    }

    atoms->haveMass    = FALSE;
    atoms->haveCharge  = FALSE;
    atoms->haveType    = FALSE;
    atoms->haveBState  = FALSE;
    atoms->havePdbInfo = (atoms->pdbinfo != nullptr);

    bCOMPND  = FALSE;
    title[0] = '\0';
    natom    = 0;
    chainnum = 0;
    static const gmx::StringToEnumValueConverter<PdbRecordType, enumValueToString, gmx::StringCompareType::Exact, gmx::StripStrings::Yes>
            sc_pdbRecordTypeIdentifier;
    while (!bStop && (fgets2(line, STRLEN, in) != nullptr))
    {
        // Convert line to a null-terminated string, then take a substring of at most 6 chars
        std::string                  lineAsString(line);
        std::optional<PdbRecordType> line_type =
                sc_pdbRecordTypeIdentifier.valueFrom(lineAsString.substr(0, 6));
        if (!line_type)
        {
            // Skip this line because it does not start with a
            // recognized leading string describing a PDB record type.
            continue;
        }

        switch (line_type.value())
        {
            case PdbRecordType::Atom:
            case PdbRecordType::Hetatm:
                natom = read_atom(symtab, line, line_type.value(), natom, atoms, x, chainnum);
                break;

            case PdbRecordType::Anisou:
                if (atoms->havePdbInfo)
                {
                    read_anisou(line, natom, atoms);
                }
                break;

            case PdbRecordType::Cryst1: read_cryst1(line, pbcType, box); break;

            case PdbRecordType::Title:
            case PdbRecordType::Header:
                if (std::strlen(line) > 6)
                {
                    c = line + 6;
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
                    d = std::strstr(c, "      ");
                    if (d)
                    {
                        d[0] = '\0';
                    }
                    if (std::strlen(c) > 0)
                    {
                        std::strcpy(title, c);
                    }
                }
                break;

            case PdbRecordType::Compound:
                if ((!std::strstr(line, ": ")) || (std::strstr(line + 6, "MOLECULE:")))
                {
                    if (!(c = std::strstr(line + 6, "MOLECULE:")))
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
                        while ((d[-1] == ';') && d > c)
                        {
                            d--;
                        }
                        d[0] = '\0';
                    }
                    if (strlen(c) > 0)
                    {
                        if (bCOMPND)
                        {
                            std::strcat(title, "; ");
                            std::strcat(title, c);
                        }
                        else
                        {
                            std::strcpy(title, c);
                        }
                    }
                    bCOMPND = TRUE;
                }
                break;

            case PdbRecordType::Ter: chainnum++; break;

            case PdbRecordType::Model:
                if (model_nr)
                {
                    sscanf(line, "%*s%d", model_nr);
                }
                break;

            case PdbRecordType::EndModel: bStop = TRUE; break;
            case PdbRecordType::Conect:
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

            default: break;
        }
    }

    return natom;
}

void get_pdb_coordnum(FILE* in, int* natoms)
{
    char line[STRLEN];

    *natoms = 0;
    while (fgets2(line, STRLEN, in))
    {
        if (std::strncmp(line, "ENDMDL", 6) == 0)
        {
            break;
        }
        if ((std::strncmp(line, "ATOM  ", 6) == 0) || (std::strncmp(line, "HETATM", 6) == 0))
        {
            (*natoms)++;
        }
    }
}

void gmx_pdb_read_conf(const std::filesystem::path& infile,
                       t_symtab*                    symtab,
                       char**                       name,
                       t_atoms*                     atoms,
                       rvec                         x[],
                       PbcType*                     pbcType,
                       matrix                       box)
{
    FILE* in = gmx_fio_fopen(infile, "r");
    char  title[STRLEN];
    read_pdbfile(in, title, nullptr, atoms, symtab, x, pbcType, box, nullptr);
    if (name != nullptr)
    {
        *name = gmx_strdup(title);
    }
    gmx_fio_fclose(in);
}

gmx_conect gmx_conect_generate(const t_topology* top)
{
    int        f, i;
    gmx_conect gc;

    /* Fill the conect records */
    gc = gmx_conect_init();

    for (f = 0; (f < F_NRE); f++)
    {
        if (IS_CHEMBOND(f))
        {
            for (i = 0; (i < top->idef.il[f].nr); i += interaction_function[f].nratoms + 1)
            {
                gmx_conect_add(gc, top->idef.il[f].iatoms[i + 1], top->idef.il[f].iatoms[i + 2]);
            }
        }
    }
    return gc;
}

int gmx_fprintf_pdb_atomline(FILE*         fp,
                             PdbRecordType record,
                             int           atom_seq_number,
                             const char*   atom_name,
                             char          alternate_location,
                             const char*   res_name,
                             char          chain_id,
                             int           res_seq_number,
                             char          res_insertion_code,
                             real          x,
                             real          y,
                             real          z,
                             real          occupancy,
                             real          b_factor,
                             const char*   element)
{
    char     tmp_atomname[6], tmp_resname[6];
    gmx_bool start_name_in_col13;
    int      n;

    if (record != PdbRecordType::Atom && record != PdbRecordType::Hetatm)
    {
        gmx_fatal(FARGS, "Can only print PDB atom lines as ATOM or HETATM records");
    }

    /* Format atom name */
    if (atom_name != nullptr)
    {
        /* If the atom name is an element name with two chars, it should start already in column 13.
         * Otherwise it should start in column 14, unless the name length is 4 chars.
         */
        if ((element != nullptr) && (std::strlen(element) >= 2)
            && (gmx_strncasecmp(atom_name, element, 2) == 0))
        {
            start_name_in_col13 = TRUE;
        }
        else
        {
            start_name_in_col13 = (std::strlen(atom_name) >= 4);
        }
        snprintf(tmp_atomname, sizeof(tmp_atomname), start_name_in_col13 ? "" : " ");
        std::strncat(tmp_atomname, atom_name, 4);
        tmp_atomname[5] = '\0';
    }
    else
    {
        tmp_atomname[0] = '\0';
    }

    /* Format residue name */
    std::strncpy(tmp_resname, (res_name != nullptr) ? res_name : "", 4);
    /* Make sure the string is terminated if strlen was > 4 */
    tmp_resname[4] = '\0';
    /* String is properly terminated, so now we can use strcat. By adding a
     * space we can write it right-justified, and if the original name was
     * three characters or less there will be a space added on the right side.
     */
    std::strcat(tmp_resname, " ");

    /* Truncate integers so they fit */
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    n = fprintf(fp,
                "%-6s%5d %-4.4s%c%4.4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                enumValueToString(record),
                atom_seq_number,
                tmp_atomname,
                alternate_location,
                tmp_resname,
                chain_id,
                res_seq_number,
                res_insertion_code,
                x,
                y,
                z,
                occupancy,
                b_factor,
                (element != nullptr) ? element : "");

    return n;
}
