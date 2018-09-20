/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "confio.h"

#include <cstdio>
#include <cstring>

#include "gromacs/fileio/espio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

void write_sto_conf_indexed(const char *outfile, const char *title,
                            const t_atoms *atoms,
                            const rvec x[], const rvec *v, int ePBC, const matrix box,
                            int nindex, int index[])
{
    FILE       *out;
    int         ftp;
    t_trxframe  fr;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            out = gmx_fio_fopen(outfile, "w");
            write_hconf_indexed_p(out, title, atoms, nindex, index, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.natoms = atoms->nr;
            fr.bAtoms = TRUE;
            fr.atoms  = const_cast<t_atoms *>(atoms);
            fr.bX     = TRUE;
            fr.x      = const_cast<rvec *>(x);
            if (v)
            {
                fr.bV = TRUE;
                fr.v  = const_cast<rvec *>(v);
            }
            fr.bBox = TRUE;
            copy_mat(box, fr.box);
            out = gmx_fio_fopen(outfile, "w");
            write_g96_conf(out, title, &fr, nindex, index);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
        case efPQR:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile_indexed(out, title, atoms, x, ePBC, box, ' ', -1, nindex, index, nullptr, TRUE, ftp == efPQR);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, nindex, index, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
        default:
            gmx_incons("Not supported in write_sto_conf_indexed");
    }
}

void write_sto_conf(const char *outfile, const char *title, const t_atoms *atoms,
                    const rvec x[], const rvec *v, int ePBC, const matrix box)
{
    FILE       *out;
    int         ftp;
    t_trxframe  fr;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            write_conf_p(outfile, title, atoms, x, v, box);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.natoms = atoms->nr;
            fr.bAtoms = TRUE;
            fr.atoms  = const_cast<t_atoms *>(atoms); // TODO check
            fr.bX     = TRUE;
            fr.x      = const_cast<rvec *>(x);
            if (v)
            {
                fr.bV = TRUE;
                fr.v  = const_cast<rvec *>(v);
            }
            fr.bBox = TRUE;
            copy_mat(box, fr.box);
            out = gmx_fio_fopen(outfile, "w");
            write_g96_conf(out, title, &fr, -1, nullptr);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile(out, title, atoms, x, ePBC, box, ' ', -1, nullptr, TRUE);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, atoms->nr, nullptr, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
        default:
            gmx_incons("Not supported in write_sto_conf");
    }
}

void write_sto_conf_mtop(const char *outfile, const char *title,
                         gmx_mtop_t *mtop,
                         const rvec x[], const rvec *v, int ePBC, const matrix box)
{
    int     ftp;
    FILE   *out;
    t_atoms atoms;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            out = gmx_fio_fopen(outfile, "w");
            write_hconf_mtop(out, title, mtop, x, v, box);
            gmx_fio_fclose(out);
            break;
        default:
            /* This is a brute force approach which requires a lot of memory.
             * We should implement mtop versions of all writing routines.
             */
            atoms = gmx_mtop_global_atoms(mtop);

            write_sto_conf(outfile, title, &atoms, x, v, ePBC, box);

            done_atom(&atoms);
            break;
    }
}

static void get_stx_coordnum(const char *infile, int *natoms)
{
    FILE      *in;
    int        ftp;
    t_trxframe fr;
    char       g96_line[STRLEN+1];

    ftp = fn2ftp(infile);
    range_check(ftp, 0, efNR);
    switch (ftp)
    {
        case efGRO:
            get_coordnum(infile, natoms);
            break;
        case efG96:
        {
            in        = gmx_fio_fopen(infile, "r");
            fr.natoms = -1;
            fr.atoms  = nullptr;
            fr.x      = nullptr;
            fr.v      = nullptr;
            fr.f      = nullptr;
            *natoms   = read_g96_conf(in, infile, nullptr, &fr, nullptr, g96_line);
            gmx_fio_fclose(in);
            break;
        }
        case efPDB:
        case efBRK:
        case efENT:
            in = gmx_fio_fopen(infile, "r");
            get_pdb_coordnum(in, natoms);
            gmx_fio_fclose(in);
            break;
        case efESP:
            *natoms = get_espresso_coordnum(infile);
            break;
        default:
            gmx_fatal(FARGS, "File type %s not supported in get_stx_coordnum",
                      ftp2ext(ftp));
    }
}

// TODO molecule index handling is suspected of being broken here
static void tpx_make_chain_identifiers(t_atoms *atoms, t_block *mols)
{
    /* We always assign a new chain number, but save the chain id characters
     * for larger molecules.
     */
    const int chainMinAtoms = 15;

    int       chainnum = 0;
    char      chainid  = 'A';
    bool      outOfIds = false;
    for (int m = 0; m < mols->nr; m++)
    {
        int a0 = mols->index[m];
        int a1 = mols->index[m+1];
        int c;
        if (a1 - a0 >= chainMinAtoms && !outOfIds)
        {
            /* Set the chain id for the output */
            c = chainid;
            /* Here we allow for the max possible 2*26+10=62 chain ids */
            if (chainid == 'Z')
            {
                chainid = 'a';
            }
            else if (chainid == 'z')
            {
                chainid = '0';
            }
            else if (chainid == '9')
            {
                outOfIds = true;
            }
            else
            {
                chainid++;
            }
        }
        else
        {
            c = ' ';
        }
        for (int a = a0; a < a1; a++)
        {
            atoms->resinfo[atoms->atom[a].resind].chainnum = chainnum;
            atoms->resinfo[atoms->atom[a].resind].chainid  = c;
        }
        chainnum++;
    }

    /* Blank out the chain id if there was only one chain */
    if (chainid == 'B')
    {
        for (int r = 0; r < atoms->nres; r++)
        {
            atoms->resinfo[r].chainid = ' ';
        }
    }
}

static void read_stx_conf(const char *infile,
                          t_symtab *symtab, char **name, t_atoms *atoms,
                          rvec x[], rvec *v, int *ePBC, matrix box)
{
    FILE       *in;
    t_trxframe  fr;
    int         ftp;
    char        g96_line[STRLEN+1];

    if (atoms->nr == 0)
    {
        fprintf(stderr, "Warning: Number of atoms in %s is 0\n", infile);
    }
    else if (atoms->atom == nullptr)
    {
        gmx_mem("Uninitialized array atom");
    }

    if (ePBC)
    {
        *ePBC = -1;
    }

    ftp = fn2ftp(infile);
    switch (ftp)
    {
        case efGRO:
            gmx_gro_read_conf(infile, symtab, name, atoms, x, v, box);
            break;
        case efG96:
            fr.natoms = atoms->nr;
            fr.atoms  = atoms;
            fr.x      = x;
            fr.v      = v;
            fr.f      = nullptr;
            in        = gmx_fio_fopen(infile, "r");
            read_g96_conf(in, infile, name, &fr, symtab, g96_line);
            gmx_fio_fclose(in);
            copy_mat(fr.box, box);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            gmx_pdb_read_conf(infile, symtab, name, atoms, x, ePBC, box);
            break;
        case efESP:
            gmx_espresso_read_conf(infile, symtab, name, atoms, x, v, box);
            break;
        default:
            gmx_incons("Not supported in read_stx_conf");
    }
}

void readConfAndAtoms(const char *infile,
                      t_symtab *symtab, char **name, t_atoms *atoms,
                      int *ePBC,
                      rvec **x, rvec **v, matrix box)
{
    int natoms;
    get_stx_coordnum(infile, &natoms);

    init_t_atoms(atoms, natoms, (fn2ftp(infile) == efPDB));

    bool xIsNull = false;
    if (x == nullptr)
    {
        snew(x, 1);
        xIsNull = true;
    }
    snew(*x, natoms);
    if (v)
    {
        snew(*v, natoms);
    }
    read_stx_conf(infile,
                  symtab, name, atoms,
                  *x, (v == nullptr) ? nullptr : *v, ePBC, box);
    if (xIsNull)
    {
        sfree(*x);
        sfree(x);
    }
}

void readConfAndTopology(const char *infile,
                         bool *haveTopology, gmx_mtop_t *mtop,
                         int *ePBC,
                         rvec **x, rvec **v, matrix box)
{
    GMX_RELEASE_ASSERT(mtop != nullptr, "readConfAndTopology requires mtop!=NULL");

    if (ePBC != nullptr)
    {
        *ePBC = -1;
    }

    *haveTopology = fn2bTPX(infile);
    if (*haveTopology)
    {
        t_tpxheader header;
        read_tpxheader(infile, &header, TRUE);
        if (x)
        {
            snew(*x, header.natoms);
        }
        if (v)
        {
            snew(*v, header.natoms);
        }
        int natoms;
        int ePBC_tmp
            = read_tpx(infile, nullptr, box, &natoms,
                       (x == nullptr) ? nullptr : *x, (v == nullptr) ? nullptr : *v, mtop);
        if (ePBC != nullptr)
        {
            *ePBC = ePBC_tmp;
        }
    }
    else
    {
        t_symtab   symtab;
        char      *name;
        t_atoms    atoms;

        open_symtab(&symtab);

        readConfAndAtoms(infile, &symtab, &name, &atoms, ePBC, x, v, box);

        convertAtomsToMtop(&symtab, put_symtab(&symtab, name), &atoms, mtop);
        sfree(name);
    }
}

gmx_bool read_tps_conf(const char *infile, t_topology *top, int *ePBC,
                       rvec **x, rvec **v, matrix box, gmx_bool requireMasses)
{
    bool        haveTopology;
    gmx_mtop_t  mtop;

    readConfAndTopology(infile, &haveTopology, &mtop, ePBC, x, v, box);

    *top = gmx_mtop_t_to_t_topology(&mtop, true);

    tpx_make_chain_identifiers(&top->atoms, &top->mols);

    if (requireMasses && !top->atoms.haveMass)
    {
        atomsSetMassesBasedOnNames(&top->atoms, TRUE);

        if (!top->atoms.haveMass)
        {
            gmx_fatal(FARGS, "Masses were requested, but for some atom(s) masses could not be found in the database. Use a tpr file as input, if possible, or add these atoms to the mass database.");
        }
    }

    return haveTopology;
}
