/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
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
            write_hconf_indexed_p(out, title, atoms, nindex, index, 3, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.bTitle = TRUE;
            fr.title  = title;
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
            write_g96_conf(out, &fr, nindex, index);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
        case efPQR:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile_indexed(out, title, atoms, x, ePBC, box, ' ', -1, nindex, index, NULL, TRUE);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, nindex, index, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
            break;
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
            write_conf_p(outfile, title, atoms, 3, x, v, box);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.bTitle = TRUE;
            fr.title  = title;
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
            write_g96_conf(out, &fr, -1, NULL);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile(out, title, atoms, x, ePBC, box, ' ', -1, NULL, TRUE);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, atoms->nr, NULL, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
            break;
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
            write_hconf_mtop(out, title, mtop, 3, x, v, box);
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
            fr.title  = NULL;
            fr.natoms = -1;
            fr.atoms  = NULL;
            fr.x      = NULL;
            fr.v      = NULL;
            fr.f      = NULL;
            *natoms   = read_g96_conf(in, infile, &fr, NULL, g96_line);
            sfree(const_cast<char *>(fr.title));
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

static void read_stx_conf(const char *infile, t_topology *top,
                          rvec x[], rvec *v, int *ePBC, matrix box)
{
    FILE       *in;
    t_trxframe  fr;
    int         ftp;
    char        g96_line[STRLEN+1];

    if (top->atoms.nr == 0)
    {
        fprintf(stderr, "Warning: Number of atoms in %s is 0\n", infile);
    }
    else if (top->atoms.atom == NULL)
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
            gmx_gro_read_conf(infile, top, x, v, box);
            break;
        case efG96:
            fr.title  = NULL;
            fr.natoms = top->atoms.nr;
            fr.atoms  = &top->atoms;
            fr.x      = x;
            fr.v      = v;
            fr.f      = NULL;
            in        = gmx_fio_fopen(infile, "r");
            read_g96_conf(in, infile, &fr, &top->symtab, g96_line);
            gmx_fio_fclose(in);
            copy_mat(fr.box, box);
            top->name = put_symtab(&top->symtab, fr.title);
            sfree(const_cast<char *>(fr.title));
            break;
        case efPDB:
        case efBRK:
        case efENT:
            gmx_pdb_read_conf(infile, top, x, ePBC, box);
            break;
        case efESP:
            gmx_espresso_read_conf(infile, top, x, v, box);
            break;
        default:
            gmx_incons("Not supported in read_stx_conf");
    }
}

static void done_gmx_groups_t(gmx_groups_t *g)
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        if (NULL != g->grps[i].nm_ind)
        {
            sfree(g->grps[i].nm_ind);
            g->grps[i].nm_ind = NULL;
        }
        if (NULL != g->grpnr[i])
        {
            sfree(g->grpnr[i]);
            g->grpnr[i] = NULL;
        }
    }
    /* The contents of this array is in symtab, don't free it here */
    sfree(g->grpname);
}

gmx_bool read_tps_conf(const char *infile, t_topology *top, int *ePBC,
                       rvec **x, rvec **v, matrix box, gmx_bool bMass)
{
    t_tpxheader      header;
    int              natoms, i;
    gmx_bool         bTop, bXNULL = FALSE;
    gmx_mtop_t      *mtop;
    gmx_atomprop_t   aps;

    bTop  = fn2bTPX(infile);
    if (ePBC != NULL)
    {
        *ePBC = -1;
    }
    if (bTop)
    {
        read_tpxheader(infile, &header, TRUE);
        if (x)
        {
            snew(*x, header.natoms);
        }
        if (v)
        {
            snew(*v, header.natoms);
        }
        snew(mtop, 1);
        int ePBC_tmp
            = read_tpx(infile, NULL, box, &natoms,
                       (x == NULL) ? NULL : *x, (v == NULL) ? NULL : *v, mtop);
        if (ePBC != NULL)
        {
            *ePBC = ePBC_tmp;
        }
        *top = gmx_mtop_t_to_t_topology(mtop);
        /* In this case we need to throw away the group data too */
        done_gmx_groups_t(&mtop->groups);
        sfree(mtop);
        tpx_make_chain_identifiers(&top->atoms, &top->mols);
    }
    else
    {
        open_symtab(&top->symtab);
        get_stx_coordnum(infile, &natoms);
        init_t_atoms(&top->atoms, natoms, (fn2ftp(infile) == efPDB));
        if (x == NULL)
        {
            snew(x, 1);
            bXNULL = TRUE;
        }
        snew(*x, natoms);
        if (v)
        {
            snew(*v, natoms);
        }
        read_stx_conf(infile, top, *x, (v == NULL) ? NULL : *v, ePBC, box);
        if (bXNULL)
        {
            sfree(*x);
            sfree(x);
        }
        if (bMass)
        {
            aps = gmx_atomprop_init();
            for (i = 0; (i < natoms); i++)
            {
                if (!gmx_atomprop_query(aps, epropMass,
                                        *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                        *top->atoms.atomname[i],
                                        &(top->atoms.atom[i].m)))
                {
                    if (debug)
                    {
                        fprintf(debug, "Can not find mass for atom %s %d %s, setting to 1\n",
                                *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                top->atoms.resinfo[top->atoms.atom[i].resind].nr,
                                *top->atoms.atomname[i]);
                    }
                }
            }
            gmx_atomprop_destroy(aps);
        }
        top->idef.ntypes = -1;
    }

    return bTop;
}
