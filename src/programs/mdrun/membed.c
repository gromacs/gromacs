/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "membed.h"

#include <signal.h>
#include <stdlib.h>

#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/* information about scaling center */
typedef struct {
    rvec      xmin;       /* smallest coordinates of all embedded molecules */
    rvec      xmax;       /* largest coordinates of all embedded molecules */
    rvec     *geom_cent;  /* scaling center of each independent molecule to embed */
    int       pieces;     /* number of molecules to embed independently */
    int      *nidx;       /* n atoms for every independent embedded molecule (index in subindex) */
    atom_id **subindex;   /* atomids for independent molecule *
                           * atoms of piece i run from subindex[i][0] to subindex[i][nidx[i]] */
} pos_ins_t;

/* variables needed in do_md */
struct membed {
    int        it_xy;     /* number of iterations (steps) used to grow something in the xy-plane */
    int        it_z;      /* same, but for z */
    real       xy_step;   /* stepsize used during resize in xy-plane */
    real       z_step;    /* same, but in z */
    rvec       fac;       /* initial scaling of the molecule to grow into the membrane */
    rvec      *r_ins;     /* final coordinates of the molecule to grow  */
    pos_ins_t *pos_ins;   /* scaling center for each piece to embed */
};

/* membrane related variables */
typedef struct {
    char      *name;     /* name of index group to embed molecule into (usually membrane) */
    t_block    mem_at;   /* list all atoms in membrane */
    int        nmol;     /* number of membrane molecules overlapping with the molecule to embed */
    int       *mol_id;   /* list of molecules in membrane that overlap with the molecule to embed */
    real       lip_area; /* average area per lipid in membrane (only correct for homogeneous bilayers)*/
    real       zmin;     /* minimum z coordinate of membrane */
    real       zmax;     /* maximum z coordinate of membrane */
    real       zmed;     /* median z coordinate of membrane */
} mem_t;

/* Lists all molecules in the membrane that overlap with the molecule to be embedded. *
 * These will then be removed from the system */
typedef struct {
    int   nr;     /* number of molecules to remove */
    int  *mol;    /* list of molecule ids to remove */
    int  *block;  /* id of the molblock that the molecule to remove is part of */
} rm_t;

/* Get the global molecule id, and the corresponding molecule type and id of the *
 * molblock from the global atom nr. */
static int get_mol_id(int at, gmx_mtop_t  *mtop, int *type, int *block)
{
    int                   mol_id = 0;
    int                   i;
    int                   atnr_mol;
    gmx_mtop_atomlookup_t alook;

    alook = gmx_mtop_atomlookup_settle_init(mtop);
    gmx_mtop_atomnr_to_molblock_ind(alook, at, block, &mol_id, &atnr_mol);
    for (i = 0; i < *block; i++)
    {
        mol_id += mtop->molblock[i].nmol;
    }
    *type = mtop->molblock[*block].type;

    gmx_mtop_atomlookup_destroy(alook);

    return mol_id;
}

/* Get the id of the molblock from a global molecule id */
static int get_molblock(int mol_id, int nmblock, gmx_molblock_t *mblock)
{
    int i;
    int nmol = 0;

    for (i = 0; i < nmblock; i++)
    {
        nmol += mblock[i].nmol;
        if (mol_id < nmol)
        {
            return i;
        }
    }

    gmx_fatal(FARGS, "mol_id %d larger than total number of molecules %d.\n", mol_id, nmol);

    return -1;
}

static int get_tpr_version(const char *infile)
{
    t_tpxheader  header;
    int          version, generation;

    read_tpxheader(infile, &header, TRUE, &version, &generation);

    return version;
}

/* Get a list of all the molecule types that are present in a group of atoms. *
 * Because all interaction within the group to embed are removed on the topology *
 * level, if the same molecule type is found in another part of the system, these *
 * would also be affected. Therefore we have to check if the embedded and rest group *
 * share common molecule types. If so, membed will stop with an error. */
static int get_mtype_list(t_block *at, gmx_mtop_t *mtop, t_block *tlist)
{
    int      i, j, nr, mol_id;
    int      type = 0, block = 0;
    gmx_bool bNEW;

    nr = 0;
    snew(tlist->index, at->nr);
    for (i = 0; i < at->nr; i++)
    {
        bNEW   = TRUE;
        mol_id = get_mol_id(at->index[i], mtop, &type, &block);
        for (j = 0; j < nr; j++)
        {
            if (tlist->index[j] == type)
            {
                bNEW = FALSE;
            }
        }

        if (bNEW == TRUE)
        {
            tlist->index[nr] = type;
            nr++;
        }
    }
    srenew(tlist->index, nr);
    return nr;
}

/* Do the actual check of the molecule types between embedded and rest group */
static void check_types(t_block *ins_at, t_block *rest_at, gmx_mtop_t *mtop)
{
    t_block        *ins_mtype, *rest_mtype;
    int             i, j;

    snew(ins_mtype, 1);
    snew(rest_mtype, 1);
    ins_mtype->nr  = get_mtype_list(ins_at, mtop, ins_mtype );
    rest_mtype->nr = get_mtype_list(rest_at, mtop, rest_mtype);

    for (i = 0; i < ins_mtype->nr; i++)
    {
        for (j = 0; j < rest_mtype->nr; j++)
        {
            if (ins_mtype->index[i] == rest_mtype->index[j])
            {
                gmx_fatal(FARGS, "Moleculetype %s is found both in the group to insert and the rest of the system.\n"
                          "1. Your *.ndx and *.top do not match\n"
                          "2. You are inserting some molecules of type %s (for example xray-solvent), while\n"
                          "the same moleculetype is also used in the rest of the system (solvent box). Because\n"
                          "we need to exclude all interactions between the atoms in the group to\n"
                          "insert, the same moleculetype can not be used in both groups. Change the\n"
                          "moleculetype of the molecules %s in the inserted group. Do not forget to provide\n"
                          "an appropriate *.itp file", *(mtop->moltype[rest_mtype->index[j]].name),
                          *(mtop->moltype[rest_mtype->index[j]].name), *(mtop->moltype[rest_mtype->index[j]].name));
            }
        }
    }

    sfree(ins_mtype->index);
    sfree(rest_mtype->index);
    sfree(ins_mtype);
    sfree(rest_mtype);
}

static void get_input(const char *membed_input, real *xy_fac, real *xy_max, real *z_fac, real *z_max,
                      int *it_xy, int *it_z, real *probe_rad, int *low_up_rm, int *maxwarn,
                      int *pieces, gmx_bool *bALLOW_ASYMMETRY)
{
    warninp_t  wi;
    t_inpfile *inp;
    int        ninp;

    wi = init_warning(TRUE, 0);

    inp = read_inpfile(membed_input, &ninp, wi);
    ITYPE ("nxy", *it_xy, 1000);
    ITYPE ("nz", *it_z, 0);
    RTYPE ("xyinit", *xy_fac, 0.5);
    RTYPE ("xyend", *xy_max, 1.0);
    RTYPE ("zinit", *z_fac, 1.0);
    RTYPE ("zend", *z_max, 1.0);
    RTYPE ("rad", *probe_rad, 0.22);
    ITYPE ("ndiff", *low_up_rm, 0);
    ITYPE ("maxwarn", *maxwarn, 0);
    ITYPE ("pieces", *pieces, 1);
    EETYPE("asymmetry", *bALLOW_ASYMMETRY, yesno_names);

    write_inpfile(membed_input, ninp, inp, FALSE, wi);
}

/* Obtain the maximum and minimum coordinates of the group to be embedded */
static int init_ins_at(t_block *ins_at, t_block *rest_at, t_state *state, pos_ins_t *pos_ins,
                       gmx_groups_t *groups, int ins_grp_id, real xy_max)
{
    int        i, gid, c = 0;
    real       x, xmin, xmax, y, ymin, ymax, z, zmin, zmax;
    const real min_memthick = 6.0;      /* minimum thickness of the bilayer that will be used to *
                                         * determine the overlap between molecule to embed and membrane */
    const real fac_inp_size = 1.000001; /* scaling factor to obtain input_size + 0.000001 (comparing reals) */
    snew(rest_at->index, state->natoms);

    xmin = xmax = state->x[ins_at->index[0]][XX];
    ymin = ymax = state->x[ins_at->index[0]][YY];
    zmin = zmax = state->x[ins_at->index[0]][ZZ];

    for (i = 0; i < state->natoms; i++)
    {
        gid = groups->grpnr[egcFREEZE][i];
        if (groups->grps[egcFREEZE].nm_ind[gid] == ins_grp_id)
        {
            x = state->x[i][XX];
            if (x < xmin)
            {
                xmin = x;
            }
            if (x > xmax)
            {
                xmax = x;
            }
            y = state->x[i][YY];
            if (y < ymin)
            {
                ymin = y;
            }
            if (y > ymax)
            {
                ymax = y;
            }
            z = state->x[i][ZZ];
            if (z < zmin)
            {
                zmin = z;
            }
            if (z > zmax)
            {
                zmax = z;
            }
        }
        else
        {
            rest_at->index[c] = i;
            c++;
        }
    }

    rest_at->nr = c;
    srenew(rest_at->index, c);

    if (xy_max > fac_inp_size)
    {
        pos_ins->xmin[XX] = xmin-((xmax-xmin)*xy_max-(xmax-xmin))/2;
        pos_ins->xmin[YY] = ymin-((ymax-ymin)*xy_max-(ymax-ymin))/2;

        pos_ins->xmax[XX] = xmax+((xmax-xmin)*xy_max-(xmax-xmin))/2;
        pos_ins->xmax[YY] = ymax+((ymax-ymin)*xy_max-(ymax-ymin))/2;
    }
    else
    {
        pos_ins->xmin[XX] = xmin;
        pos_ins->xmin[YY] = ymin;

        pos_ins->xmax[XX] = xmax;
        pos_ins->xmax[YY] = ymax;
    }

    if ( (zmax-zmin) < min_memthick)
    {
        pos_ins->xmin[ZZ] = zmin+(zmax-zmin)/2.0-0.5*min_memthick;
        pos_ins->xmax[ZZ] = zmin+(zmax-zmin)/2.0+0.5*min_memthick;
    }
    else
    {
        pos_ins->xmin[ZZ] = zmin;
        pos_ins->xmax[ZZ] = zmax;
    }

    return c;
}

/* Estimate the area of the embedded molecule by projecting all coordinates on a grid in the *
 * xy-plane and counting the number of occupied grid points */
static real est_prot_area(pos_ins_t *pos_ins, rvec *r, t_block *ins_at, mem_t *mem_p)
{
    real x, y, dx = 0.15, dy = 0.15, area = 0.0;
    real add, memmin, memmax;
    int  c, at;

    /* min and max membrane coordinate are altered to reduce the influence of the *
     * boundary region */
    memmin = mem_p->zmin+0.1*(mem_p->zmax-mem_p->zmin);
    memmax = mem_p->zmax-0.1*(mem_p->zmax-mem_p->zmin);

    for (x = pos_ins->xmin[XX]; x < pos_ins->xmax[XX]; x += dx)
    {
        for (y = pos_ins->xmin[YY]; y < pos_ins->xmax[YY]; y += dy)
        {
            c   = 0;
            add = 0.0;
            do
            {
                at = ins_at->index[c];
                if ( (r[at][XX] >= x) && (r[at][XX] < x+dx) &&
                     (r[at][YY] >= y) && (r[at][YY] < y+dy) &&
                     (r[at][ZZ] > memmin) && (r[at][ZZ] < memmax) )
                {
                    add = 1.0;
                }
                c++;
            }
            while ( (c < ins_at->nr) && (add < 0.5) );
            area += add;
        }
    }
    area = area*dx*dy;

    return area;
}

static int init_mem_at(mem_t *mem_p, gmx_mtop_t *mtop, rvec *r, matrix box, pos_ins_t *pos_ins)
{
    int      i, j, at, mol, nmol, nmolbox, count;
    t_block *mem_a;
    real     z, zmin, zmax, mem_area;
    gmx_bool bNew;
    atom_id *mol_id;
    int      type = 0, block = 0;

    nmol  = count = 0;
    mem_a = &(mem_p->mem_at);
    snew(mol_id, mem_a->nr);
    zmin = pos_ins->xmax[ZZ];
    zmax = pos_ins->xmin[ZZ];
    for (i = 0; i < mem_a->nr; i++)
    {
        at = mem_a->index[i];
        if ( (r[at][XX] > pos_ins->xmin[XX]) && (r[at][XX] < pos_ins->xmax[XX]) &&
             (r[at][YY] > pos_ins->xmin[YY]) && (r[at][YY] < pos_ins->xmax[YY]) &&
             (r[at][ZZ] > pos_ins->xmin[ZZ]) && (r[at][ZZ] < pos_ins->xmax[ZZ]) )
        {
            mol  = get_mol_id(at, mtop, &type, &block);
            bNew = TRUE;
            for (j = 0; j < nmol; j++)
            {
                if (mol == mol_id[j])
                {
                    bNew = FALSE;
                }
            }

            if (bNew)
            {
                mol_id[nmol] = mol;
                nmol++;
            }

            z = r[at][ZZ];
            if (z < zmin)
            {
                zmin = z;
            }

            if (z > zmax)
            {
                zmax = z;
            }

            count++;
        }
    }

    mem_p->nmol = nmol;
    srenew(mol_id, nmol);
    mem_p->mol_id = mol_id;

    if ((zmax-zmin) > (box[ZZ][ZZ]-0.5))
    {
        gmx_fatal(FARGS, "Something is wrong with your membrane. Max and min z values are %f and %f.\n"
                  "Maybe your membrane is not centered in the box, but located at the box edge in the z-direction,\n"
                  "so that one membrane is distributed over two periodic box images. Another possibility is that\n"
                  "your water layer is not thick enough.\n", zmax, zmin);
    }
    mem_p->zmin = zmin;
    mem_p->zmax = zmax;
    mem_p->zmed = (zmax-zmin)/2+zmin;

    /*number of membrane molecules in protein box*/
    nmolbox = count/mtop->molblock[block].natoms_mol;
    /*membrane area within the box defined by the min and max coordinates of the embedded molecule*/
    mem_area = (pos_ins->xmax[XX]-pos_ins->xmin[XX])*(pos_ins->xmax[YY]-pos_ins->xmin[YY]);
    /*rough estimate of area per lipid, assuming there is only one type of lipid in the membrane*/
    mem_p->lip_area = 2.0*mem_area/(double)nmolbox;

    return mem_p->mem_at.nr;
}

static void init_resize(t_block *ins_at, rvec *r_ins, pos_ins_t *pos_ins, mem_t *mem_p, rvec *r,
                        gmx_bool bALLOW_ASYMMETRY)
{
    int i, j, at, c, outsidesum, gctr = 0;
    int idxsum = 0;

    /*sanity check*/
    for (i = 0; i < pos_ins->pieces; i++)
    {
        idxsum += pos_ins->nidx[i];
    }

    if (idxsum != ins_at->nr)
    {
        gmx_fatal(FARGS, "Piecewise sum of inserted atoms not same as size of group selected to insert.");
    }

    snew(pos_ins->geom_cent, pos_ins->pieces);
    for (i = 0; i < pos_ins->pieces; i++)
    {
        c          = 0;
        outsidesum = 0;
        for (j = 0; j < DIM; j++)
        {
            pos_ins->geom_cent[i][j] = 0;
        }

        for (j = 0; j < pos_ins->nidx[i]; j++)
        {
            at = pos_ins->subindex[i][j];
            copy_rvec(r[at], r_ins[gctr]);
            if ( (r_ins[gctr][ZZ] < mem_p->zmax) && (r_ins[gctr][ZZ] > mem_p->zmin) )
            {
                rvec_inc(pos_ins->geom_cent[i], r_ins[gctr]);
                c++;
            }
            else
            {
                outsidesum++;
            }
            gctr++;
        }

        if (c > 0)
        {
            svmul(1/(double)c, pos_ins->geom_cent[i], pos_ins->geom_cent[i]);
        }

        if (!bALLOW_ASYMMETRY)
        {
            pos_ins->geom_cent[i][ZZ] = mem_p->zmed;
        }

        fprintf(stderr, "Embedding piece %d with center of geometry: %f %f %f\n",
                i, pos_ins->geom_cent[i][XX], pos_ins->geom_cent[i][YY], pos_ins->geom_cent[i][ZZ]);
    }
    fprintf(stderr, "\n");
}

/* resize performed in the md loop */
static void resize(rvec *r_ins, rvec *r, pos_ins_t *pos_ins, rvec fac)
{
    int i, j, k, at, c = 0;
    for (k = 0; k < pos_ins->pieces; k++)
    {
        for (i = 0; i < pos_ins->nidx[k]; i++)
        {
            at = pos_ins->subindex[k][i];
            for (j = 0; j < DIM; j++)
            {
                r[at][j] = pos_ins->geom_cent[k][j]+fac[j]*(r_ins[c][j]-pos_ins->geom_cent[k][j]);
            }
            c++;
        }
    }
}

/* generate the list of membrane molecules that overlap with the molecule to be embedded. *
 * The molecule to be embedded is already reduced in size. */
static int gen_rm_list(rm_t *rm_p, t_block *ins_at, t_block *rest_at, t_pbc *pbc, gmx_mtop_t *mtop,
                       rvec *r, mem_t *mem_p, pos_ins_t *pos_ins, real probe_rad,
                       int low_up_rm, gmx_bool bALLOW_ASYMMETRY)
{
    int      i, j, k, l, at, at2, mol_id;
    int      type = 0, block = 0;
    int      nrm, nupper, nlower;
    real     r_min_rad, z_lip, min_norm;
    gmx_bool bRM;
    rvec     dr, dr_tmp;
    real    *dist;
    int     *order;

    r_min_rad = probe_rad*probe_rad;
    snew(rm_p->mol, mtop->mols.nr);
    snew(rm_p->block, mtop->mols.nr);
    nrm    = nupper = 0;
    nlower = low_up_rm;
    for (i = 0; i < ins_at->nr; i++)
    {
        at = ins_at->index[i];
        for (j = 0; j < rest_at->nr; j++)
        {
            at2 = rest_at->index[j];
            pbc_dx(pbc, r[at], r[at2], dr);

            if (norm2(dr) < r_min_rad)
            {
                mol_id = get_mol_id(at2, mtop, &type, &block);
                bRM    = TRUE;
                for (l = 0; l < nrm; l++)
                {
                    if (rm_p->mol[l] == mol_id)
                    {
                        bRM = FALSE;
                    }
                }

                if (bRM)
                {
                    rm_p->mol[nrm]   = mol_id;
                    rm_p->block[nrm] = block;
                    nrm++;
                    z_lip = 0.0;
                    for (l = 0; l < mem_p->nmol; l++)
                    {
                        if (mol_id == mem_p->mol_id[l])
                        {
                            for (k = mtop->mols.index[mol_id]; k < mtop->mols.index[mol_id+1]; k++)
                            {
                                z_lip += r[k][ZZ];
                            }
                            z_lip /= mtop->molblock[block].natoms_mol;
                            if (z_lip < mem_p->zmed)
                            {
                                nlower++;
                            }
                            else
                            {
                                nupper++;
                            }
                        }
                    }
                }
            }
        }
    }

    /*make sure equal number of lipids from upper and lower layer are removed */
    if ( (nupper != nlower) && (!bALLOW_ASYMMETRY) )
    {
        snew(dist, mem_p->nmol);
        snew(order, mem_p->nmol);
        for (i = 0; i < mem_p->nmol; i++)
        {
            at = mtop->mols.index[mem_p->mol_id[i]];
            pbc_dx(pbc, r[at], pos_ins->geom_cent[0], dr);
            if (pos_ins->pieces > 1)
            {
                /*minimum dr value*/
                min_norm = norm2(dr);
                for (k = 1; k < pos_ins->pieces; k++)
                {
                    pbc_dx(pbc, r[at], pos_ins->geom_cent[k], dr_tmp);
                    if (norm2(dr_tmp) < min_norm)
                    {
                        min_norm = norm2(dr_tmp);
                        copy_rvec(dr_tmp, dr);
                    }
                }
            }
            dist[i] = dr[XX]*dr[XX]+dr[YY]*dr[YY];
            j       = i-1;
            while (j >= 0 && dist[i] < dist[order[j]])
            {
                order[j+1] = order[j];
                j--;
            }
            order[j+1] = i;
        }

        i = 0;
        while (nupper != nlower)
        {
            mol_id = mem_p->mol_id[order[i]];
            block  = get_molblock(mol_id, mtop->nmolblock, mtop->molblock);
            bRM    = TRUE;
            for (l = 0; l < nrm; l++)
            {
                if (rm_p->mol[l] == mol_id)
                {
                    bRM = FALSE;
                }
            }

            if (bRM)
            {
                z_lip = 0;
                for (k = mtop->mols.index[mol_id]; k < mtop->mols.index[mol_id+1]; k++)
                {
                    z_lip += r[k][ZZ];
                }
                z_lip /= mtop->molblock[block].natoms_mol;
                if (nupper > nlower && z_lip < mem_p->zmed)
                {
                    rm_p->mol[nrm]   = mol_id;
                    rm_p->block[nrm] = block;
                    nrm++;
                    nlower++;
                }
                else if (nupper < nlower && z_lip > mem_p->zmed)
                {
                    rm_p->mol[nrm]   = mol_id;
                    rm_p->block[nrm] = block;
                    nrm++;
                    nupper++;
                }
            }
            i++;

            if (i > mem_p->nmol)
            {
                gmx_fatal(FARGS, "Trying to remove more lipid molecules than there are in the membrane");
            }
        }
        sfree(dist);
        sfree(order);
    }

    rm_p->nr = nrm;
    srenew(rm_p->mol, nrm);
    srenew(rm_p->block, nrm);

    return nupper+nlower;
}

/*remove all lipids and waters overlapping and update all important structures (e.g. state and mtop)*/
static void rm_group(gmx_groups_t *groups, gmx_mtop_t *mtop, rm_t *rm_p, t_state *state,
                     t_block *ins_at, pos_ins_t *pos_ins)
{
    int             i, j, k, n, rm, mol_id, at, block;
    rvec           *x_tmp, *v_tmp;
    atom_id        *list, *new_mols;
    unsigned char  *new_egrp[egcNR];
    gmx_bool        bRM;
    int             RMmolblock;

    snew(list, state->natoms);
    n = 0;
    for (i = 0; i < rm_p->nr; i++)
    {
        mol_id = rm_p->mol[i];
        at     = mtop->mols.index[mol_id];
        block  = rm_p->block[i];
        mtop->molblock[block].nmol--;
        for (j = 0; j < mtop->molblock[block].natoms_mol; j++)
        {
            list[n] = at+j;
            n++;
        }
        mtop->mols.index[mol_id] = -1;
    }

    mtop->mols.nr           -= rm_p->nr;
    mtop->mols.nalloc_index -= rm_p->nr;
    snew(new_mols, mtop->mols.nr);
    for (i = 0; i < mtop->mols.nr+rm_p->nr; i++)
    {
        j = 0;
        if (mtop->mols.index[i] != -1)
        {
            new_mols[j] = mtop->mols.index[i];
            j++;
        }
    }
    sfree(mtop->mols.index);
    mtop->mols.index = new_mols;
    mtop->natoms    -= n;
    state->natoms   -= n;
    state->nalloc    = state->natoms;
    snew(x_tmp, state->nalloc);
    snew(v_tmp, state->nalloc);

    for (i = 0; i < egcNR; i++)
    {
        if (groups->grpnr[i] != NULL)
        {
            groups->ngrpnr[i] = state->natoms;
            snew(new_egrp[i], state->natoms);
        }
    }

    rm = 0;
    for (i = 0; i < state->natoms+n; i++)
    {
        bRM = FALSE;
        for (j = 0; j < n; j++)
        {
            if (i == list[j])
            {
                bRM = TRUE;
                rm++;
            }
        }

        if (!bRM)
        {
            for (j = 0; j < egcNR; j++)
            {
                if (groups->grpnr[j] != NULL)
                {
                    new_egrp[j][i-rm] = groups->grpnr[j][i];
                }
            }
            copy_rvec(state->x[i], x_tmp[i-rm]);
            copy_rvec(state->v[i], v_tmp[i-rm]);
            for (j = 0; j < ins_at->nr; j++)
            {
                if (i == ins_at->index[j])
                {
                    ins_at->index[j] = i-rm;
                }
            }

            for (j = 0; j < pos_ins->pieces; j++)
            {
                for (k = 0; k < pos_ins->nidx[j]; k++)
                {
                    if (i == pos_ins->subindex[j][k])
                    {
                        pos_ins->subindex[j][k] = i-rm;
                    }
                }
            }
        }
    }
    sfree(state->x);
    state->x = x_tmp;
    sfree(state->v);
    state->v = v_tmp;

    for (i = 0; i < egcNR; i++)
    {
        if (groups->grpnr[i] != NULL)
        {
            sfree(groups->grpnr[i]);
            groups->grpnr[i] = new_egrp[i];
        }
    }

    /* remove empty molblocks */
    RMmolblock = 0;
    for (i = 0; i < mtop->nmolblock; i++)
    {
        if (mtop->molblock[i].nmol == 0)
        {
            RMmolblock++;
        }
        else
        {
            mtop->molblock[i-RMmolblock] = mtop->molblock[i];
        }
    }
    mtop->nmolblock -= RMmolblock;
}

/* remove al bonded interactions from mtop for the molecule to be embedded */
int rm_bonded(t_block *ins_at, gmx_mtop_t *mtop)
{
    int       i, j, m;
    int       type, natom, nmol, at, atom1 = 0, rm_at = 0;
    gmx_bool *bRM, bINS;
    /*this routine lives dangerously by assuming that all molecules of a given type are in order in the structure*/
    /*this routine does not live as dangerously as it seems. There is namely a check in init_membed to make *
     * sure that g_membed exits with a warning when there are molecules of the same type not in the *
     * ins_at index group. MGWolf 050710 */


    snew(bRM, mtop->nmoltype);
    for (i = 0; i < mtop->nmoltype; i++)
    {
        bRM[i] = TRUE;
    }

    for (i = 0; i < mtop->nmolblock; i++)
    {
        /*loop over molecule blocks*/
        type         = mtop->molblock[i].type;
        natom        = mtop->molblock[i].natoms_mol;
        nmol         = mtop->molblock[i].nmol;

        for (j = 0; j < natom*nmol && bRM[type] == TRUE; j++)
        {
            /*loop over atoms in the block*/
            at   = j+atom1; /*atom index = block index + offset*/
            bINS = FALSE;

            for (m = 0; (m < ins_at->nr) && (bINS == FALSE); m++)
            {
                /*loop over atoms in insertion index group to determine if we're inserting one*/
                if (at == ins_at->index[m])
                {
                    bINS = TRUE;
                }
            }
            bRM[type] = bINS;
        }
        atom1 += natom*nmol; /*update offset*/
        if (bRM[type])
        {
            rm_at += natom*nmol; /*increment bonded removal counter by # atoms in block*/
        }
    }

    for (i = 0; i < mtop->nmoltype; i++)
    {
        if (bRM[i])
        {
            for (j = 0; j < F_LJ; j++)
            {
                mtop->moltype[i].ilist[j].nr = 0;
            }

            for (j = F_POSRES; j <= F_VSITEN; j++)
            {
                mtop->moltype[i].ilist[j].nr = 0;
            }
        }
    }
    sfree(bRM);

    return rm_at;
}

/* Write a topology where the number of molecules is correct for the system after embedding */
static void top_update(const char *topfile, rm_t *rm_p, gmx_mtop_t *mtop)
{
    int        bMolecules         = 0;
    FILE      *fpin, *fpout;
    char       buf[STRLEN], buf2[STRLEN], *temp;
    int        i, *nmol_rm, nmol, line;
    char       temporary_filename[STRLEN];

    fpin  = gmx_ffopen(topfile, "r");
    strncpy(temporary_filename, "temp.topXXXXXX", STRLEN);
    gmx_tmpnam(temporary_filename);
    fpout = gmx_ffopen(temporary_filename, "w");

    snew(nmol_rm, mtop->nmoltype);
    for (i = 0; i < rm_p->nr; i++)
    {
        nmol_rm[rm_p->block[i]]++;
    }

    line = 0;
    while (fgets(buf, STRLEN, fpin))
    {
        line++;
        if (buf[0] != ';')
        {
            strcpy(buf2, buf);
            if ((temp = strchr(buf2, '\n')) != NULL)
            {
                temp[0] = '\0';
            }
            ltrim(buf2);
            if (buf2[0] == '[')
            {
                buf2[0] = ' ';
                if ((temp = strchr(buf2, '\n')) != NULL)
                {
                    temp[0] = '\0';
                }
                rtrim(buf2);
                if (buf2[strlen(buf2)-1] == ']')
                {
                    buf2[strlen(buf2)-1] = '\0';
                    ltrim(buf2);
                    rtrim(buf2);
                    if (gmx_strcasecmp(buf2, "molecules") == 0)
                    {
                        bMolecules = 1;
                    }
                }
                fprintf(fpout, "%s", buf);
            }
            else if (bMolecules == 1)
            {
                for (i = 0; i < mtop->nmolblock; i++)
                {
                    nmol = mtop->molblock[i].nmol;
                    sprintf(buf, "%-15s %5d\n", *(mtop->moltype[mtop->molblock[i].type].name), nmol);
                    fprintf(fpout, "%s", buf);
                }
                bMolecules = 2;
            }
            else if (bMolecules == 2)
            {
                /* print nothing */
            }
            else
            {
                fprintf(fpout, "%s", buf);
            }
        }
        else
        {
            fprintf(fpout, "%s", buf);
        }
    }

    gmx_ffclose(fpout);
    /* use gmx_ffopen to generate backup of topinout */
    fpout = gmx_ffopen(topfile, "w");
    gmx_ffclose(fpout);
    rename(temporary_filename, topfile);
}

void rescale_membed(int step_rel, gmx_membed_t membed, rvec *x)
{
    /* Set new positions for the group to embed */
    if (step_rel <= membed->it_xy)
    {
        membed->fac[0] += membed->xy_step;
        membed->fac[1] += membed->xy_step;
    }
    else if (step_rel <= (membed->it_xy+membed->it_z))
    {
        membed->fac[2] += membed->z_step;
    }
    resize(membed->r_ins, x, membed->pos_ins, membed->fac);
}

/* We would like gn to be const as well, but C doesn't allow this */
/* TODO this is utility functionality (search for the index of a
   string in a collection), so should be refactored and located more
   centrally. Also, it nearly duplicates the same string in readir.c */
static int search_string(const char *s, int ng, char *gn[])
{
    int i;

    for (i = 0; (i < ng); i++)
    {
        if (gmx_strcasecmp(s, gn[i]) == 0)
        {
            return i;
        }
    }

    gmx_fatal(FARGS,
              "Group %s selected for embedding was not found in the index file.\n"
              "Group names must match either [moleculetype] names or custom index group\n"
              "names, in which case you must supply an index file to the '-n' option\n"
              "of grompp.",
              s);

    return -1;
}

gmx_membed_t init_membed(FILE *fplog, int nfile, const t_filenm fnm[], gmx_mtop_t *mtop,
                         t_inputrec *inputrec, t_state *state, t_commrec *cr, real *cpt)
{
    char                     *ins, **gnames;
    int                       i, rm_bonded_at, fr_id, fr_i = 0, tmp_id, warn = 0;
    int                       ng, j, max_lip_rm, ins_grp_id, ins_nat, mem_nat, ntype, lip_rm, tpr_version;
    real                      prot_area;
    rvec                     *r_ins = NULL;
    t_block                  *ins_at, *rest_at;
    pos_ins_t                *pos_ins;
    mem_t                    *mem_p;
    rm_t                     *rm_p;
    gmx_groups_t             *groups;
    gmx_bool                  bExcl = FALSE;
    t_atoms                   atoms;
    t_pbc                    *pbc;
    char                    **piecename = NULL;
    gmx_membed_t              membed    = NULL;

    /* input variables */
    const char *membed_input;
    real        xy_fac           = 0.5;
    real        xy_max           = 1.0;
    real        z_fac            = 1.0;
    real        z_max            = 1.0;
    int         it_xy            = 1000;
    int         it_z             = 0;
    real        probe_rad        = 0.22;
    int         low_up_rm        = 0;
    int         maxwarn          = 0;
    int         pieces           = 1;
    gmx_bool    bALLOW_ASYMMETRY = FALSE;

    /* sanity check constants */         /* Issue a warning when: */
    const int  membed_version = 58;        /* tpr version is smaller */
    const real min_probe_rad  = 0.2199999; /* A probe radius for overlap between embedded molecule *
                                            * and rest smaller than this value is probably too small */
    const real min_xy_init    = 0.0999999; /* the initial shrinking of the molecule to embed is smaller */
    const int  min_it_xy      = 1000;      /* the number of steps to embed in xy-plane is smaller */
    const int  min_it_z       = 100;       /* the number of steps to embed in z is smaller */
    const real prot_vs_box    = 7.5;       /* molecule to embed is large (more then prot_vs_box) with respect */
    const real box_vs_prot    = 50;        /* to the box size (less than box_vs_prot) */

    snew(membed, 1);
    snew(ins_at, 1);
    snew(pos_ins, 1);

    if (MASTER(cr))
    {
        /* get input data out membed file */
        membed_input = opt2fn("-membed", nfile, fnm);
        get_input(membed_input, &xy_fac, &xy_max, &z_fac, &z_max, &it_xy, &it_z, &probe_rad, &low_up_rm,
                  &maxwarn, &pieces, &bALLOW_ASYMMETRY);

        tpr_version = get_tpr_version(ftp2fn(efTPR, nfile, fnm));
        if (tpr_version < membed_version)
        {
            gmx_fatal(FARGS, "Version of *.tpr file to old (%d). "
                      "Rerun grompp with GROMACS version 4.0.3 or newer.\n", tpr_version);
        }

        if (!EI_DYNAMICS(inputrec->eI) )
        {
            gmx_input("Change integrator to a dynamics integrator in mdp file (e.g. md or sd).");
        }

        if (PAR(cr))
        {
            gmx_input("Sorry, parallel g_membed is not yet fully functional.");
        }

        if (*cpt >= 0)
        {
            fprintf(stderr, "\nSetting -cpt to -1, because embedding cannot be restarted from cpt-files.\n");
            *cpt = -1;
        }
        groups = &(mtop->groups);
        snew(gnames, groups->ngrpname);
        for (i = 0; i < groups->ngrpname; i++)
        {
            gnames[i] = *(groups->grpname[i]);
        }

        atoms = gmx_mtop_global_atoms(mtop);
        snew(mem_p, 1);
        fprintf(stderr, "\nSelect a group to embed in the membrane:\n");
        get_index(&atoms, opt2fn_null("-mn", nfile, fnm), 1, &(ins_at->nr), &(ins_at->index), &ins);
        ins_grp_id = search_string(ins, groups->ngrpname, gnames);
        fprintf(stderr, "\nSelect a group to embed %s into (e.g. the membrane):\n", ins);
        get_index(&atoms, opt2fn_null("-mn", nfile, fnm), 1, &(mem_p->mem_at.nr), &(mem_p->mem_at.index), &(mem_p->name));

        pos_ins->pieces = pieces;
        snew(pos_ins->nidx, pieces);
        snew(pos_ins->subindex, pieces);
        snew(piecename, pieces);
        if (pieces > 1)
        {
            fprintf(stderr, "\nSelect pieces to embed:\n");
            get_index(&atoms, opt2fn_null("-mn", nfile, fnm), pieces, pos_ins->nidx, pos_ins->subindex, piecename);
        }
        else
        {
            /*use whole embedded group*/
            snew(pos_ins->nidx, 1);
            snew(pos_ins->subindex, 1);
            pos_ins->nidx[0]     = ins_at->nr;
            pos_ins->subindex[0] = ins_at->index;
        }

        if (probe_rad < min_probe_rad)
        {
            warn++;
            fprintf(stderr, "\nWarning %d:\nA probe radius (-rad) smaller than 0.2 nm can result "
                    "in overlap between waters and the group to embed, which will result "
                    "in Lincs errors etc.\n\n", warn);
        }

        if (xy_fac < min_xy_init)
        {
            warn++;
            fprintf(stderr, "\nWarning %d:\nThe initial size of %s is probably too smal.\n\n", warn, ins);
        }

        if (it_xy < min_it_xy)
        {
            warn++;
            fprintf(stderr, "\nWarning %d;\nThe number of steps used to grow the xy-coordinates of %s (%d)"
                    " is probably too small.\nIncrease -nxy or.\n\n", warn, ins, it_xy);
        }

        if ( (it_z < min_it_z) && ( z_fac < 0.99999999 || z_fac > 1.0000001) )
        {
            warn++;
            fprintf(stderr, "\nWarning %d;\nThe number of steps used to grow the z-coordinate of %s (%d)"
                    " is probably too small.\nIncrease -nz or the maxwarn setting in the membed input file.\n\n", warn, ins, it_z);
        }

        if (it_xy+it_z > inputrec->nsteps)
        {
            warn++;
            fprintf(stderr, "\nWarning %d:\nThe number of growth steps (-nxy + -nz) is larger than the "
                    "number of steps in the tpr.\n"
                    "(increase maxwarn in the membed input file to override)\n\n", warn);
        }

        fr_id = -1;
        if (inputrec->opts.ngfrz == 1)
        {
            gmx_fatal(FARGS, "You did not specify \"%s\" as a freezegroup.", ins);
        }

        for (i = 0; i < inputrec->opts.ngfrz; i++)
        {
            tmp_id = mtop->groups.grps[egcFREEZE].nm_ind[i];
            if (ins_grp_id == tmp_id)
            {
                fr_id = tmp_id;
                fr_i  = i;
            }
        }

        if (fr_id == -1)
        {
            gmx_fatal(FARGS, "\"%s\" not as freezegroup defined in the mdp-file.", ins);
        }

        for (i = 0; i < DIM; i++)
        {
            if (inputrec->opts.nFreeze[fr_i][i] != 1)
            {
                gmx_fatal(FARGS, "freeze dimensions for %s are not Y Y Y\n", ins);
            }
        }

        ng = groups->grps[egcENER].nr;
        if (ng == 1)
        {
            gmx_input("No energy groups defined. This is necessary for energy exclusion in the freeze group");
        }

        for (i = 0; i < ng; i++)
        {
            for (j = 0; j < ng; j++)
            {
                if (inputrec->opts.egp_flags[ng*i+j] == EGP_EXCL)
                {
                    bExcl = TRUE;
                    if ( (groups->grps[egcENER].nm_ind[i] != ins_grp_id) ||
                         (groups->grps[egcENER].nm_ind[j] != ins_grp_id) )
                    {
                        gmx_fatal(FARGS, "Energy exclusions \"%s\" and  \"%s\" do not match the group "
                                  "to embed \"%s\"", *groups->grpname[groups->grps[egcENER].nm_ind[i]],
                                  *groups->grpname[groups->grps[egcENER].nm_ind[j]], ins);
                    }
                }
            }
        }

        if (!bExcl)
        {
            gmx_input("No energy exclusion groups defined. This is necessary for energy exclusion in "
                      "the freeze group");
        }

        /* Obtain the maximum and minimum coordinates of the group to be embedded */
        snew(rest_at, 1);
        ins_nat = init_ins_at(ins_at, rest_at, state, pos_ins, groups, ins_grp_id, xy_max);
        /* Check that moleculetypes in insertion group are not part of the rest of the system */
        check_types(ins_at, rest_at, mtop);

        mem_nat = init_mem_at(mem_p, mtop, state->x, state->box, pos_ins);

        prot_area = est_prot_area(pos_ins, state->x, ins_at, mem_p);
        if ( (prot_area > prot_vs_box) && ( (state->box[XX][XX]*state->box[YY][YY]-state->box[XX][YY]*state->box[YY][XX]) < box_vs_prot) )
        {
            warn++;
            fprintf(stderr, "\nWarning %d:\nThe xy-area is very small compared to the area of the protein.\n"
                    "This might cause pressure problems during the growth phase. Just try with\n"
                    "current setup and increase 'maxwarn' in your membed settings file, but lower the\n"
                    "compressibility in the mdp-file or disable pressure coupling if problems occur.\n\n", warn);
        }

        if (warn > maxwarn)
        {
            gmx_fatal(FARGS, "Too many warnings (override by setting maxwarn in the membed input file)\n");
        }

        printf("The estimated area of the protein in the membrane is %.3f nm^2\n", prot_area);
        printf("\nThere are %d lipids in the membrane part that overlaps the protein.\n"
               "The area per lipid is %.4f nm^2.\n", mem_p->nmol, mem_p->lip_area);

        /* Maximum number of lipids to be removed*/
        max_lip_rm = (int)(2*prot_area/mem_p->lip_area);
        printf("Maximum number of lipids that will be removed is %d.\n", max_lip_rm);

        printf("\nWill resize the protein by a factor of %.3f in the xy plane and %.3f in the z direction.\n"
               "This resizing will be done with respect to the geometrical center of all protein atoms\n"
               "that span the membrane region, i.e. z between %.3f and %.3f\n\n",
               xy_fac, z_fac, mem_p->zmin, mem_p->zmax);

        /* resize the protein by xy and by z if necessary*/
        snew(r_ins, ins_at->nr);
        init_resize(ins_at, r_ins, pos_ins, mem_p, state->x, bALLOW_ASYMMETRY);
        membed->fac[0] = membed->fac[1] = xy_fac;
        membed->fac[2] = z_fac;

        membed->xy_step = (xy_max-xy_fac)/(double)(it_xy);
        membed->z_step  = (z_max-z_fac)/(double)(it_z-1);

        resize(r_ins, state->x, pos_ins, membed->fac);

        /* remove overlapping lipids and water from the membrane box*/
        /*mark molecules to be removed*/
        snew(pbc, 1);
        set_pbc(pbc, inputrec->ePBC, state->box);

        snew(rm_p, 1);
        lip_rm = gen_rm_list(rm_p, ins_at, rest_at, pbc, mtop, state->x, mem_p, pos_ins,
                             probe_rad, low_up_rm, bALLOW_ASYMMETRY);
        lip_rm -= low_up_rm;

        if (fplog)
        {
            for (i = 0; i < rm_p->nr; i++)
            {
                fprintf(fplog, "rm mol %d\n", rm_p->mol[i]);
            }
        }

        for (i = 0; i < mtop->nmolblock; i++)
        {
            ntype = 0;
            for (j = 0; j < rm_p->nr; j++)
            {
                if (rm_p->block[j] == i)
                {
                    ntype++;
                }
            }
            printf("Will remove %d %s molecules\n", ntype, *(mtop->moltype[mtop->molblock[i].type].name));
        }

        if (lip_rm > max_lip_rm)
        {
            warn++;
            fprintf(stderr, "\nWarning %d:\nTrying to remove a larger lipid area than the estimated "
                    "protein area\nTry making the -xyinit resize factor smaller or increase "
                    "maxwarn in the membed input file.\n\n", warn);
        }

        /*remove all lipids and waters overlapping and update all important structures*/
        rm_group(groups, mtop, rm_p, state, ins_at, pos_ins);

        rm_bonded_at = rm_bonded(ins_at, mtop);
        if (rm_bonded_at != ins_at->nr)
        {
            fprintf(stderr, "Warning: The number of atoms for which the bonded interactions are removed is %d, "
                    "while %d atoms are embedded. Make sure that the atoms to be embedded are not in the same"
                    "molecule type as atoms that are not to be embedded.\n", rm_bonded_at, ins_at->nr);
        }

        if (warn > maxwarn)
        {
            gmx_fatal(FARGS, "Too many warnings.\nIf you are sure these warnings are harmless,\n"
                      "you can increase the maxwarn setting in the membed input file.");
        }

        if (ftp2bSet(efTOP, nfile, fnm))
        {
            top_update(opt2fn("-mp", nfile, fnm), rm_p, mtop);
        }

        sfree(pbc);
        sfree(rest_at);
        if (pieces > 1)
        {
            sfree(piecename);
        }

        membed->it_xy   = it_xy;
        membed->it_z    = it_z;
        membed->pos_ins = pos_ins;
        membed->r_ins   = r_ins;
    }

    return membed;
}
