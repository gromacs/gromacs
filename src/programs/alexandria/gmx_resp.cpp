/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/legacyheaders/nmsimplex.h"
#include "coulombintegrals/coulombintegrals.h"
#include "poldata.h"
#include "gmx_resp.h"
#include "gentop_qgen.h"
#include "stringutil.h"

bool gmx_ra_init(gmx_ra *ra, int atomnumber, int atype,
                 const char *atomtype, gmx_poldata_t pd,
                 ChargeGenerationModel iModel, char **dzatoms)
{
    int  k, zz;
    bool bRestr;

    bRestr = false;
    if (NULL != dzatoms)
    {
        k = 0;
        while ((NULL != dzatoms[k]) && !bRestr)
        {
            bRestr = (strcasecmp(atomtype, dzatoms[k]) == 0);
            k++;
        }
    }
    ra->atomnumber  = atomnumber;
    ra->atype       = atype;
    ra->atomtype    = strdup(atomtype);
    ra->nZeta       = gmx_poldata_get_nzeta(pd, iModel, ra->atomtype);
    if (ra->nZeta <= 0)
    {
        return false;
    }

    ra->bRestrained = bRestr;

    snew(ra->zeta, ra->nZeta);
    snew(ra->zeta_ref, ra->nZeta);
    snew(ra->q, ra->nZeta);
    snew(ra->iz, ra->nZeta);
    snew(ra->iq, ra->nZeta);
    snew(ra->row, ra->nZeta);

    for (zz = 0; (zz < ra->nZeta); zz++)
    {
        ra->iq[zz]       = -1;
        ra->q[zz]        = gmx_poldata_get_q(pd, iModel, ra->atomtype, zz);
        ra->iz[zz]       = -1;
        ra->zeta_ref[zz] =
            ra->zeta[zz] = gmx_poldata_get_zeta(pd, iModel, ra->atomtype, zz);
        ra->row[zz]      = gmx_poldata_get_row(pd, iModel, ra->atomtype, zz);
    }
    return true;
}

void gmx_ra_done(gmx_ra *ra)
{
    sfree(ra->iz);
    sfree(ra->iq);
    sfree(ra->zeta);
    sfree(ra->q);
    sfree(ra->row);
    sfree(ra->atomtype);
}

real gmx_ra_get_q(gmx_ra *ra)
{
    int  i;
    real q = 0;

    for (i = 0; (i < ra->nZeta); i++)
    {
        q += ra->q[i];
    }

    return q;
}

gmx_resp_t gmx_resp_init(ChargeGenerationModel iModel,
                         bool bAXpRESP, real qfac, real b_hyper, real qtot,
                         real zmin, real zmax, real delta_z, bool bZatype,
                         real watoms, real rDecrZeta, bool bRandZeta,
                         bool bRandQ, real penalty_fac, bool bFitZeta, bool bEntropy,
                         const char *dzatoms, unsigned int seed)
{
    gmx_resp_t gr;

    snew(gr, 1);
    gr->qtot      = qtot;
    gr->qsum      = qtot;
    gr->bAXpRESP  = bAXpRESP;
    gr->qfac      = qfac;
    gr->b_hyper   = b_hyper;
    gr->wtot      = 0;
    gr->iModel    = iModel;
    gr->zmin      = zmin;
    gr->zmax      = zmax;
    gr->delta_z   = delta_z;
    std::vector<std::string> ptr = split(dzatoms, ' ');
    snew(gr->dzatoms, ptr.size());
    for (unsigned int i = 0; (i < ptr.size()); ++i)
    {
        gr->dzatoms[i] = strdup(ptr[i].c_str());
    }
    gr->pfac      = penalty_fac;
    gr->qmin      = -2;
    gr->qmax      = 2; /* e */
    gr->nesp      = 0;
    gr->nrho      = 0;
    gr->natom     = 0;
    gr->natype    = 0;
    gr->seed      = seed;
    gr->bEntropy  = bEntropy;
    gr->bZatype   = bZatype;
    gr->rDecrZeta = rDecrZeta;
    gr->bRandZeta = bRandZeta;
    gr->bRandQ    = bRandQ;
    gr->bFitZeta  = bFitZeta;
    gr->watoms    = watoms;
    gr->nparam    = 0;

    return gr;
}

void gmx_resp_get_atom_info(gmx_resp_t gr, t_atoms *atoms,
                            t_symtab *symtab, rvec **x)
{
    int          i;
    const char  *rnm;

    init_t_atoms(atoms, gr->natom, true);
    if (NULL == (*x))
    {
        snew((*x), atoms->nr);
    }
    if (NULL == atoms->atomtype)
    {
        snew(atoms->atomtype, atoms->nr);
    }
    for (i = 0; (i < gr->natom); i++)
    {
        atoms->atom[i].atomnumber = gr->ra[i].atomnumber;
        atoms->atom[i].q          = gmx_ra_get_q(&gr->ra[i]);
        atoms->atomname[i]        = put_symtab(symtab, gr->ra[i].atomtype);
        atoms->atomtype[i]        = put_symtab(symtab, gr->ra[i].atomtype);
        atoms->atom[i].resind     = 0;

        strncpy(atoms->atom[i].elem, gr->ra[i].atomtype,
                sizeof(atoms->atom[i].elem)-1);
        copy_rvec(gr->x[i], (*x)[i]);
    }
    rnm = (NULL != gr->stoichiometry) ? gr->stoichiometry : (const char *)"BOE";
    t_atoms_set_resinfo(atoms, 0, symtab, rnm, 1, ' ', 1, ' ');

    atoms->nres = 1;
}

void gmx_resp_update_atomtypes(gmx_resp_t gr, t_atoms *atoms)
{
    int i, j;

    for (i = 0; (i < gr->natom); i++)
    {
        gr->ra[i].atomtype = strdup(*atoms->atomtype[i]);
        for (j = 0; (j < i); j++)
        {
            if (0 == strcmp(*atoms->atomtype[i], *atoms->atomtype[j]))
            {
                break;
            }
        }
        if (j == i)
        {
            gr->natype++;
        }
        gr->ra[i].atype = j;
    }
}

void gmx_resp_add_atom_coords(gmx_resp_t gr, rvec *x)
{
    int        i;

    srenew(gr->x, gr->natom);
    for (i = 0; (i < gr->natom); i++)
    {
        copy_rvec(x[i], gr->x[i]);
    }
}

void gmx_resp_fill_zeta(gmx_resp_t gr, gmx_poldata_t pd)
{
    int i, zz;

    for (i = 0; (i < gr->natom); i++)
    {
        for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
        {
            gmx_resp_set_zeta(gr, i, zz, gmx_poldata_get_zeta(pd, gr->iModel,
                                                              gr->ra[i].atomtype, zz));
        }
    }
}

void gmx_resp_fill_q(gmx_resp_t gr, t_atoms *atoms)
{
    int    i, zz;
    double q;

    for (i = 0; (i < gr->natom); i++)
    {
        q = 0;
        for (zz = 0; (zz < gr->ra[i].nZeta-1); zz++)
        {
            q -= gr->ra[i].q[zz];
        }
        q += atoms->atom[i].q;
        gmx_resp_set_q(gr, i, gr->ra[i].nZeta-1, q);
    }
}

bool gmx_resp_add_atom_info(gmx_resp_t gr, t_atoms *atoms, gmx_poldata_t pd)
{
    int  i;

    gr->natom    = atoms->nr;
    snew(gr->ra, gr->natom);

    for (i = 0; (i < gr->natom); i++)
    {
        if (!gmx_ra_init(&gr->ra[i],
                         atoms->atom[i].atomnumber, atoms->atom[i].type,
                         *(atoms->atomtype[i]), pd, gr->iModel, gr->dzatoms))
        {
            return false;
        }
    }
    return true;
}

void gmx_resp_summary(FILE *fp, gmx_resp_t gr, 
                      std::vector<int> &symmetric_atoms)
{
    int i;

    if (NULL != fp)
    {
        fprintf(fp, "There are %d atoms, %d atomtypes %d parameters for (R)ESP fitting.\n",
                gr->natom, gr->natype, gr->nparam);
        for (i = 0; (i < gr->natom); i++)
        {
            fprintf(fp, " %d", symmetric_atoms[i]);
        }
        fprintf(fp, "\n");
    }
}

enum eParm {
    eparmQ, eparmZ, eparmNR
};

void gmx_resp_add_param(gmx_resp_t gr, int atom, eParm eparm, int zz)
{
    range_check(atom, 0, gr->natom);
    if ((zz >= 0) && (zz < gr->ra[atom].nZeta))
    {
        if (eparm == eparmQ)
        {
            gr->ra[atom].iq[zz] = gr->nparam++;
            if (debug)
            {
                fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
            }
        }
        else if (gr->bFitZeta)
        {
            if (gr->ra[atom].zeta[zz] != 0)
            {
                gr->ra[atom].iz[zz] = gr->nparam++;
                if (debug)
                {
                    fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
                }
            }
            else
            {
                gr->ra[atom].iz[zz] = -1;
            }
        }
    }
}

void gmx_resp_add_atom_symmetry(gmx_resp_t gr, 
                                std::vector<int> &symmetric_atoms)
{
    int        i, k, zz;

    if (NULL == gr->ra)
    {
        gmx_fatal(FARGS, "resp_atom struct not initialized");
    }

    /* Map the symmetric atoms */
    for (i = 0; (i < gr->natom); i++)
    {
        if (0 == i)
        {
            /* The first charge is not a free variable, it follows from the total charge.
             * Only add the zeta values here.
             */
            for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
            {
                gmx_resp_add_param(gr, i, eparmZ, zz);
            }
        }
        else if (symmetric_atoms[i] == i)
        {
            gmx_resp_add_param(gr, i, eparmQ, gr->ra[i].nZeta-1);

            if (gr->bZatype)
            {
                for (k = 0; (k < i); k++)
                {
                    if (gr->ra[i].atype == gr->ra[k].atype)
                    {
                        break;
                    }
                }
                if (k == i)
                {
                    for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                    {
                        gmx_resp_add_param(gr, i, eparmZ, zz);
                    }
                }
                else
                {
                    for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                    {
                        gr->ra[i].iz[zz] = gr->ra[k].iz[zz];
                    }
                }
            }
            else
            {
                for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                {
                    gmx_resp_add_param(gr, i, eparmZ, zz);
                }
            }
        }
        else if (symmetric_atoms[i] > i)
        {
            gmx_fatal(FARGS, "The symmetric_atoms array can not point to larger atom numbers");
        }
        else if (gr->ra[i].nZeta > 0)
        {
            gr->ra[i].iq[gr->ra[i].nZeta-1] =
                gr->ra[symmetric_atoms[i]].iq[gr->ra[i].nZeta-1];
            for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
            {
                gr->ra[i].iz[zz] = gr->ra[symmetric_atoms[i]].iz[zz];
            }
            if (debug)
            {
                fprintf(debug, "Atom %d is a copy of atom %d\n",
                        i+1, symmetric_atoms[i]+1);
            }
        }
        else
        {
            if (0 == i)
            {
                gr->ra[i].iq[gr->ra[i].nZeta-1] = -1;
            }
            else
            {
                gmx_resp_add_param(gr, i, eparmQ, gr->ra[i].nZeta-1);
            }
            
            for (k = 0; (k < gr->ra[i].nZeta); k++)
            {
                gmx_resp_add_param(gr, i, eparmZ, k);
            }
        }
    }
    if (debug)
    {
        int maxz = 0;
        for (i = 0; (i < gr->natom); i++)
        {
            if (gr->ra[i].nZeta > maxz)
            {
                maxz = gr->ra[i].nZeta;
            }
        }
        
        fprintf(debug, "GRQ: %3s %5s", "nr", "type");
        for (i = 0; (i < maxz); i++)
        {
            fprintf(debug, " %8s %4s %8s %4s\n", "q", "iq", "zeta", "iz");
        }
        fprintf(debug, "\n");
        for (i = 0; (i < gr->natom); i++)
        {
            fprintf(debug, "GRQ: %3d %5s", i+1, gr->ra[i].atomtype);
            for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
            {
                fprintf(debug, " %8.4f %4d %8.4f %4d\n",
                        gr->ra[i].q[zz], gr->ra[i].iq[zz],
                        gr->ra[i].zeta[zz], gr->ra[i].iz[zz]);
            }
            fprintf(debug, "\n");
        }
        fprintf(debug, "gr->qsum = %g\n", gr->qsum);
    }
}

void gmx_resp_write_histo(gmx_resp_t gr, const char *fn, char *title, output_env_t oenv)
{
    FILE       *fp;
    gmx_stats_t gs;
    real       *x, *y;
    int         i, nbin = 100;

    if (NULL == fn)
    {
        return;
    }
    gs = gmx_stats_init();
    for (i = 0; (i < gr->nesp); i++)
    {
        gmx_stats_add_point(gs, i, gmx2convert(gr->pot_calc[i], eg2cHartree_e), 0, 0);
    }

    gmx_stats_make_histogram(gs, 0, &nbin, ehistoY, 1, &x, &y);

    fp = xvgropen(fn, title, "Pot (1/a.u.)", "()", oenv);
    for (i = 0; (i < nbin); i++)
    {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
    }
    sfree(x);
    sfree(y);
    fclose(fp);
    gmx_stats_done(gs);
}

void gmx_resp_write_diff_cube(gmx_resp_t grref, gmx_resp_t gr, const char *cube_fn,
                              const char *hist_fn, char *title, output_env_t oenv,
                              int rho)
{
    FILE       *fp;
    int         i, m, ix, iy, iz, zz;
    real        pp, q, r, rmin;
    rvec        dx;
    gmx_stats_t gst = NULL, ppcorr = NULL;

    if (NULL != hist_fn)
    {
        gst    = gmx_stats_init();
        ppcorr = gmx_stats_init();
    }
    if (NULL != cube_fn)
    {
        fp = gmx_ffopen(cube_fn, "w");
        fprintf(fp, "%s\n", title);
        fprintf(fp, "POTENTIAL\n");
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n",
                gr->natom,
                gmx2convert(gr->origin[XX], eg2cBohr),
                gmx2convert(gr->origin[YY], eg2cBohr),
                gmx2convert(gr->origin[ZZ], eg2cBohr));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", gr->nxyz[XX],
                gmx2convert(gr->space[XX], eg2cBohr), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", gr->nxyz[YY],
                0.0, gmx2convert(gr->space[YY], eg2cBohr), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", gr->nxyz[ZZ],
                0.0, 0.0, gmx2convert(gr->space[ZZ], eg2cBohr));

        for (m = 0; (m < gr->natom); m++)
        {
            q = 0;
            for (zz = 0; (zz < gr->ra[m].nZeta); zz++)
            {
                q += gr->ra[m].q[zz];
            }
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    gr->ra[m].atomnumber, q,
                    gmx2convert(gr->x[m][XX], eg2cBohr),
                    gmx2convert(gr->x[m][YY], eg2cBohr),
                    gmx2convert(gr->x[m][ZZ], eg2cBohr));
        }

        for (ix = m = 0; ix < gr->nxyz[XX]; ix++)
        {
            for (iy = 0; iy < gr->nxyz[YY]; iy++)
            {
                for (iz = 0; iz < gr->nxyz[ZZ]; iz++, m++)
                {
                    if (NULL != grref)
                    {
                        pp = gr->pot_calc[m] - grref->pot[m];
                        if (NULL != ppcorr)
                        {
                            gmx_stats_add_point(ppcorr,
                                                gmx2convert(grref->pot[m], eg2cHartree_e),
                                                gmx2convert(gr->pot_calc[m], eg2cHartree_e), 0, 0);
                        }
                    }
                    else
                    {
                        if (rho == 0)
                        {
                            pp = gmx2convert(gr->pot_calc[m], eg2cHartree_e);
                        }
                        else
                        {
                            pp = gr->rho[m]*pow(BOHR2NM, 3);
                        }
                    }
                    fprintf(fp, "%13.5e", pp);
                    if (iz % 6 == 5)
                    {
                        fprintf(fp, "\n");
                    }
                    if (NULL != gst)
                    {
                        rmin = 1000;
                        /* Add point to histogram! */
                        for (i = 0; (i < gr->natom); i++)
                        {
                            rvec_sub(gr->x[i], gr->esp[m], dx);
                            r = norm(dx);
                            if (r < rmin)
                            {
                                rmin = r;
                            }
                        }
                        gmx_stats_add_point(gst, rmin, pp, 0, 0);
                    }
                }
                if ((iz % 6) != 0)
                {
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
    }
    if (NULL != gst)
    {
        int   nb = 0;
        real *x  = NULL, *y = NULL;

        fp = xvgropen(hist_fn, "Absolute deviation from QM", "Distance (nm)",
                      "Potential", oenv);
        gmx_stats_dump_xy(gst, fp);
        if (0)
        {
            gmx_stats_make_histogram(gst, 0.01, &nb, ehistoX, 0, &x, &y);
            gmx_stats_done(gst);
            for (i = 0; (i < nb); i++)
            {
                fprintf(fp, "%10g  %10g\n", x[i], y[i]);
            }
            sfree(x);
            sfree(y);
        }
        fclose(fp);
        fp = xvgropen("diff-pot.xvg", "Correlation between QM and Calc", "Pot (QM)",
                      "Pot (Calc)", oenv);
        gmx_stats_dump_xy(ppcorr, fp);
        fclose(fp);
    }
}

void gmx_resp_write_cube(gmx_resp_t gr, const char *fn, char *title)
{
    gmx_resp_write_diff_cube(NULL, gr, fn, NULL, title, NULL, 0);
}

void gmx_resp_write_rho(gmx_resp_t gr, const char *fn, char *title)
{
    gmx_resp_write_diff_cube(NULL, gr, fn, NULL, title, NULL, 1);
}

void gmx_resp_read_cube(gmx_resp_t gr, const char *fn, bool bESPonly)
{
    char         **strings;
    bool           bOK;
    double         lx, ly, lz, pp, qq;
    int            nlines, line = 0, m, ix, iy, iz, n, anr, nxyz[DIM];
    double         origin[DIM], space[DIM];
    const  char   *forms[] = {
        "%lf", "%*s%lf", "%*s%*s%lf", "%*s%*s%*s%lf",
        "%*s%*s%*s%*s%lf", "%*s%*s%*s%*s%*s%lf"
    };
    if (NULL == fn)
    {
        return;
    }

    nlines = get_file(fn, &strings);
    bOK    = (nlines > 100);
    if (bOK)
    {
        printf("%s\n", strings[line++]);
    }
    bOK = (line < nlines) && (strcmp(strings[line++], "POTENTIAL") != 0);
    if (bOK)
    {
        bOK = (line < nlines) && (4 == sscanf(strings[line++], "%d%lf%lf%lf",
                                              &n, &origin[XX], &origin[YY], &origin[ZZ]));
    }
    if (bOK && !bESPonly)
    {
        gr->natom      = n;
        gr->origin[XX] = origin[XX];
        gr->origin[YY] = origin[YY];
        gr->origin[ZZ] = origin[ZZ];
    }
    if (bOK)
    {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%lf",
                                              &nxyz[XX], &space[XX]));
    }
    if (bOK)
    {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%*s%lf",
                                              &nxyz[YY], &space[YY]));
    }
    if (bOK)
    {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%*s%*s%lf",
                                              &nxyz[ZZ], &space[ZZ]));
    }
    if (bOK)
    {
        for (m = 0; (m < DIM); m++)
        {
            gr->nxyz[m]  = nxyz[m];
            gr->space[m] = space[m];
        }
        for (m = 0; (m < DIM); m++)
        {
            gr->origin[m] = convert2gmx(gr->origin[m], eg2cBohr);
            gr->space[m]  = convert2gmx(gr->space[m], eg2cBohr);
        }
    }
    if (bOK && ((line+gr->natom) < nlines))
    {
        snew(gr->x, gr->natom);
        for (m = 0; (m < gr->natom); m++)
        {
            bOK = (5 == sscanf(strings[line++], "%d%lf%lf%lf%lf",
                               &anr, &qq, &lx, &ly, &lz));
            if (bOK)
            {
                if (!bESPonly)
                {
                    gr->ra[m].atomnumber = anr;
                    if (gr->ra[m].nZeta > 0)
                    {
                        gr->ra[m].q[gr->ra[m].nZeta-1] = qq;
                    }
                }
                gr->x[m][XX] = convert2gmx(lx, eg2cBohr);
                gr->x[m][YY] = convert2gmx(ly, eg2cBohr);
                gr->x[m][ZZ] = convert2gmx(lz, eg2cBohr);
            }
        }
    }
    if (bOK)
    {
        gr->nesp = gr->nxyz[XX]*gr->nxyz[YY]*gr->nxyz[ZZ];
        snew(gr->pot, gr->nesp);
        snew(gr->esp, gr->nesp);
        for (ix = m = 0; ix < gr->nxyz[XX]; ix++)
        {
            for (iy = 0; iy < gr->nxyz[YY]; iy++)
            {
                for (iz = 0; iz < gr->nxyz[ZZ]; iz++, m++)
                {
                    gr->esp[m][XX] = gr->origin[XX] + ix*gr->space[XX];
                    gr->esp[m][YY] = gr->origin[YY] + iy*gr->space[YY];
                    gr->esp[m][ZZ] = gr->origin[ZZ] + iz*gr->space[ZZ];
                    bOK            = (1 == sscanf(strings[line], forms[iz % 6], &pp));
                    if (bOK)
                    {
                        gr->pot[m] = convert2gmx(pp, eg2cHartree_e);
                    }
                    if (iz % 6 == 5)
                    {
                        line++;
                    }
                }
                if ((iz % 6) != 0)
                {
                    line++;
                }
            }
        }
    }
    bOK = (line == nlines);
    if (!bOK)
    {
        gmx_fatal(FARGS, "Error reading %s, line %d out of %d", fn, line, nlines);
    }

    for (m = 0; (m < nlines); m++)
    {
        sfree(strings[m]);
    }
    sfree(strings);
}

void gmx_resp_copy_grid(gmx_resp_t dest, gmx_resp_t src)
{
    int m;

    for (m = 0; (m < DIM); m++)
    {
        dest->origin[m] = src->origin[m];
        dest->space[m]  = src->space[m];
        dest->nxyz[m]   = src->nxyz[m];
    }
    dest->nesp = src->nesp;
    snew(dest->esp, dest->nesp);
    snew(dest->pot, dest->nesp);
    snew(dest->pot_calc, dest->nesp);
    for (m = 0; (m < dest->nesp); m++)
    {
        copy_rvec(src->esp[m], dest->esp[m]);
    }
}

gmx_resp_t gmx_resp_copy(gmx_resp_t src)
{
    gmx_resp_t dest;

    snew(dest, 1);
    memcpy(dest, src, sizeof(*src));

    return dest;
}

void gmx_resp_make_grid(gmx_resp_t gr, real spacing, matrix box, rvec x[])
{
    int  i, j, k, m, n;
    rvec xyz;

    if (0 != gr->nesp)
    {
        fprintf(stderr, "Overwriting existing ESP grid\n");
    }
    if (0 <= spacing)
    {
        spacing = 0.1;
        fprintf(stderr, "spacing too small, setting it to %g\n", spacing);
    }
    snew(gr->x, gr->natom);
    for (i = 0; (i < gr->natom); i++)
    {
        copy_rvec(x[i], gr->x[i]);
    }
    gr->nesp = 1;
    for (m = 0; (m < DIM); m++)
    {
        gr->nxyz[m]  = 1+(int) (box[m][m]/spacing);
        gr->space[m] = box[m][m]/gr->nxyz[m];
        gr->nesp    *= gr->nxyz[m];
    }
    n = 0;
    snew(gr->esp, gr->nesp);
    snew(gr->pot_calc, gr->nesp);
    for (i = 0; (i < gr->nxyz[XX]); i++)
    {
        xyz[XX] = (i-0.5*gr->nxyz[XX])*gr->space[XX];
        for (j = 0; (j < gr->nxyz[YY]); j++)
        {
            xyz[YY] = (j-0.5*gr->nxyz[YY])*gr->space[YY];
            for (k = 0; (k < gr->nxyz[ZZ]); k++)
            {
                xyz[ZZ] = (k-0.5*gr->nxyz[ZZ])*gr->space[ZZ];
                copy_rvec(xyz, gr->esp[n]);
                n++;
            }
        }
    }
}

void gmx_resp_calc_rho(gmx_resp_t gr)
{
    int  i, j, k;
    real r, z, V, vv, pi32;
    rvec dx;

    pi32 = pow(M_PI, -1.5);
    if (gr->nrho < gr->nesp)
    {
        srenew(gr->rho, gr->nesp);
        gr->nrho = gr->nesp;
    }
    for (i = 0; (i < gr->nrho); i++)
    {
        V = 0;
        for (j = 0; (j < gr->natom); j++)
        {
            vv = 0;
            rvec_sub(gr->esp[i], gr->x[j], dx);
            r = norm(dx);
            switch (gr->iModel)
            {
                case eqgBultinck:
                case eqgAXp:
                    return;
                case eqgAXs:
                    vv = 0;
                    break;
                case eqgYang:
                case eqgRappe:
                    vv = gr->ra[j].q[0]*Nuclear_SS(r, gr->ra[j].row[0],
                                                   gr->ra[j].zeta[0]);
                    break;
                case eqgAXg:
                    vv = 0;
                    for (k = 0; (k < gr->ra[j].nZeta); k++)
                    {
                        z = gr->ra[j].zeta[k];
                        if (z > 0)
                        {
                            vv -= (gr->ra[j].q[k]*pi32*exp(-sqr(r*z))*
                                   pow(z, 3));
                        }
                    }
                    break;
                default:
                    gmx_fatal(FARGS, "Krijg nou wat, iModel = %d!",
                              gr->iModel);
            }
            V  += vv;
        }
        gr->rho[i] = V;
    }
}

void gmx_resp_calc_pot(gmx_resp_t gr)
{
    int    i, j, k, m;
    double r, r2, dx, V, vv;

    for (i = 0; (i < gr->nesp); i++)
    {
        V = 0;
        for (j = 0; (j < gr->natom); j++)
        {
            vv = 0;
            r2 = 0;
            for (m = 0; (m < DIM); m++)
            {
                dx  = gr->esp[i][m]-gr->x[j][m];
                r2 += dx*dx;
            }
            r = sqrt(r2);
            switch (gr->iModel)
            {
                case eqgBultinck:
                case eqgAXp:
                    if (r > 0.01)
                    {
                        vv = gr->ra[j].q[0]/r;
                    }
                    break;
                case eqgAXs:
                    vv = 0;
                    for (k = 0; (k < gr->ra[j].nZeta); k++)
                    {
                        vv += gr->ra[j].q[k]*Nuclear_SS(r, gr->ra[j].row[k],
                                                        gr->ra[j].zeta[k]);
                    }
                    break;
                case eqgYang:
                case eqgRappe:
                    vv = gr->ra[j].q[0]*Nuclear_SS(r, gr->ra[j].row[0],
                                                   gr->ra[j].zeta[0]);
                    break;
                case eqgAXg:
                    vv = 0;
                    for (k = 0; (k < gr->ra[j].nZeta); k++)
                    {
                        vv += gr->ra[j].q[k]*Nuclear_GG(r, gr->ra[j].zeta[k]);
                    }
                    break;
                default:
                    gmx_fatal(FARGS, "Krijg nou wat, iModel = %s!",
                              get_eemtype_name(gr->iModel));
            }
            V  += vv;
        }
        gr->pot_calc[i] = V*ONE_4PI_EPS0;
    }
}

static void gmx_resp_warning(const char *fn, int line)
{
    fprintf(stderr, "WARNING: It seems like you have two sets of ESP data in your file\n         %s\n", fn);
    fprintf(stderr, "         using the second set, starting at line %d\n", line);
}

const char *gmx_resp_get_stoichiometry(gmx_resp_t gr)
{
    return gr->stoichiometry;
}

static void get_set_vector(gmx_resp_t   gr,
                           bool         bSet,
                           bool         bRandQ,
                           bool         bRandZeta,
                           unsigned int seed,
                           double      *nmx)
{
    int       i, n, zz, zzz, nrest;
    double    qtot, dq, qi, zeta;
    gmx_rng_t rnd = NULL;

    if (bSet && (bRandQ || bRandZeta))
    {
        rnd = gmx_rng_init(seed);
    }
    gr->penalty = 0;
    n           = 0;
    qtot        = 0;
    nrest       = 0;
    for (i = 0; (i < gr->natom); i++)
    {
        if (bSet)
        {
            /* First do charges */
            qi = 0;
            for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
            {
                if (gr->ra[i].iq[zz] == n)
                {
                    if (gr->ra[i].q[zz] == 0)
                    {
                        nmx[n] = -qi;
                        if (bRandQ)
                        {
                            nmx[n] += 0.2*(gmx_rng_uniform_real(rnd)-0.5);
                        }
                    }
                    else
                    {
                        nmx[n] = gr->ra[i].q[zz];
                    }
                    n++;
                }
                qi += gr->ra[i].q[zz];
            }
            /* Then do zeta */
            if (gr->bFitZeta)
            {
                for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                {
                    if (gr->ra[i].iz[zz] == n)
                    {
                        real zmin = gr->zmin;
                        real zmax = gr->zmax;

                        if ((gr->delta_z > 0) && (gr->ra[i].bRestrained))
                        {
                            zmax = gr->ra[i].zeta_ref[zz]+gr->delta_z;
                            zmin = gr->ra[i].zeta_ref[zz]-gr->delta_z;
                        }
                        if ((zz > 1) && (gr->rDecrZeta >= 0))
                        {
                            zmax = gr->ra[i].zeta[zz-1]-gr->rDecrZeta;
                            if (zmax < zmin)
                            {
                                zmax = (zmin+gr->ra[i].zeta[zz-1])/2;
                            }
                        }
                        if (bRandZeta)
                        {
                            nmx[n] = zmin + (zmax-zmin)*gmx_rng_uniform_real(rnd);
                        }
                        else
                        {
                            nmx[n] = gr->ra[i].zeta[zz];
                        }
                        gr->ra[i].zeta[zz] = nmx[n];
                        n++;
                    }
                }
            }
        }
        else
        {
            /* Initialize to something strange */
            if (gr->bFitZeta)
            {
                for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                {
                    if (gr->ra[i].zeta[zz] != 0)
                    {
                        gr->ra[i].zeta[zz] = NOTSET;
                    }
                }
            }
            gr->ra[i].q[gr->ra[i].nZeta-1] = NOTSET;

            /* First do charges */
            for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
            {
                if (gr->ra[i].iq[zz] == n)
                {
                    gr->ra[i].q[zz] = nmx[n];
                    qtot           += gr->ra[i].q[zz];
                    n++;
                }
                else if ((gr->ra[i].iq[zz] < n) && (gr->ra[i].iq[zz] >= 0))
                {
                    for (zzz = 0; (zzz < i); zzz++)
                    {
                        if (gr->ra[zzz].iq[zz] == gr->ra[i].iq[zz])
                        {
                            gr->ra[i].q[zz] = gr->ra[zzz].q[zz];
                            break;
                        }
                    }
                    if (zzz == i)
                    {
                        gmx_fatal(FARGS, "Can not find a previous atom with iq[%d] = %d", zz, n);
                    }

                    /* Only sum those atoms to qtot, that are not part of
                       the "rest" charge */
                    if (gr->ra[i].iq[zz] != -1)
                    {
                        qtot += gr->ra[i].q[zz];
                    }
                }
                else if (zz == gr->ra[i].nZeta-1)
                {
                    nrest++;
                }
                else
                {
                    qtot += gr->ra[i].q[zz];
                }
            }

            if (gr->bFitZeta)
            {
                /* Then do zeta */
                for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
                {
                    if (gr->ra[i].iz[zz] == n)
                    {
                        zeta               = nmx[n];
                        gr->ra[i].zeta[zz] = zeta;
                        if (gr->delta_z >= 0)
                        {
                            real zmin = gr->ra[i].zeta_ref[zz]-gr->delta_z;
                            real zmax = gr->ra[i].zeta_ref[zz]+gr->delta_z;
                            if (zeta <= zmin)
                            {
                                gr->penalty += sqr(zeta-zmin);
                            }
                            else if (zeta >= zmax)
                            {
                                gr->penalty += sqr(zmax-zeta);
                            }
                        }
                        else
                        {
                            if (zeta <= gr->zmin)
                            {
                                gr->penalty += sqr(gr->zmin-zeta);
                            }
                            else if (zeta >= gr->zmax)
                            {
                                gr->penalty += sqr(gr->zmax-zeta);
                            }
                            if ((gr->rDecrZeta >= 0) && (zz > 0) &&
                                (gr->ra[i].zeta[zz-1] != 0) &&
                                ((gr->ra[i].zeta[zz-1] - zeta) < gr->rDecrZeta))
                            {
                                gr->penalty += sqr(gr->ra[i].zeta[zz-1] - zeta - gr->rDecrZeta);
                            }
                        }
                        n++;
                    }
                    else if ((gr->ra[i].iz[zz] < n) && (gr->ra[i].iz[zz] >= 0))
                    {
                        for (zzz = 0; (zzz < i); zzz++)
                        {
                            if (gr->ra[zzz].iz[zz] == gr->ra[i].iz[zz])
                            {
                                gr->ra[i].zeta[zz] = gr->ra[zzz].zeta[zz];
                                break;
                            }
                        }
                        if (zzz == i)
                        {
                            gmx_fatal(FARGS, "Can not find a previous atom with iz[%d] = %d", zz, n);
                        }
                    }
                    else if ((gr->ra[i].iz[zz] == -1) && (gr->ra[i].zeta[zz] != 0))
                    {
                        gmx_fatal(FARGS, "ra[%d].iz[%d] = %d whereas ra[%d].zeta[%d] = %g", i, zz, gr->ra[i].iz[zz], i, zz, gr->ra[i].zeta[zz]);
                    }
                }
            }
        }
    }
    if (NULL != rnd)
    {
        gmx_rng_destroy(rnd);
    }
    if (n != gr->nparam)
    {
        gmx_fatal(FARGS, "Whoopsydaisies! n = %d, should be %d. bSet = %d", n, gr->nparam, bSet);
    }

    if (nrest > 0)
    {
        dq = (gr->qtot-qtot)/nrest;
        if (debug)
        {
            fprintf(debug, "gr->qtot = %g, qtot = %g, nrest = %d, dq = %g\n",
                    gr->qtot, qtot, nrest, dq);
        }
        for (i = 0; (i < gr->natom); i++)
        {
            if (gr->ra[i].iq[gr->ra[i].nZeta-1] == -1)
            {
                gr->ra[i].q[gr->ra[i].nZeta-1] = dq;
            }
        }
    }
    /* Check for excessive charges */
    for (i = 0; (i < gr->natom); i++)
    {
        qi = 0;
        for (zz = 0; (zz < gr->ra[i].nZeta); zz++)
        {
            qi += gr->ra[i].q[zz];
        }
        if (qi < gr->qmin)
        {
            gr->penalty += sqr(gr->qmin-qi);
        }
        else if (qi > gr->qmax)
        {
            gr->penalty += sqr(gr->qmax-qi);
        }
        else if ((qi < -0.02) && (gr->ra[i].atomnumber == 1))
        {
            gr->penalty += qi*qi;
        }
    }
    gr->penalty *= gr->pfac;
}

void gmx_resp_add_point(gmx_resp_t gr, double x, double y,
                        double z, double V)
{
    int i;

    i = gr->nesp++;
    srenew(gr->esp, gr->nesp);
    srenew(gr->pot, gr->nesp);
    srenew(gr->pot_calc, gr->nesp);
    gr->esp[i][XX]  = x;
    gr->esp[i][YY]  = y;
    gr->esp[i][ZZ]  = z;
    gr->pot[i]      = V;
    gr->pot_calc[i] = 0;
}

static real my_weight(gmx_resp_t gr, int iatom)
{
    if (iatom < gr->natom)
    {
        return gr->watoms;
    }
    else
    {
        return 1.0;
    }
}

void gmx_resp_pot_lsq(gmx_resp_t gr, gmx_stats_t lsq)
{
    int    i;
    double w;

    for (i = 0; (i < gr->nesp); i++)
    {
        w = my_weight(gr, i);
        if (w > 0)
        {
            gmx_stats_add_point(lsq,
                                gmx2convert(gr->pot[i], eg2cHartree_e),
                                gmx2convert(gr->pot_calc[i], eg2cHartree_e), 0, 0);
        }
    }
}

void gmx_resp_calc_rms(gmx_resp_t gr)
{
    int    i;
    double pot2, s2, sum2, w, wtot, entropy;
    char   buf[STRLEN];

    pot2 = sum2 = wtot = entropy = 0;
    sprintf(buf, " - weight %g in fit", gr->watoms);
    for (i = 0; (i < gr->nesp); i++)
    {
        w = my_weight(gr, i);
        if ((NULL != debug) && (i < 2*gr->natom))
        {
            fprintf(debug, "ESP %d QM: %g EEM: %g DIFF: %g%s\n",
                    i, gr->pot[i], gr->pot_calc[i],
                    gr->pot[i]-gr->pot_calc[i],
                    (i < gr->natom)  ? buf : "");
        }
        s2    = w*sqr(gr->pot[i]-gr->pot_calc[i]);
        if ((s2 > 0) && (gr->bEntropy))
        {
            entropy += s2*log(s2);
        }
        sum2 += s2;
        pot2 += w*sqr(gr->pot[i]);
        wtot += w;
    }
    gr->wtot = wtot;
    if (wtot > 0)
    {
        gr->rms     = gmx2convert(sqrt(sum2/wtot), eg2cHartree_e);
        gr->entropy = gmx2convert(entropy/wtot, eg2cHartree_e);
    }
    else
    {
        gr->rms     = 0;
        gr->entropy = 0;
    }
    gr->rrms = sqrt(sum2/pot2);
}

double gmx_resp_get_rms(gmx_resp_t gr, real *wtot)
{
    gmx_resp_calc_rms(gr);
    *wtot = gr->wtot;
    if (gr->bEntropy)
    {
        return gr->entropy;
    }
    else
    {
        return gr->rms;
    }
}

static void gmx_resp_calc_penalty(gmx_resp_t gr)
{
    int    i;
    double p, b2;

    p = 0;
    if (gr->bAXpRESP && (gr->iModel == eqgAXp))
    {
        b2 = sqr(gr->b_hyper);
        for (i = 0; (i < gr->natom); i++)
        {
            p += sqrt(sqr(gr->ra[i].q[0]) + b2) - gr->b_hyper;
        }
        p = (gr->qfac * p);
    }
    gr->penalty += p;
}

static double charge_function(void *params, double v[])
{
    gmx_resp_t gr  = (gmx_resp_t) params;
    double     rms = 0;
    real       wtot;

    get_set_vector(gr, false, false, false, gr->seed, v);
    gmx_resp_calc_pot(gr);
    gmx_resp_calc_penalty(gr);
    rms = gmx_resp_get_rms(gr, &wtot);

    return rms; // + gr->penalty;
}

void gmx_resp_statistics(gmx_resp_t gr, int len, char buf[])
{
    if (len >= 100)
    {
        sprintf(buf, "RMS: %10e [Hartree/e] RRMS: %10e Entropy: %10e Penalty: %10e",
                gr->rms, gr->rrms, gr->entropy, gr->penalty);
    }
    else
    {
        fprintf(stderr, "buflen too small (%d) in gmx_resp_statistics\n", len);
    }
}

int gmx_resp_optimize_charges(FILE *fp, gmx_resp_t gr, int maxiter,
                              real toler, real *rms)
{
    double *param;
    double  ccc;
    int     bConv;
    char    buf[STRLEN];

    snew(param, gr->nparam);

    get_set_vector(gr, true, gr->bRandQ, gr->bRandZeta, gr->seed, param);

    bConv = nmsimplex(fp, (void *)gr, charge_function, param, gr->nparam,
                      toler, 1, maxiter, &ccc);
    if (bConv)
    {
        gmx_resp_statistics(gr, STRLEN-1, buf);
    }
    else
    {
        printf("NM Simplex did not converge\n\n");
    }

    if (gr->bEntropy)
    {
        *rms = gr->entropy;
    }
    else
    {
        *rms = gr->rms;
    }

    get_set_vector(gr, false, false, false, gr->seed, param);

    sfree(param);

    if (bConv)
    {
        return eQGEN_OK;
    }
    else
    {
        return eQGEN_NOTCONVERGED;
    }
}

void gmx_resp_destroy(gmx_resp_t gr)
{
    int i;

    sfree(gr->x);
    sfree(gr->esp);
    sfree(gr->pot);
    i = 0;
    while (NULL != gr->dzatoms[i])
    {
        sfree(gr->dzatoms[i]);
        i++;
    }
    if (NULL != gr->dzatoms)
    {
        sfree(gr->dzatoms);
    }
    for (i = 0; (i < gr->natom); i++)
    {
        gmx_ra_done(&(gr->ra[i]));
    }
    sfree(gr->ra);
}

void gmx_resp_potcomp(gmx_resp_t gr, const char *potcomp,
                      const char *pdbdiff, output_env_t oenv)
{
    int     i;
    double  pp, exp, eem;
    FILE   *fp;
    int     unit = eg2cHartree_e;

    if (NULL != potcomp)
    {
        const char *pcleg[2] = { "Atoms", "ESP points" };
        fp = xvgropen(potcomp, "Electrostatic potential", unit2string(unit), unit2string(unit), oenv);
        xvgr_legend(fp, 2, pcleg, oenv);
        fprintf(fp, "@type xy\n");
        for (i = 0; (i < gr->nesp); i++)
        {
            /* Conversion may or may not be in vain depending on unit */
            exp = gmx2convert(gr->pot[i], unit);
            eem = gmx2convert(gr->pot_calc[i], unit);
            if (i == gr->natom)
            {
                fprintf(fp, "&\n");
                fprintf(fp, "@type xy\n");
            }
            fprintf(fp, "%10.5e  %10.5e\n", exp, eem);
        }
        fprintf(fp, "&\n");
        fclose(fp);
    }
    if (NULL != pdbdiff)
    {
        fp = fopen(pdbdiff, "w");
        fprintf(fp, "REMARK All distances are scaled by a factor of two.\n");
        for (i = 0; (i < gr->nesp); i++)
        {
            exp = gmx2convert(gr->pot[i], eg2cHartree_e);
            eem = gmx2convert(gr->pot_calc[i], eg2cHartree_e);
            pp  = gr->pot[i]-gr->pot_calc[i];
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, "HE", "HE", ' ', i+1, ' ', 20*gr->esp[i][XX],
                    20*gr->esp[i][YY], 20*gr->esp[i][ZZ], 0.0, pp);
        }
        fclose(fp);
    }
}

double gmx_resp_get_qtot(gmx_resp_t grt, int atom)
{
    int    i;
    double q = 0;

    range_check(atom, 0, grt->natom);
    for (i = 0; (i < grt->ra[atom].nZeta); i++)
    {
        q += grt->ra[atom].q[i];
    }
    return q;
}

double gmx_resp_get_q(gmx_resp_t grt, int atom, int zz)
{
    range_check(atom, 0, grt->natom);
    range_check(zz, 0, grt->ra[atom].nZeta);

    return grt->ra[atom].q[zz];
}

double gmx_resp_get_zeta(gmx_resp_t grt, int atom, int zz)
{
    range_check(atom, 0, grt->natom);
    range_check(zz, 0, grt->ra[atom].nZeta);

    return grt->ra[atom].zeta[zz];
}

void gmx_resp_set_q(gmx_resp_t grt, int atom, int zz, double q)
{
    range_check(atom, 0, grt->natom);
    range_check(zz, 0, grt->ra[atom].nZeta);

    grt->ra[atom].q[zz] = q;
}

void gmx_resp_set_zeta(gmx_resp_t grt, int atom, int zz, double zeta)
{
    range_check(atom, 0, grt->natom);
    range_check(zz, 0, grt->ra[atom].nZeta);

    grt->ra[atom].zeta[zz] = zeta;
}
