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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "physics.h"
#include "index.h"
#include "gromacs/utility/smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "gstat.h"
#include "gromacs/fileio/matio.h"
#include "gmx_ana.h"
#include "names.h"

#include "gromacs/legacyheaders/gmx_fatal.h"

static void check_box_c(matrix box)
{
    if (fabs(box[ZZ][XX]) > GMX_REAL_EPS*box[ZZ][ZZ] ||
        fabs(box[ZZ][YY]) > GMX_REAL_EPS*box[ZZ][ZZ])
    {
        gmx_fatal(FARGS,
                  "The last box vector is not parallel to the z-axis: %f %f %f",
                  box[ZZ][XX], box[ZZ][YY], box[ZZ][ZZ]);
    }
}

static void calc_comg(int is, int *coi, int *index, gmx_bool bMass, t_atom *atom,
                      rvec *x, rvec *x_comg)
{
    int  c, i, d;
    rvec xc;
    real mtot, m;

    if (bMass && atom == NULL)
    {
        gmx_fatal(FARGS, "No masses available while mass weighting was requested");
    }

    for (c = 0; c < is; c++)
    {
        clear_rvec(xc);
        mtot = 0;
        for (i = coi[c]; i < coi[c+1]; i++)
        {
            if (bMass)
            {
                m = atom[index[i]].m;
                for (d = 0; d < DIM; d++)
                {
                    xc[d] += m*x[index[i]][d];
                }
                mtot += m;
            }
            else
            {
                rvec_inc(xc, x[index[i]]);
                mtot += 1.0;
            }
        }
        svmul(1/mtot, xc, x_comg[c]);
    }
}

static void split_group(int isize, int *index, char *grpname,
                        t_topology *top, char type,
                        int *is_out, int **coi_out)
{
    t_block *mols = NULL;
    t_atom  *atom = NULL;
    int      is, *coi;
    int      cur, mol, res, i, a, i1;

    /* Split up the group in molecules or residues */
    switch (type)
    {
        case 'm':
            mols = &top->mols;
            break;
        case 'r':
            atom = top->atoms.atom;
            break;
        default:
            gmx_fatal(FARGS, "Unknown rdf option '%s'", type);
    }
    snew(coi, isize+1);
    is  = 0;
    cur = -1;
    mol = 0;
    for (i = 0; i < isize; i++)
    {
        a = index[i];
        if (type == 'm')
        {
            /* Check if the molecule number has changed */
            i1 = mols->index[mol+1];
            while (a >= i1)
            {
                mol++;
                i1 = mols->index[mol+1];
            }
            if (mol != cur)
            {
                coi[is++] = i;
                cur       = mol;
            }
        }
        else if (type == 'r')
        {
            /* Check if the residue index has changed */
            res = atom[a].resind;
            if (res != cur)
            {
                coi[is++] = i;
                cur       = res;
            }
        }
    }
    coi[is] = i;
    srenew(coi, is+1);
    printf("Group '%s' of %d atoms consists of %d %s\n",
           grpname, isize, is,
           (type == 'm' ? "molecules" : "residues"));

    *is_out  = is;
    *coi_out = coi;
}

static void do_rdf(const char *fnNDX, const char *fnTPS, const char *fnTRX,
                   const char *fnRDF, const char *fnCNRDF, const char *fnHQ,
                   gmx_bool bCM, const char *close,
                   const char **rdft, gmx_bool bXY, gmx_bool bPBC, gmx_bool bNormalize,
                   real cutoff, real binwidth, real fade, int ng,
                   const output_env_t oenv)
{
    FILE          *fp;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, i, ii, j, k, nbin, j0, j1, n, nframes;
    int          **count;
    char         **grpname;
    int           *isize, isize_cm = 0, nrdf = 0, max_i, isize0, isize_g;
    atom_id      **index, *index_cm = NULL;
    gmx_int64_t   *sum;
    real           t, rmax2, cut2, r, r2, r2ii, invhbinw, normfac;
    real           segvol, spherevol, prev_spherevol, **rdf;
    rvec          *x, dx, *x0 = NULL, *x_i1, xi;
    real          *inv_segvol, invvol, invvol_sum, rho;
    gmx_bool       bClose, *bExcl, bTop, bNonSelfExcl;
    matrix         box, box_pbc;
    int          **npairs;
    atom_id        ix, jx, ***pairs;
    t_topology    *top  = NULL;
    int            ePBC = -1, ePBCrdf = -1;
    t_block       *mols = NULL;
    t_blocka      *excl;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    int           *is   = NULL, **coi = NULL, cur, mol, i1, res, a;

    excl = NULL;

    bClose = (close[0] != 'n');

    if (fnTPS)
    {
        snew(top, 1);
        bTop = read_tps_conf(fnTPS, title, top, &ePBC, &x, NULL, box, TRUE);
        if (bTop && !(bCM || bClose))
        {
            /* get exclusions from topology */
            excl = &(top->excls);
        }
    }
    snew(grpname, ng+1);
    snew(isize, ng+1);
    snew(index, ng+1);
    fprintf(stderr, "\nSelect a reference group and %d group%s\n",
            ng, ng == 1 ? "" : "s");
    if (fnTPS)
    {
        get_index(&(top->atoms), fnNDX, ng+1, isize, index, grpname);
        atom = top->atoms.atom;
    }
    else
    {
        rd_index(fnNDX, ng+1, isize, index, grpname);
    }

    if (bCM || bClose)
    {
        snew(is, 1);
        snew(coi, 1);
        if (bClose)
        {
            split_group(isize[0], index[0], grpname[0], top, close[0], &is[0], &coi[0]);
        }
    }
    if (rdft[0][0] != 'a')
    {
        /* Split up all the groups in molecules or residues */
        srenew(is, ng+1);
        srenew(coi, ng+1);
        for (g = ((bCM || bClose) ? 1 : 0); g < ng+1; g++)
        {
            split_group(isize[g], index[g], grpname[g], top, rdft[0][0], &is[g], &coi[g]);
        }
    }

    if (bCM)
    {
        is[0] = 1;
        snew(coi[0], is[0]+1);
        coi[0][0] = 0;
        coi[0][1] = isize[0];
        isize0    = is[0];
        snew(x0, isize0);
    }
    else if (bClose || rdft[0][0] != 'a')
    {
        isize0 = is[0];
        snew(x0, isize0);
    }
    else
    {
        isize0 = isize[0];
    }

    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (!natoms)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }
    if (fnTPS)
    {
        /* check with topology */
        if (natoms > top->atoms.nr)
        {
            gmx_fatal(FARGS, "Trajectory (%d atoms) does not match topology (%d atoms)",
                      natoms, top->atoms.nr);
        }
    }
    /* check with index groups */
    for (i = 0; i < ng+1; i++)
    {
        for (j = 0; j < isize[i]; j++)
        {
            if (index[i][j] >= natoms)
            {
                gmx_fatal(FARGS, "Atom index (%d) in index group %s (%d atoms) larger "
                          "than number of atoms in trajectory (%d atoms)",
                          index[i][j], grpname[i], isize[i], natoms);
            }
        }
    }

    /* initialize some handy things */
    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }
    copy_mat(box, box_pbc);
    if (bXY)
    {
        check_box_c(box);
        switch (ePBC)
        {
            case epbcXYZ:
            case epbcXY:   ePBCrdf = epbcXY;   break;
            case epbcNONE: ePBCrdf = epbcNONE; break;
            default:
                gmx_fatal(FARGS, "xy-rdf's are not supported for pbc type'%s'",
                          EPBC(ePBC));
                break;
        }
        /* Make sure the z-height does not influence the cut-off */
        box_pbc[ZZ][ZZ] = 2*max(box[XX][XX], box[YY][YY]);
    }
    else
    {
        ePBCrdf = ePBC;
    }
    if (bPBC)
    {
        rmax2   = 0.99*0.99*max_cutoff2(bXY ? epbcXY : epbcXYZ, box_pbc);
    }
    else
    {
        rmax2   = sqr(3*max(box[XX][XX], max(box[YY][YY], box[ZZ][ZZ])));
    }
    if (debug)
    {
        fprintf(debug, "rmax2 = %g\n", rmax2);
    }

    /* We use the double amount of bins, so we can correctly
     * write the rdf and rdf_cn output at i*binwidth values.
     */
    nbin     = (int)(sqrt(rmax2) * 2 / binwidth);
    invhbinw = 2.0 / binwidth;
    cut2     = sqr(cutoff);

    snew(count, ng);
    snew(pairs, ng);
    snew(npairs, ng);

    snew(bExcl, natoms);
    max_i = 0;
    for (g = 0; g < ng; g++)
    {
        if (isize[g+1] > max_i)
        {
            max_i = isize[g+1];
        }

        /* this is THE array */
        snew(count[g], nbin+1);

        /* make pairlist array for groups and exclusions */
        snew(pairs[g], isize[0]);
        snew(npairs[g], isize[0]);
        for (i = 0; i < isize[0]; i++)
        {
            /* We can only have exclusions with atomic rdfs */
            if (!(bCM || bClose || rdft[0][0] != 'a'))
            {
                ix = index[0][i];
                for (j = 0; j < natoms; j++)
                {
                    bExcl[j] = FALSE;
                }
                /* exclusions? */
                if (excl)
                {
                    for (j = excl->index[ix]; j < excl->index[ix+1]; j++)
                    {
                        bExcl[excl->a[j]] = TRUE;
                    }
                }
                k = 0;
                snew(pairs[g][i], isize[g+1]);
                bNonSelfExcl = FALSE;
                for (j = 0; j < isize[g+1]; j++)
                {
                    jx = index[g+1][j];
                    if (!bExcl[jx])
                    {
                        pairs[g][i][k++] = jx;
                    }
                    else if (ix != jx)
                    {
                        /* Check if we have exclusions other than self exclusions */
                        bNonSelfExcl = TRUE;
                    }
                }
                if (bNonSelfExcl)
                {
                    npairs[g][i] = k;
                    srenew(pairs[g][i], npairs[g][i]);
                }
                else
                {
                    /* Save a LOT of memory and some cpu cycles */
                    npairs[g][i] = -1;
                    sfree(pairs[g][i]);
                }
            }
            else
            {
                npairs[g][i] = -1;
            }
        }
    }
    sfree(bExcl);

    snew(x_i1, max_i);
    nframes    = 0;
    invvol_sum = 0;
    if (bPBC && (NULL != top))
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    do
    {
        /* Must init pbc every step because of pressure coupling */
        copy_mat(box, box_pbc);
        if (bPBC)
        {
            if (top != NULL)
            {
                gmx_rmpbc(gpbc, natoms, box, x);
            }
            if (bXY)
            {
                check_box_c(box);
                clear_rvec(box_pbc[ZZ]);
            }
            set_pbc(&pbc, ePBCrdf, box_pbc);

            if (bXY)
            {
                /* Set z-size to 1 so we get the surface iso the volume */
                box_pbc[ZZ][ZZ] = 1;
            }
        }
        invvol      = 1/det(box_pbc);
        invvol_sum += invvol;

        if (bCM)
        {
            /* Calculate center of mass of the whole group */
            calc_comg(is[0], coi[0], index[0], TRUE, atom, x, x0);
        }
        else if (!bClose && rdft[0][0] != 'a')
        {
            calc_comg(is[0], coi[0], index[0], rdft[0][6] == 'm', atom, x, x0);
        }

        for (g = 0; g < ng; g++)
        {
            if (rdft[0][0] == 'a')
            {
                /* Copy the indexed coordinates to a continuous array */
                for (i = 0; i < isize[g+1]; i++)
                {
                    copy_rvec(x[index[g+1][i]], x_i1[i]);
                }
            }
            else
            {
                /* Calculate the COMs/COGs and store in x_i1 */
                calc_comg(is[g+1], coi[g+1], index[g+1], rdft[0][6] == 'm', atom, x, x_i1);
            }

            for (i = 0; i < isize0; i++)
            {
                if (bClose)
                {
                    /* Special loop, since we need to determine the minimum distance
                     * over all selected atoms in the reference molecule/residue.
                     */
                    if (rdft[0][0] == 'a')
                    {
                        isize_g = isize[g+1];
                    }
                    else
                    {
                        isize_g = is[g+1];
                    }
                    for (j = 0; j < isize_g; j++)
                    {
                        r2 = 1e30;
                        /* Loop over the selected atoms in the reference molecule */
                        for (ii = coi[0][i]; ii < coi[0][i+1]; ii++)
                        {
                            if (bPBC)
                            {
                                pbc_dx(&pbc, x[index[0][ii]], x_i1[j], dx);
                            }
                            else
                            {
                                rvec_sub(x[index[0][ii]], x_i1[j], dx);
                            }
                            if (bXY)
                            {
                                r2ii = dx[XX]*dx[XX] + dx[YY]*dx[YY];
                            }
                            else
                            {
                                r2ii = iprod(dx, dx);
                            }
                            if (r2ii < r2)
                            {
                                r2 = r2ii;
                            }
                        }
                        if (r2 > cut2 && r2 <= rmax2)
                        {
                            count[g][(int)(sqrt(r2)*invhbinw)]++;
                        }
                    }
                }
                else
                {
                    /* Real rdf between points in space */
                    if (bCM || rdft[0][0] != 'a')
                    {
                        copy_rvec(x0[i], xi);
                    }
                    else
                    {
                        copy_rvec(x[index[0][i]], xi);
                    }
                    if (rdft[0][0] == 'a' && npairs[g][i] >= 0)
                    {
                        /* Expensive loop, because of indexing */
                        for (j = 0; j < npairs[g][i]; j++)
                        {
                            jx = pairs[g][i][j];
                            if (bPBC)
                            {
                                pbc_dx(&pbc, xi, x[jx], dx);
                            }
                            else
                            {
                                rvec_sub(xi, x[jx], dx);
                            }

                            if (bXY)
                            {
                                r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
                            }
                            else
                            {
                                r2 = iprod(dx, dx);
                            }
                            if (r2 > cut2 && r2 <= rmax2)
                            {
                                count[g][(int)(sqrt(r2)*invhbinw)]++;
                            }
                        }
                    }
                    else
                    {
                        /* Cheaper loop, no exclusions */
                        if (rdft[0][0] == 'a')
                        {
                            isize_g = isize[g+1];
                        }
                        else
                        {
                            isize_g = is[g+1];
                        }
                        for (j = 0; j < isize_g; j++)
                        {
                            if (bPBC)
                            {
                                pbc_dx(&pbc, xi, x_i1[j], dx);
                            }
                            else
                            {
                                rvec_sub(xi, x_i1[j], dx);
                            }
                            if (bXY)
                            {
                                r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
                            }
                            else
                            {
                                r2 = iprod(dx, dx);
                            }
                            if (r2 > cut2 && r2 <= rmax2)
                            {
                                count[g][(int)(sqrt(r2)*invhbinw)]++;
                            }
                        }
                    }
                }
            }
        }
        nframes++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");

    if (bPBC && (NULL != top))
    {
        gmx_rmpbc_done(gpbc);
    }

    close_trj(status);

    sfree(x);

    /* Average volume */
    invvol = invvol_sum/nframes;

    /* Calculate volume of sphere segments or length of circle segments */
    snew(inv_segvol, (nbin+1)/2);
    prev_spherevol = 0;
    for (i = 0; (i < (nbin+1)/2); i++)
    {
        r = (i + 0.5)*binwidth;
        if (bXY)
        {
            spherevol = M_PI*r*r;
        }
        else
        {
            spherevol = (4.0/3.0)*M_PI*r*r*r;
        }
        segvol         = spherevol-prev_spherevol;
        inv_segvol[i]  = 1.0/segvol;
        prev_spherevol = spherevol;
    }

    snew(rdf, ng);
    for (g = 0; g < ng; g++)
    {
        /* We have to normalize by dividing by the number of frames */
        if (rdft[0][0] == 'a')
        {
            normfac = 1.0/(nframes*invvol*isize0*isize[g+1]);
        }
        else
        {
            normfac = 1.0/(nframes*invvol*isize0*is[g+1]);
        }

        /* Do the normalization */
        nrdf = max((nbin+1)/2, 1+2*fade/binwidth);
        snew(rdf[g], nrdf);
        for (i = 0; i < (nbin+1)/2; i++)
        {
            r = i*binwidth;
            if (i == 0)
            {
                j = count[g][0];
            }
            else
            {
                j = count[g][i*2-1] + count[g][i*2];
            }
            if ((fade > 0) && (r >= fade))
            {
                rdf[g][i] = 1 + (j*inv_segvol[i]*normfac-1)*exp(-16*sqr(r/fade-1));
            }
            else
            {
                if (bNormalize)
                {
                    rdf[g][i] = j*inv_segvol[i]*normfac;
                }
                else
                {
                    rdf[g][i] = j/(binwidth*isize0*nframes);
                }
            }
        }
        for (; (i < nrdf); i++)
        {
            rdf[g][i] = 1.0;
        }
    }

    if (rdft[0][0] == 'a')
    {
        sprintf(gtitle, "Radial distribution");
    }
    else
    {
        sprintf(gtitle, "Radial distribution of %s %s",
                rdft[0][0] == 'm' ? "molecule" : "residue",
                rdft[0][6] == 'm' ? "COM" : "COG");
    }
    fp = xvgropen(fnRDF, gtitle, "r", "", oenv);
    if (bCM)
    {
        sprintf(refgt, " %s", "COM");
    }
    else if (bClose)
    {
        sprintf(refgt, " closest atom in %s.", close);
    }
    else
    {
        sprintf(refgt, "%s", "");
    }
    if (ng == 1)
    {
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fp, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
        }
    }
    else
    {
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fp, "@ subtitle \"reference %s%s\"\n", grpname[0], refgt);
        }
        xvgr_legend(fp, ng, (const char**)(grpname+1), oenv);
    }
    for (i = 0; (i < nrdf); i++)
    {
        fprintf(fp, "%10g", i*binwidth);
        for (g = 0; g < ng; g++)
        {
            fprintf(fp, " %10g", rdf[g][i]);
        }
        fprintf(fp, "\n");
    }
    gmx_ffclose(fp);

    do_view(oenv, fnRDF, NULL);

    /* h(Q) function: fourier transform of rdf */
    if (fnHQ)
    {
        int   nhq = 401;
        real *hq, *integrand, Q;

        /* Get a better number density later! */
        rho = isize[1]*invvol;
        snew(hq, nhq);
        snew(integrand, nrdf);
        for (i = 0; (i < nhq); i++)
        {
            Q            = i*0.5;
            integrand[0] = 0;
            for (j = 1; (j < nrdf); j++)
            {
                r             = j*binwidth;
                integrand[j]  = (Q == 0) ? 1.0 : sin(Q*r)/(Q*r);
                integrand[j] *= 4.0*M_PI*rho*r*r*(rdf[0][j]-1.0);
            }
            hq[i] = print_and_integrate(debug, nrdf, binwidth, integrand, NULL, 0);
        }
        fp = xvgropen(fnHQ, "h(Q)", "Q(/nm)", "h(Q)", oenv);
        for (i = 0; (i < nhq); i++)
        {
            fprintf(fp, "%10g %10g\n", i*0.5, hq[i]);
        }
        gmx_ffclose(fp);
        do_view(oenv, fnHQ, NULL);
        sfree(hq);
        sfree(integrand);
    }

    if (fnCNRDF)
    {
        normfac = 1.0/(isize0*nframes);
        fp      = xvgropen(fnCNRDF, "Cumulative Number RDF", "r", "number", oenv);
        if (ng == 1)
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"%s-%s\"\n", grpname[0], grpname[1]);
            }
        }
        else
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"reference %s\"\n", grpname[0]);
            }
            xvgr_legend(fp, ng, (const char**)(grpname+1), oenv);
        }
        snew(sum, ng);
        for (i = 0; (i <= nbin/2); i++)
        {
            fprintf(fp, "%10g", i*binwidth);
            for (g = 0; g < ng; g++)
            {
                fprintf(fp, " %10g", (real)((double)sum[g]*normfac));
                if (i*2+1 < nbin)
                {
                    sum[g] += count[g][i*2] + count[g][i*2+1];
                }
            }
            fprintf(fp, "\n");
        }
        gmx_ffclose(fp);
        sfree(sum);

        do_view(oenv, fnCNRDF, NULL);
    }

    for (g = 0; g < ng; g++)
    {
        sfree(rdf[g]);
    }
    sfree(rdf);
}


int gmx_rdf(int argc, char *argv[])
{
    const char        *desc[] = {
        "The structure of liquids can be studied by either neutron or X-ray",
        "scattering. The most common way to describe liquid structure is by a",
        "radial distribution function. However, this is not easy to obtain from",
        "a scattering experiment.[PAR]",
        "[THISMODULE] calculates radial distribution functions in different ways.",
        "The normal method is around a (set of) particle(s), the other methods",
        "are around the center of mass of a set of particles ([TT]-com[tt])",
        "or to the closest particle in a set ([TT]-surf[tt]).",
        "With all methods, the RDF can also be calculated around axes parallel",
        "to the [IT]z[it]-axis with option [TT]-xy[tt].",
        "With option [TT]-surf[tt] normalization can not be used.[PAR]",
        "The option [TT]-rdf[tt] sets the type of RDF to be computed.",
        "Default is for atoms or particles, but one can also select center",
        "of mass or geometry of molecules or residues. In all cases, only",
        "the atoms in the index groups are taken into account.",
        "For molecules and/or the center of mass option, a run input file",
        "is required.",
        "Weighting other than COM or COG can currently only be achieved",
        "by providing a run input file with different masses.",
        "Options [TT]-com[tt] and [TT]-surf[tt] also work in conjunction",
        "with [TT]-rdf[tt].[PAR]",
        "If a run input file is supplied ([TT]-s[tt]) and [TT]-rdf[tt] is set",
        "to [TT]atom[tt], exclusions defined",
        "in that file are taken into account when calculating the RDF.",
        "The option [TT]-cut[tt] is meant as an alternative way to avoid",
        "intramolecular peaks in the RDF plot.",
        "It is however better to supply a run input file with a higher number of",
        "exclusions. For e.g. benzene a topology, setting nrexcl to 5",
        "would eliminate all intramolecular contributions to the RDF.",
        "Note that all atoms in the selected groups are used, also the ones",
        "that don't have Lennard-Jones interactions.[PAR]",
        "Option [TT]-cn[tt] produces the cumulative number RDF,",
        "i.e. the average number of particles within a distance r.[PAR]"
    };
    static gmx_bool    bCM     = FALSE, bXY = FALSE, bPBC = TRUE, bNormalize = TRUE;
    static real        cutoff  = 0, binwidth = 0.002, fade = 0.0;
    static int         ngroups = 1;

    static const char *closet[] = { NULL, "no", "mol", "res", NULL };
    static const char *rdft[]   = { NULL, "atom", "mol_com", "mol_cog", "res_com", "res_cog", NULL };

    t_pargs            pa[] = {
        { "-bin",      FALSE, etREAL, {&binwidth},
          "Binwidth (nm)" },
        { "-com",      FALSE, etBOOL, {&bCM},
          "RDF with respect to the center of mass of first group" },
        { "-surf",     FALSE, etENUM, {closet},
          "RDF with respect to the surface of the first group" },
        { "-rdf",   FALSE, etENUM, {rdft},
          "RDF type" },
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-norm",     FALSE, etBOOL, {&bNormalize},
          "Normalize for volume and density" },
        { "-xy",       FALSE, etBOOL, {&bXY},
          "Use only the x and y components of the distance" },
        { "-cut",      FALSE, etREAL, {&cutoff},
          "Shortest distance (nm) to be considered"},
        { "-ng",       FALSE, etINT, {&ngroups},
          "Number of secondary groups to compute RDFs around a central group" },
        { "-fade",     FALSE, etREAL, {&fade},
          "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/fade-1)^2 to make it go to 1 smoothly. If fade is 0.0 nothing is done." }
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX;
    output_env_t       oenv;

    t_filenm           fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD },
        { efTPS, NULL,  NULL,     ffOPTRD },
        { efNDX, NULL,  NULL,     ffOPTRD },
        { efXVG, "-o",  "rdf",    ffWRITE },
        { efXVG, "-cn", "rdf_cn", ffOPTWR },
        { efXVG, "-hq", "hq",     ffOPTWR },
    };
#define NFILE asize(fnm)
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    if (bCM || closet[0][0] != 'n' || rdft[0][0] == 'm' || rdft[0][6] == 'm')
    {
        fnTPS = ftp2fn(efTPS, NFILE, fnm);
    }
    else
    {
        fnTPS = ftp2fn_null(efTPS, NFILE, fnm);
    }
    fnNDX = ftp2fn_null(efNDX, NFILE, fnm);

    if (!fnTPS && !fnNDX)
    {
        gmx_fatal(FARGS, "Neither index file nor topology file specified\n"
                  "Nothing to do!");
    }

    if (closet[0][0] != 'n')
    {
        if (bCM)
        {
            gmx_fatal(FARGS, "Can not have both -com and -surf");
        }
        if (bNormalize)
        {
            fprintf(stderr, "Turning of normalization because of option -surf\n");
            bNormalize = FALSE;
        }
    }

    do_rdf(fnNDX, fnTPS, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn_null("-cn", NFILE, fnm),
           opt2fn_null("-hq", NFILE, fnm),
           bCM, closet[0], rdft, bXY, bPBC, bNormalize, cutoff, binwidth, fade, ngroups,
           oenv);

    return 0;
}
