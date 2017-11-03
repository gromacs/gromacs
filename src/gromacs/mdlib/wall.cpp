/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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

#include <string.h>

#include <algorithm>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void make_wall_tables(FILE *fplog,
                      const t_inputrec *ir, const char *tabfn,
                      const gmx_groups_t *groups,
                      t_forcerec *fr)
{
    int           negp_pp;
    int          *nm_ind;
    char          buf[STRLEN];

    negp_pp = ir->opts.ngener - ir->nwall;
    nm_ind  = groups->grps[egcENER].nm_ind;

    if (fplog)
    {
        fprintf(fplog, "Reading user tables for %d energy groups with %d walls\n",
                negp_pp, ir->nwall);
    }

    snew(fr->wall_tab, ir->nwall);
    for (int w = 0; w < ir->nwall; w++)
    {
        snew(fr->wall_tab[w], negp_pp);
        for (int egp = 0; egp < negp_pp; egp++)
        {
            /* If the energy group pair is excluded, we don't need a table */
            if (!(fr->egp_flags[egp*ir->opts.ngener+negp_pp+w] & EGP_EXCL))
            {
                sprintf(buf, "%s", tabfn);
                sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1, "_%s_%s.%s",
                        *groups->grpname[nm_ind[egp]],
                        *groups->grpname[nm_ind[negp_pp+w]],
                        ftp2ext(efXVG));
                fr->wall_tab[w][egp] = make_tables(fplog, fr->ic, buf, 0,
                                                   GMX_MAKETABLES_FORCEUSER);

                /* Since wall have no charge, we can compress the table */
                for (int i = 0; i <= fr->wall_tab[w][egp]->n; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        fr->wall_tab[w][egp]->data[8*i+j] =
                            fr->wall_tab[w][egp]->data[12*i+4+j];
                    }
                }
            }
        }
    }
}

static void wall_error(int a, rvec *x, real r)
{
    gmx_fatal(FARGS,
              "An atom is beyond the wall: coordinates %f %f %f, distance %f\n"
              "You might want to use the mdp option wall_r_linpot",
              x[a][XX], x[a][YY], x[a][ZZ], r);
}

real do_walls(t_inputrec *ir, t_forcerec *fr, matrix box, t_mdatoms *md,
              rvec x[], rvec f[], real lambda, real Vlj[], t_nrnb *nrnb)
{
    int             nwall;
    int             ntw[2], at, ntype, ngid, ggid, *egp_flags, *type;
    real           *nbfp, lamfac, fac_d[2], fac_r[2], Cd, Cr;
    real            wall_z[2], r, mr, r1, r2, r4, Vd, Vr, V = 0, Fd, Fr, F = 0, dvdlambda;
    dvec            xf_z;
    int             n0, nnn;
    real            tabscale, *VFtab, rt, eps, eps2, Yt, Ft, Geps, Heps2, Fp, VV, FF;
    unsigned short *gid = md->cENER;
    t_forcetable   *tab;

    nwall     = ir->nwall;
    ngid      = ir->opts.ngener;
    ntype     = fr->ntype;
    nbfp      = fr->nbfp;
    egp_flags = fr->egp_flags;

    for (int w = 0; w < nwall; w++)
    {
        ntw[w] = 2*ntype*ir->wall_atomtype[w];
        switch (ir->wall_type)
        {
            case ewt93:
                fac_d[w] = ir->wall_density[w]*M_PI/6;
                fac_r[w] = ir->wall_density[w]*M_PI/45;
                break;
            case ewt104:
                fac_d[w] = ir->wall_density[w]*M_PI/2;
                fac_r[w] = ir->wall_density[w]*M_PI/5;
                break;
            default:
                break;
        }
    }
    wall_z[0] = 0;
    wall_z[1] = box[ZZ][ZZ];

    dvdlambda = 0;
    clear_dvec(xf_z);
    for (int lam = 0; lam < (md->nPerturbed ? 2 : 1); lam++)
    {
        if (md->nPerturbed)
        {
            if (lam == 0)
            {
                lamfac = 1 - lambda;
                type   = md->typeA;
            }
            else
            {
                lamfac = lambda;
                type   = md->typeB;
            }
        }
        else
        {
            lamfac = 1;
            type   = md->typeA;
        }

        real Vlambda = 0;
        for (int i = 0; i < md->homenr; i++)
        {
            for (int w = 0; w < std::min(nwall, 2); w++)
            {
                /* The wall energy groups are always at the end of the list */
                ggid = gid[i]*ngid + ngid - nwall + w;
                at   = type[i];
                /* nbfp now includes the 6/12 derivative prefactors */
                Cd = nbfp[ntw[w]+2*at]/6;
                Cr = nbfp[ntw[w]+2*at+1]/12;
                if (!((Cd == 0 && Cr == 0) || (egp_flags[ggid] & EGP_EXCL)))
                {
                    if (w == 0)
                    {
                        r = x[i][ZZ];
                    }
                    else
                    {
                        r = wall_z[1] - x[i][ZZ];
                    }
                    if (r < ir->wall_r_linpot)
                    {
                        mr = ir->wall_r_linpot - r;
                        r  = ir->wall_r_linpot;
                    }
                    else
                    {
                        mr = 0;
                    }
                    switch (ir->wall_type)
                    {
                        case ewtTABLE:
                            if (r < 0)
                            {
                                wall_error(i, x, r);
                            }
                            tab      = fr->wall_tab[w][gid[i]];
                            tabscale = tab->scale;
                            VFtab    = tab->data;

                            rt    = r*tabscale;
                            n0    = static_cast<int>(rt);
                            if (n0 >= tab->n)
                            {
                                /* Beyond the table range, set V and F to zero */
                                V     = 0;
                                F     = 0;
                            }
                            else
                            {
                                eps   = rt - n0;
                                eps2  = eps*eps;
                                /* Dispersion */
                                nnn   = 8*n0;
                                Yt    = VFtab[nnn];
                                Ft    = VFtab[nnn+1];
                                Geps  = VFtab[nnn+2]*eps;
                                Heps2 = VFtab[nnn+3]*eps2;
                                Fp    = Ft + Geps + Heps2;
                                VV    = Yt + Fp*eps;
                                FF    = Fp + Geps + 2.0*Heps2;
                                Vd    = 6*Cd*VV;
                                Fd    = 6*Cd*FF;
                                /* Repulsion */
                                nnn   = nnn + 4;
                                Yt    = VFtab[nnn];
                                Ft    = VFtab[nnn+1];
                                Geps  = VFtab[nnn+2]*eps;
                                Heps2 = VFtab[nnn+3]*eps2;
                                Fp    = Ft + Geps + Heps2;
                                VV    = Yt + Fp*eps;
                                FF    = Fp + Geps + 2.0*Heps2;
                                Vr    = 12*Cr*VV;
                                Fr    = 12*Cr*FF;
                                V     = Vd + Vr;
                                F     = -lamfac*(Fd + Fr)*tabscale;
                            }
                            break;
                        case ewt93:
                            if (r <= 0)
                            {
                                wall_error(i, x, r);
                            }
                            r1 = 1/r;
                            r2 = r1*r1;
                            r4 = r2*r2;
                            Vd = fac_d[w]*Cd*r2*r1;
                            Vr = fac_r[w]*Cr*r4*r4*r1;
                            V  = Vr - Vd;
                            F  = lamfac*(9*Vr - 3*Vd)*r1;
                            break;
                        case ewt104:
                            if (r <= 0)
                            {
                                wall_error(i, x, r);
                            }
                            r1 = 1/r;
                            r2 = r1*r1;
                            r4 = r2*r2;
                            Vd = fac_d[w]*Cd*r4;
                            Vr = fac_r[w]*Cr*r4*r4*r2;
                            V  = Vr - Vd;
                            F  = lamfac*(10*Vr - 4*Vd)*r1;
                            break;
                        case ewt126:
                            if (r <= 0)
                            {
                                wall_error(i, x, r);
                            }
                            r1 = 1/r;
                            r2 = r1*r1;
                            r4 = r2*r2;
                            Vd = Cd*r4*r2;
                            Vr = Cr*r4*r4*r4;
                            V  = Vr - Vd;
                            F  = lamfac*(12*Vr - 6*Vd)*r1;
                            break;
                        default:
                            break;
                    }
                    if (mr > 0)
                    {
                        V += mr*F;
                    }
                    if (w == 1)
                    {
                        F = -F;
                    }
                    Vlj[ggid] += lamfac*V;
                    Vlambda   += V;
                    f[i][ZZ]  += F;
                    /* Because of the single sum virial calculation we need
                     * to add  the full virial contribution of the walls.
                     * Since the force only has a z-component, there is only
                     * a contribution to the z component of the virial tensor.
                     * We could also determine the virial contribution directly,
                     * which would be cheaper here, but that would require extra
                     * communication for f_novirsum for with virtual sites
                     * in parallel.
                     */
                    xf_z[XX]  -= x[i][XX]*F;
                    xf_z[YY]  -= x[i][YY]*F;
                    xf_z[ZZ]  -= wall_z[w]*F;
                }
            }
        }
        if (md->nPerturbed)
        {
            dvdlambda += (lam == 0 ? -1 : 1)*Vlambda;
        }

        inc_nrnb(nrnb, eNR_WALLS, md->homenr);
    }

    for (int i = 0; i < DIM; i++)
    {
        fr->vir_wall_z[i] = -0.5*xf_z[i];
    }

    return dvdlambda;
}
