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

#include "wall.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <vector>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void make_wall_tables(FILE*                   fplog,
                      const t_inputrec&       ir,
                      const char*             tabfn,
                      const SimulationGroups* groups,
                      t_forcerec*             fr)
{
    int  negp_pp;
    char buf[STRLEN];

    negp_pp                         = ir.opts.ngener - ir.nwall;
    gmx::ArrayRef<const int> nm_ind = groups->groups[SimulationAtomGroupType::EnergyOutput];

    if (fplog)
    {
        fprintf(fplog, "Reading user tables for %d energy groups with %d walls\n", negp_pp, ir.nwall);
    }

    fr->wall_tab.resize(ir.nwall);
    for (int w = 0; w < ir.nwall; w++)
    {
        fr->wall_tab[w].resize(negp_pp);
        for (int egp = 0; egp < negp_pp; egp++)
        {
            /* If the energy group pair is excluded, we don't need a table */
            if (!(fr->egp_flags[egp * ir.opts.ngener + negp_pp + w] & EGP_EXCL))
            {
                sprintf(buf, "%s", tabfn);
                sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1,
                        "_%s_%s.%s",
                        *groups->groupNames[nm_ind[egp]],
                        *groups->groupNames[nm_ind[negp_pp + w]],
                        ftp2ext(efXVG));
                fr->wall_tab[w][egp] = make_tables(fplog, fr->ic.get(), buf, 0, GMX_MAKETABLES_FORCEUSER);

                /* Since wall have no charge, we can compress the table */
                for (int i = 0; i <= fr->wall_tab[w][egp]->numTablePoints; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        fr->wall_tab[w][egp]->data[8 * i + j] =
                                fr->wall_tab[w][egp]->data[12 * i + 4 + j];
                    }
                }
            }
        }
    }
}

[[noreturn]] static void wall_error(int a, gmx::ArrayRef<const gmx::RVec> x, real r)
{
    gmx_fatal(FARGS,
              "An atom is beyond the wall: coordinates %f %f %f, distance %f\n"
              "You might want to use the mdp option wall_r_linpot",
              x[a][XX],
              x[a][YY],
              x[a][ZZ],
              r);
}

static void tableForce(real r, const t_forcetable& tab, real Cd, real Cr, real* V, real* F)
{
    const real  tabscale = tab.scale;
    const real* VFtab    = tab.data.data();

    real rt = r * tabscale;
    int  n0 = static_cast<int>(rt);
    if (n0 >= tab.numTablePoints)
    {
        /* Beyond the table range, set V and F to zero */
        *V = 0;
        *F = 0;
    }
    else
    {
        real eps  = rt - n0;
        real eps2 = eps * eps;
        /* Dispersion */
        int  nnn   = 8 * n0;
        real Yt    = VFtab[nnn];
        real Ft    = VFtab[nnn + 1];
        real Geps  = VFtab[nnn + 2] * eps;
        real Heps2 = VFtab[nnn + 3] * eps2;
        real Fp    = Ft + Geps + Heps2;
        real VV    = Yt + Fp * eps;
        real FF    = Fp + Geps + 2.0 * Heps2;
        real Vd    = 6 * Cd * VV;
        real Fd    = 6 * Cd * FF;
        /* Repulsion */
        nnn     = nnn + 4;
        Yt      = VFtab[nnn];
        Ft      = VFtab[nnn + 1];
        Geps    = VFtab[nnn + 2] * eps;
        Heps2   = VFtab[nnn + 3] * eps2;
        Fp      = Ft + Geps + Heps2;
        VV      = Yt + Fp * eps;
        FF      = Fp + Geps + 2.0 * Heps2;
        real Vr = 12 * Cr * VV;
        real Fr = 12 * Cr * FF;
        *V      = Vd + Vr;
        *F      = -(Fd + Fr) * tabscale;
    }
}

real do_walls(const t_inputrec&                   ir,
              const t_forcerec&                   fr,
              const matrix                        box,
              gmx::ArrayRef<const int>            typeA,
              gmx::ArrayRef<const int>            typeB,
              gmx::ArrayRef<const unsigned short> cENER,
              const int                           homenr,
              const int                           numPerturbedAtoms,
              gmx::ArrayRef<const gmx::RVec>      x,
              gmx::ForceWithVirial*               forceWithVirial,
              real                                lambda,
              gmx::ArrayRef<real>                 Vlj,
              t_nrnb*                             nrnb)
{
    constexpr real sixth   = 1.0 / 6.0;
    constexpr real twelfth = 1.0 / 12.0;

    int  ntw[2];
    real fac_d[2], fac_r[2];

    const int   nwall     = ir.nwall;
    const int   ngid      = ir.opts.ngener;
    const int   ntype     = fr.ntype;
    const real* nbfp      = fr.nbfp.data();
    const int*  egp_flags = fr.egp_flags;

    for (int w = 0; w < nwall; w++)
    {
        ntw[w] = 2 * ntype * ir.wall_atomtype[w];
        switch (ir.wall_type)
        {
            case WallType::NineThree:
                fac_d[w] = ir.wall_density[w] * M_PI / 6;
                fac_r[w] = ir.wall_density[w] * M_PI / 45;
                break;
            case WallType::TenFour:
                fac_d[w] = ir.wall_density[w] * M_PI / 2;
                fac_r[w] = ir.wall_density[w] * M_PI / 5;
                break;
            default: break;
        }
    }
    const real wall_z[2] = { 0, box[ZZ][ZZ] };

    rvec* gmx_restrict f = as_rvec_array(forceWithVirial->force_.data());

    real   dvdlambda = 0;
    double sumRF     = 0;
    for (int lam = 0; lam < (numPerturbedAtoms ? 2 : 1); lam++)
    {
        real                     lamfac;
        gmx::ArrayRef<const int> type;
        if (numPerturbedAtoms != 0)
        {
            if (lam == 0)
            {
                lamfac = 1 - lambda;
                type   = typeA;
            }
            else
            {
                lamfac = lambda;
                type   = typeB;
            }
        }
        else
        {
            lamfac = 1;
            type   = typeA;
        }

        real Vlambda = 0;
        for (int i = 0; i < homenr; i++)
        {
            for (int w = 0; w < std::min(nwall, 2); w++)
            {
                /* The wall energy groups are always at the end of the list */
                const int ggid = cENER[i] * ngid + ngid - nwall + w;
                const int at   = type[i];
                /* nbfp now includes the 6/12 derivative prefactors */
                const real Cd = nbfp[ntw[w] + 2 * at] * sixth;
                const real Cr = nbfp[ntw[w] + 2 * at + 1] * twelfth;
                if (!((Cd == 0 && Cr == 0) || (egp_flags[ggid] & EGP_EXCL)))
                {
                    real r, mr;
                    if (w == 0)
                    {
                        r = x[i][ZZ];
                    }
                    else
                    {
                        r = wall_z[1] - x[i][ZZ];
                    }
                    if (r < ir.wall_r_linpot)
                    {
                        mr = ir.wall_r_linpot - r;
                        r  = ir.wall_r_linpot;
                    }
                    else
                    {
                        mr = 0;
                    }
                    if (r <= 0)
                    {
                        wall_error(i, x, r);
                    }

                    real V, F;
                    real r1, r2, r4, Vd, Vr;
                    switch (ir.wall_type)
                    {
                        case WallType::Table:
                            tableForce(r, *fr.wall_tab[w][cENER[i]], Cd, Cr, &V, &F);
                            F *= lamfac;
                            break;
                        case WallType::NineThree:
                            r1 = 1 / r;
                            r2 = r1 * r1;
                            r4 = r2 * r2;
                            Vd = fac_d[w] * Cd * r2 * r1;
                            Vr = fac_r[w] * Cr * r4 * r4 * r1;
                            V  = Vr - Vd;
                            F  = lamfac * (9 * Vr - 3 * Vd) * r1;
                            break;
                        case WallType::TenFour:
                            r1 = 1 / r;
                            r2 = r1 * r1;
                            r4 = r2 * r2;
                            Vd = fac_d[w] * Cd * r4;
                            Vr = fac_r[w] * Cr * r4 * r4 * r2;
                            V  = Vr - Vd;
                            F  = lamfac * (10 * Vr - 4 * Vd) * r1;
                            break;
                        case WallType::TwelveSix:
                            r1 = 1 / r;
                            r2 = r1 * r1;
                            r4 = r2 * r2;
                            Vd = Cd * r4 * r2;
                            Vr = Cr * r4 * r4 * r4;
                            V  = Vr - Vd;
                            F  = lamfac * (12 * Vr - 6 * Vd) * r1;
                            break;
                        default:
                            V = 0;
                            F = 0;
                            break;
                    }
                    if (mr > 0)
                    {
                        V += mr * F;
                    }
                    sumRF += r * F;
                    if (w == 1)
                    {
                        F = -F;
                    }
                    Vlj[ggid] += lamfac * V;
                    Vlambda += V;
                    f[i][ZZ] += F;
                }
            }
        }
        if (numPerturbedAtoms != 0)
        {
            dvdlambda += (lam == 0 ? -1 : 1) * Vlambda;
        }

        inc_nrnb(nrnb, eNR_WALLS, homenr);
    }

    if (forceWithVirial->computeVirial_)
    {
        rvec virial = { 0, 0, static_cast<real>(-0.5 * sumRF) };
        forceWithVirial->addVirialContribution(virial);
    }

    return dvdlambda;
}
