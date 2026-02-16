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

#include "mk_angndx.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

static int calc_ntype(int nft, const InteractionFunction* ft, const t_idef* idef)
{
    int i, f, nf = 0;

    for (i = 0; (i < idef->ntypes); i++)
    {
        for (f = 0; f < nft; f++)
        {
            if (idef->functype[i] == ft[f])
            {
                nf++;
            }
        }
    }

    return nf;
}

static void fill_ft_ind(int nft, const InteractionFunction* ft, const t_idef* idef, int ft_ind[], char* grpnames[])
{
    char buf[125];
    int  i, f, ind = 0;

    /* Loop over all the function types in the topology */
    for (i = 0; (i < idef->ntypes); i++)
    {
        ft_ind[i] = -1;
        /* Check all the selected function types */
        for (f = 0; f < nft; f++)
        {
            const InteractionFunction ftype = ft[f];
            if (idef->functype[i] == ftype)
            {
                ft_ind[i] = ind;
                switch (ftype)
                {
                    case InteractionFunction::Angles:
                        sprintf(buf,
                                "Theta=%.1f_%.2f",
                                idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case InteractionFunction::GROMOS96Angles:
                        sprintf(buf,
                                "Cos_th=%.1f_%.2f",
                                idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case InteractionFunction::UreyBradleyPotential:
                        sprintf(buf,
                                "UB_th=%.1f_%.2f2f",
                                idef->iparams[i].u_b.thetaA,
                                idef->iparams[i].u_b.kthetaA);
                        break;
                    case InteractionFunction::QuarticAngles:
                        sprintf(buf,
                                "Q_th=%.1f_%.2f_%.2f",
                                idef->iparams[i].qangle.theta,
                                idef->iparams[i].qangle.c[0],
                                idef->iparams[i].qangle.c[1]);
                        break;
                    case InteractionFunction::TabulatedAngles:
                        sprintf(buf,
                                "Table=%d_%.2f",
                                idef->iparams[i].tab.table,
                                idef->iparams[i].tab.kA);
                        break;
                    case InteractionFunction::ProperDihedrals:
                        sprintf(buf,
                                "Phi=%.1f_%d_%.2f",
                                idef->iparams[i].pdihs.phiA,
                                idef->iparams[i].pdihs.mult,
                                idef->iparams[i].pdihs.cpA);
                        break;
                    case InteractionFunction::ImproperDihedrals:
                        sprintf(buf,
                                "Xi=%.1f_%.2f",
                                idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case InteractionFunction::RyckaertBellemansDihedrals:
                        sprintf(buf, "RB-A1=%.2f", idef->iparams[i].rbdihs.rbcA[1]);
                        break;
                    case InteractionFunction::RestrictedBendingPotential:
                        // Fall through intended
                    case InteractionFunction::RestrictedTorsionPotential:
                        sprintf(buf,
                                "Theta=%.1f_%.2f",
                                idef->iparams[i].harmonic.rA,
                                idef->iparams[i].harmonic.krA);
                        break;
                    case InteractionFunction::CombinedBendingTorsionPotential:
                        sprintf(buf, "CBT-A1=%.2f", idef->iparams[i].cbtdihs.cbtcA[1]);
                        break;

                    default:
                        gmx_fatal(FARGS,
                                  "Unsupported function type '%s' selected",
                                  interaction_function[ftype].longname);
                }
                grpnames[ind] = gmx_strdup(buf);
                ind++;
            }
        }
    }
}

static void fill_ang(int                        nft,
                     const InteractionFunction* ft,
                     int                        fac,
                     int                        nr[],
                     int*                       index[],
                     const int                  ft_ind[],
                     const t_topology*          top,
                     gmx_bool                   bNoH,
                     real                       hq)
{
    int           f, i, j, indg, nr_fac;
    gmx_bool      bUse;
    const t_idef* idef;
    t_atom*       atom;
    t_iatom*      ia;


    idef = &top->idef;
    atom = top->atoms.atom;

    for (f = 0; f < nft; f++)
    {
        InteractionFunction ftype = ft[f];
        ia                        = idef->il[ftype].iatoms;
        for (i = 0; (i < idef->il[ftype].nr);)
        {
            indg = ft_ind[ia[0]];
            if (indg == -1)
            {
                gmx_incons("Routine fill_ang");
            }
            bUse = TRUE;
            if (bNoH)
            {
                for (j = 0; j < fac; j++)
                {
                    if (atom[ia[1 + j]].m < 1.5)
                    {
                        bUse = FALSE;
                    }
                }
            }
            if (hq != 0.0F)
            {
                for (j = 0; j < fac; j++)
                {
                    if (atom[ia[1 + j]].m < 1.5 && std::abs(atom[ia[1 + j]].q) < hq)
                    {
                        bUse = FALSE;
                    }
                }
            }
            if (bUse)
            {
                if (nr[indg] % 1000 == 0)
                {
                    srenew(index[indg], fac * (nr[indg] + 1000));
                }
                nr_fac = fac * nr[indg];
                for (j = 0; (j < fac); j++)
                {
                    index[indg][nr_fac + j] = ia[j + 1];
                }
                nr[indg]++;
            }
            ia += interaction_function[ftype].nratoms + 1;
            i += interaction_function[ftype].nratoms + 1;
        }
    }
}

static InteractionFunction* select_ftype(const char* opt, int* nft, int* mult)
{
    InteractionFunction* ft = nullptr;

    if (opt[0] == 'a')
    {
        *mult = 3;
        for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
        {
            if ((interaction_function[ftype].flags & IF_ATYPE) || ftype == InteractionFunction::TabulatedAngles)
            {
                (*nft)++;
                srenew(ft, *nft);
                ft[*nft - 1] = ftype;
            }
        }
    }
    else
    {
        *mult = 4;
        *nft  = 1;
        snew(ft, *nft);
        switch (opt[0])
        {
            case 'd': ft[0] = InteractionFunction::ProperDihedrals; break;
            case 'i': ft[0] = InteractionFunction::ImproperDihedrals; break;
            case 'r': ft[0] = InteractionFunction::RyckaertBellemansDihedrals; break;
            default: break;
        }
    }

    return ft;
}

int gmx_mk_angndx(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] makes an index file for calculation of",
        "angle distributions etc. It uses a run input file ([REF].tpx[ref]) for the",
        "definitions of the angles, dihedrals etc."
    };
    static const char* opt[] = { nullptr, "angle", "dihedral", "improper", "ryckaert-bellemans",
                                 nullptr };
    static gmx_bool    bH    = TRUE;
    static real        hq    = -1;
    t_pargs            pa[]  = { { "-type", FALSE, etENUM, { opt }, "Type of angle" },
                                 { "-hyd", FALSE, etBOOL, { &bH }, "Include angles with atoms with mass < 1.5" },
                                 { "-hq",
                                   FALSE,
                                   etREAL,
                                   { &hq },
                                   "Ignore angles with atoms with mass < 1.5 and magnitude of their charge "
                                               "less than this value" } };

    gmx_output_env_t* oenv;
    FILE*             out;
    t_topology*       top;
    int               i, j, ntype;
    int               nft = 0, mult = 0;
    int**             index;
    int*              ft_ind;
    int*              nr;
    char**            grpnames;
    t_filenm fnm[] = { { efTPR, nullptr, nullptr, ffREAD }, { efNDX, nullptr, "angle", ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    GMX_RELEASE_ASSERT(opt[0] != nullptr, "Options inconsistency; opt[0] is NULL");

    InteractionFunction* ft = select_ftype(opt[0], &nft, &mult);

    top = read_top(ftp2fn(efTPR, NFILE, fnm), nullptr);

    ntype = calc_ntype(nft, ft, &(top->idef));
    snew(grpnames, ntype);
    snew(ft_ind, top->idef.ntypes);
    fill_ft_ind(nft, ft, &top->idef, ft_ind, grpnames);

    snew(nr, ntype);
    snew(index, ntype);
    fill_ang(nft, ft, mult, nr, index, ft_ind, top, !bH, hq);

    out = ftp2FILE(efNDX, NFILE, fnm, "w");
    for (i = 0; (i < ntype); i++)
    {
        if (nr[i] > 0)
        {
            fprintf(out, "[ %s ]\n", grpnames[i]);
            for (j = 0; (j < nr[i] * mult); j++)
            {
                fprintf(out, " %5d", index[i][j] + 1);
                if ((j % 12) == 11)
                {
                    fprintf(out, "\n");
                }
            }
            fprintf(out, "\n");
        }
    }
    gmx_ffclose(out);

    return 0;
}
