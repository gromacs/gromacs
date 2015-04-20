/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <math.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

int gmx_genpr(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] produces an #include file for a topology containing",
        "a list of atom numbers and three force constants for the",
        "[IT]x[it]-, [IT]y[it]-, and [IT]z[it]-direction based on",
        "the contents of the [TT]-f[tt] file. A single isotropic force constant may",
        "be given on the command line instead of three components.[PAR]",
        "WARNING: Position restraints are interactions within molecules, therefore",
        "they must be included within the correct [TT][ moleculetype ][tt]",
        "block in the topology. The atom indices within the",
        "[TT][ position_restraints ][tt] block must be within the range of the",
        "atom indices for that molecule type. Since the atom numbers in every",
        "moleculetype in the topology start at 1 and the numbers in the input file",
        "for [THISMODULE] number consecutively from 1, [THISMODULE] will only",
        "produce a useful file for the first molecule. You may wish to",
        "edit the resulting index file to remove the lines for later atoms,",
        "or construct a suitable index group to provide",
        "as input to [THISMODULE].[PAR]",
        "The [TT]-of[tt] option produces an index file that can be used for",
        "freezing atoms. In this case, the input file must be a [REF].pdb[ref] file.[PAR]",
        "With the [TT]-disre[tt] option, half a matrix of distance restraints",
        "is generated instead of position restraints. With this matrix, that",
        "one typically would apply to C[GRK]alpha[grk] atoms in a protein, one can",
        "maintain the overall conformation of a protein without tieing it to",
        "a specific position (as with position restraints)."
    };
    static rvec        fc           = {1000.0, 1000.0, 1000.0};
    static real        freeze_level = 0.0;
    static real        disre_dist   = 0.1;
    static real        disre_frac   = 0.0;
    static real        disre_up2    = 1.0;
    static gmx_bool    bDisre       = FALSE;
    static gmx_bool    bConstr      = FALSE;
    static real        cutoff       = -1.0;

    t_pargs            pa[] = {
        { "-fc", FALSE, etRVEC, {fc},
          "Force constants (kJ/mol nm^2)" },
        { "-freeze", FALSE, etREAL, {&freeze_level},
          "If the [TT]-of[tt] option or this one is given an index file will be written containing atom numbers of all atoms that have a B-factor less than the level given here" },
        { "-disre", FALSE, etBOOL, {&bDisre},
          "Generate a distance restraint matrix for all the atoms in index" },
        { "-disre_dist", FALSE, etREAL, {&disre_dist},
          "Distance range around the actual distance for generating distance restraints" },
        { "-disre_frac", FALSE, etREAL, {&disre_frac},
          "Fraction of distance to be used as interval rather than a fixed distance. If the fraction of the distance that you specify here is less than the distance given in the previous option, that one is used instead." },
        { "-disre_up2", FALSE, etREAL, {&disre_up2},
          "Distance between upper bound for distance restraints, and the distance at which the force becomes constant (see manual)" },
        { "-cutoff", FALSE, etREAL, {&cutoff},
          "Only generate distance restraints for atoms pairs within cutoff (nm)" },
        { "-constr", FALSE, etBOOL, {&bConstr},
          "Generate a constraint matrix rather than distance restraints. Constraints of type 2 will be generated that do generate exclusions." }
    };
#define npargs asize(pa)

    output_env_t     oenv;
    t_atoms         *atoms = NULL;
    int              i, j, k;
    FILE            *out;
    int              igrp;
    real             d, dd, lo, hi;
    atom_id         *ind_grp;
    const char      *xfn, *nfn;
    char            *gn_grp;
    char             title[STRLEN];
    matrix           box;
    gmx_bool         bFreeze;
    rvec             dx, *x = NULL, *v = NULL;

    t_filenm         fnm[] = {
        { efSTX, "-f",  NULL,    ffREAD },
        { efNDX, "-n",  NULL,    ffOPTRD },
        { efITP, "-o",  "posre", ffWRITE },
        { efNDX, "-of", "freeze",    ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, npargs, pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    bFreeze = opt2bSet("-of", NFILE, fnm) || opt2parg_bSet("-freeze", asize(pa), pa);
    bDisre  = bDisre || opt2parg_bSet("-disre_dist", npargs, pa);
    xfn     = opt2fn_null("-f", NFILE, fnm);
    nfn     = opt2fn_null("-n", NFILE, fnm);

    if (( nfn == NULL ) && ( xfn == NULL))
    {
        gmx_fatal(FARGS, "no index file and no structure file supplied");
    }

    if ((disre_frac < 0) || (disre_frac >= 1))
    {
        gmx_fatal(FARGS, "disre_frac should be between 0 and 1");
    }
    if (disre_dist < 0)
    {
        gmx_fatal(FARGS, "disre_dist should be >= 0");
    }

    if (xfn != NULL)
    {
        snew(atoms, 1);
        get_stx_coordnum(xfn, &(atoms->nr));
        init_t_atoms(atoms, atoms->nr, TRUE);
        snew(x, atoms->nr);
        snew(v, atoms->nr);
        fprintf(stderr, "\nReading structure file\n");
        read_stx_conf(xfn, title, atoms, x, v, NULL, box);
    }

    if (bFreeze)
    {
        if (!atoms || !atoms->pdbinfo)
        {
            gmx_fatal(FARGS, "No B-factors in input file %s, use a pdb file next time.",
                      xfn);
        }

        out = opt2FILE("-of", NFILE, fnm, "w");
        fprintf(out, "[ freeze ]\n");
        for (i = 0; (i < atoms->nr); i++)
        {
            if (atoms->pdbinfo[i].bfac <= freeze_level)
            {
                fprintf(out, "%d\n", i+1);
            }
        }
        gmx_ffclose(out);
    }
    else if ((bDisre || bConstr) && x)
    {
        printf("Select group to generate %s matrix from\n",
               bConstr ? "constraint" : "distance restraint");
        get_index(atoms, nfn, 1, &igrp, &ind_grp, &gn_grp);

        out = ftp2FILE(efITP, NFILE, fnm, "w");
        if (bConstr)
        {
            fprintf(out, "; constraints for %s of %s\n\n", gn_grp, title);
            fprintf(out, "[ constraints ]\n");
            fprintf(out, ";%4s %5s %1s %10s\n", "i", "j", "tp", "dist");
        }
        else
        {
            fprintf(out, "; distance restraints for %s of %s\n\n", gn_grp, title);
            fprintf(out, "[ distance_restraints ]\n");
            fprintf(out, ";%4s %5s %1s %5s %10s %10s %10s %10s %10s\n", "i", "j", "?",
                    "label", "funct", "lo", "up1", "up2", "weight");
        }
        for (i = k = 0; i < igrp; i++)
        {
            for (j = i+1; j < igrp; j++, k++)
            {
                rvec_sub(x[ind_grp[i]], x[ind_grp[j]], dx);
                d = norm(dx);
                if (bConstr)
                {
                    fprintf(out, "%5d %5d %1d %10g\n", ind_grp[i]+1, ind_grp[j]+1, 2, d);
                }
                else
                {
                    if (cutoff < 0 || d < cutoff)
                    {
                        if (disre_frac > 0)
                        {
                            dd = min(disre_dist, disre_frac*d);
                        }
                        else
                        {
                            dd = disre_dist;
                        }
                        lo = max(0, d-dd);
                        hi = d+dd;
                        fprintf(out, "%5d %5d %1d %5d %10d %10g %10g %10g %10g\n",
                                ind_grp[i]+1, ind_grp[j]+1, 1, k, 1,
                                lo, hi, hi+disre_up2, 1.0);
                    }
                }
            }
        }
        gmx_ffclose(out);
    }
    else
    {
        printf("Select group to position restrain\n");
        get_index(atoms, nfn, 1, &igrp, &ind_grp, &gn_grp);

        out = ftp2FILE(efITP, NFILE, fnm, "w");
        fprintf(out, "; position restraints for %s of %s\n\n", gn_grp, title);
        fprintf(out, "[ position_restraints ]\n");
        fprintf(out, ";%3s %5s %9s %10s %10s\n", "i", "funct", "fcx", "fcy", "fcz");
        for (i = 0; i < igrp; i++)
        {
            fprintf(out, "%4d %4d %10g %10g %10g\n",
                    ind_grp[i]+1, 1, fc[XX], fc[YY], fc[ZZ]);
        }
        gmx_ffclose(out);
    }
    if (xfn)
    {
        sfree(x);
        sfree(v);
    }

    return 0;
}
