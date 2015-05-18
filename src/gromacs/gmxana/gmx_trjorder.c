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
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    atom_id i;
    real    d2;
} t_order;

t_order *order;

static int ocomp(const void *a, const void *b)
{
    t_order *oa, *ob;

    oa = (t_order *)a;
    ob = (t_order *)b;

    if (oa->d2 < ob->d2)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

int gmx_trjorder(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] orders molecules according to the smallest distance",
        "to atoms in a reference group",
        "or on z-coordinate (with option [TT]-z[tt]).",
        "With distance ordering, it will ask for a group of reference",
        "atoms and a group of molecules. For each frame of the trajectory",
        "the selected molecules will be reordered according to the shortest",
        "distance between atom number [TT]-da[tt] in the molecule and all the",
        "atoms in the reference group. The center of mass of the molecules can",
        "be used instead of a reference atom by setting [TT]-da[tt] to 0.",
        "All atoms in the trajectory are written",
        "to the output trajectory.[PAR]",
        "[THISMODULE] can be useful for e.g. analyzing the n waters closest to a",
        "protein.",
        "In that case the reference group would be the protein and the group",
        "of molecules would consist of all the water atoms. When an index group",
        "of the first n waters is made, the ordered trajectory can be used",
        "with any GROMACS program to analyze the n closest waters.",
        "[PAR]",
        "If the output file is a [REF].pdb[ref] file, the distance to the reference target",
        "will be stored in the B-factor field in order to color with e.g. Rasmol.",
        "[PAR]",
        "With option [TT]-nshell[tt] the number of molecules within a shell",
        "of radius [TT]-r[tt] around the reference group are printed."
    };
    static int      na   = 3, ref_a = 1;
    static real     rcut = 0;
    static gmx_bool bCOM = FALSE, bZ = FALSE;
    t_pargs         pa[] = {
        { "-na", FALSE, etINT,  {&na},
          "Number of atoms in a molecule" },
        { "-da", FALSE, etINT,  {&ref_a},
          "Atom used for the distance calculation, 0 is COM" },
        { "-com", FALSE, etBOOL, {&bCOM},
          "Use the distance to the center of mass of the reference group" },
        { "-r",  FALSE, etREAL, {&rcut},
          "Cutoff used for the distance calculation when computing the number of molecules in a shell around e.g. a protein" },
        { "-z", FALSE, etBOOL, {&bZ},
          "Order molecules on z-coordinate" }
    };
    FILE           *fp;
    t_trxstatus    *out;
    t_trxstatus    *status;
    gmx_bool        bNShell, bPDBout;
    t_topology      top;
    int             ePBC;
    rvec           *x, *xsol, xcom, dx;
    matrix          box;
    t_pbc           pbc;
    gmx_rmpbc_t     gpbc;
    real            t, totmass, mass, rcut2 = 0, n2;
    int             natoms, nwat, ncut;
    char          **grpname, title[256];
    int             i, j, d, *isize, isize_ref = 0, isize_sol;
    atom_id         sa, sr, *swi, **index, *ind_ref = NULL, *ind_sol;
    output_env_t    oenv;
    t_filenm        fnm[] = {
        { efTRX, "-f", NULL, ffREAD  },
        { efTPS, NULL, NULL, ffREAD  },
        { efNDX, NULL, NULL, ffOPTRD },
        { efTRO, "-o", "ordered", ffOPTWR },
        { efXVG, "-nshell", "nshell", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &x, NULL, box, TRUE);
    sfree(x);

    /* get index groups */
    printf("Select %sa group of molecules to be ordered:\n",
           bZ ? "" : "a group of reference atoms and ");
    snew(grpname, 2);
    snew(index, 2);
    snew(isize, 2);
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), bZ ? 1 : 2,
              isize, index, grpname);

    if (!bZ)
    {
        isize_ref = isize[0];
        isize_sol = isize[1];
        ind_ref   = index[0];
        ind_sol   = index[1];
    }
    else
    {
        isize_sol = isize[0];
        ind_sol   = index[0];
    }

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    if (natoms > top.atoms.nr)
    {
        gmx_fatal(FARGS, "Number of atoms in the run input file is larger than in the trjactory");
    }
    for (i = 0; (i < 2); i++)
    {
        for (j = 0; (j < isize[i]); j++)
        {
            if (index[i][j] > natoms)
            {
                gmx_fatal(FARGS, "An atom number in group %s is larger than the number of atoms in the trajectory");
            }
        }
    }

    if ((isize_sol % na) != 0)
    {
        gmx_fatal(FARGS, "Number of atoms in the molecule group (%d) is not a multiple of na (%d)",
                  isize[1], na);
    }

    nwat = isize_sol/na;
    if (ref_a > na)
    {
        gmx_fatal(FARGS, "The reference atom can not be larger than the number of atoms in a molecule");
    }
    ref_a--;
    snew(xsol, nwat);
    snew(order, nwat);
    snew(swi, natoms);
    for (i = 0; (i < natoms); i++)
    {
        swi[i] = i;
    }

    out     = NULL;
    fp      = NULL;
    bNShell = ((opt2bSet("-nshell", NFILE, fnm)) ||
               (opt2parg_bSet("-r", asize(pa), pa)));
    bPDBout = FALSE;
    if (bNShell)
    {
        rcut2   = rcut*rcut;
        fp      = xvgropen(opt2fn("-nshell", NFILE, fnm), "Number of molecules",
                           "Time (ps)", "N", oenv);
        printf("Will compute the number of molecules within a radius of %g\n",
               rcut);
    }
    if (!bNShell || opt2bSet("-o", NFILE, fnm))
    {
        bPDBout = (fn2ftp(opt2fn("-o", NFILE, fnm)) == efPDB);
        if (bPDBout && !top.atoms.pdbinfo)
        {
            fprintf(stderr, "Creating pdbfino records\n");
            snew(top.atoms.pdbinfo, top.atoms.nr);
        }
        out = open_trx(opt2fn("-o", NFILE, fnm), "w");
    }
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms);
    do
    {
        gmx_rmpbc(gpbc, natoms, box, x);
        set_pbc(&pbc, ePBC, box);

        if (ref_a == -1)
        {
            /* Calculate the COM of all solvent molecules */
            for (i = 0; i < nwat; i++)
            {
                totmass = 0;
                clear_rvec(xsol[i]);
                for (j = 0; j < na; j++)
                {
                    sa       = ind_sol[i*na+j];
                    mass     = top.atoms.atom[sa].m;
                    totmass += mass;
                    for (d = 0; d < DIM; d++)
                    {
                        xsol[i][d] += mass*x[sa][d];
                    }
                }
                svmul(1/totmass, xsol[i], xsol[i]);
            }
        }
        else
        {
            /* Copy the reference atom of all solvent molecules */
            for (i = 0; i < nwat; i++)
            {
                copy_rvec(x[ind_sol[i*na+ref_a]], xsol[i]);
            }
        }

        if (bZ)
        {
            for (i = 0; (i < nwat); i++)
            {
                sa           = ind_sol[na*i];
                order[i].i   = sa;
                order[i].d2  = xsol[i][ZZ];
            }
        }
        else if (bCOM)
        {
            totmass = 0;
            clear_rvec(xcom);
            for (i = 0; i < isize_ref; i++)
            {
                mass     = top.atoms.atom[ind_ref[i]].m;
                totmass += mass;
                for (j = 0; j < DIM; j++)
                {
                    xcom[j] += mass*x[ind_ref[i]][j];
                }
            }
            svmul(1/totmass, xcom, xcom);
            for (i = 0; (i < nwat); i++)
            {
                sa = ind_sol[na*i];
                pbc_dx(&pbc, xcom, xsol[i], dx);
                order[i].i   = sa;
                order[i].d2  = norm2(dx);
            }
        }
        else
        {
            /* Set distance to first atom */
            for (i = 0; (i < nwat); i++)
            {
                sa = ind_sol[na*i];
                pbc_dx(&pbc, x[ind_ref[0]], xsol[i], dx);
                order[i].i   = sa;
                order[i].d2  = norm2(dx);
            }
            for (j = 1; (j < isize_ref); j++)
            {
                sr = ind_ref[j];
                for (i = 0; (i < nwat); i++)
                {
                    sa = ind_sol[na*i];
                    pbc_dx(&pbc, x[sr], xsol[i], dx);
                    n2 = norm2(dx);
                    if (n2 < order[i].d2)
                    {
                        order[i].d2  = n2;
                    }
                }
            }
        }

        if (bNShell)
        {
            ncut = 0;
            for (i = 0; (i < nwat); i++)
            {
                if (order[i].d2 <= rcut2)
                {
                    ncut++;
                }
            }
            fprintf(fp, "%10.3f  %8d\n", t, ncut);
        }
        if (out)
        {
            qsort(order, nwat, sizeof(*order), ocomp);
            for (i = 0; (i < nwat); i++)
            {
                for (j = 0; (j < na); j++)
                {
                    swi[ind_sol[na*i]+j] = order[i].i+j;
                }
            }

            /* Store the distance as the B-factor */
            if (bPDBout)
            {
                for (i = 0; (i < nwat); i++)
                {
                    for (j = 0; (j < na); j++)
                    {
                        top.atoms.pdbinfo[order[i].i+j].bfac = sqrt(order[i].d2);
                    }
                }
            }
            write_trx(out, natoms, swi, &top.atoms, 0, t, box, x, NULL, NULL);
        }
    }
    while (read_next_x(oenv, status, &t, x, box));
    close_trj(status);
    if (out)
    {
        close_trx(out);
    }
    if (fp)
    {
        xvgrclose(fp);
    }
    gmx_rmpbc_done(gpbc);

    return 0;
}
