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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static real find_pdb_bfac(t_atoms *atoms, t_resinfo *ri, char *atomnm)
{
    char rresnm[8];
    int  i;

    strcpy(rresnm, *ri->name);
    rresnm[3] = '\0';
    for (i = 0; (i < atoms->nr); i++)
    {
        if ((ri->nr == atoms->resinfo[atoms->atom[i].resind].nr) &&
            (ri->ic == atoms->resinfo[atoms->atom[i].resind].ic) &&
            (strcmp(*atoms->resinfo[atoms->atom[i].resind].name, rresnm) == 0) &&
            (strstr(*atoms->atomname[i], atomnm) != NULL))
        {
            break;
        }
    }
    if (i == atoms->nr)
    {
        fprintf(stderr, "\rCan not find %s%d-%s in pdbfile\n",
                rresnm, ri->nr, atomnm);
        return 0.0;
    }

    return atoms->pdbinfo[i].bfac;
}

void correlate_aniso(const char *fn, t_atoms *ref, t_atoms *calc,
                     const output_env_t oenv)
{
    FILE *fp;
    int   i, j;

    fp = xvgropen(fn, "Correlation between X-Ray and Computed Uij", "X-Ray",
                  "Computed", oenv);
    for (i = 0; (i < ref->nr); i++)
    {
        if (ref->pdbinfo[i].bAnisotropic)
        {
            for (j = U11; (j <= U23); j++)
            {
                fprintf(fp, "%10d  %10d\n", ref->pdbinfo[i].uij[j], calc->pdbinfo[i].uij[j]);
            }
        }
    }
    xvgrclose(fp);
}

static void average_residues(double f[], double **U, int uind,
                             int isize, atom_id index[], real w_rls[],
                             t_atoms *atoms)
{
    int    i, j, start;
    double av, m;

    start = 0;
    av    = 0;
    m     = 0;
    for (i = 0; i < isize; i++)
    {
        av += w_rls[index[i]]*(f != NULL ? f[i] : U[i][uind]);
        m  += w_rls[index[i]];
        if (i+1 == isize ||
            atoms->atom[index[i]].resind != atoms->atom[index[i+1]].resind)
        {
            av /= m;
            if (f != NULL)
            {
                for (j = start; j <= i; j++)
                {
                    f[i] = av;
                }
            }
            else
            {
                for (j = start; j <= i; j++)
                {
                    U[j][uind] = av;
                }
            }
            start = i+1;
            av    = 0;
            m     = 0;
        }
    }
}

void print_dir(FILE *fp, real *Uaver)
{
    real eigvec[DIM*DIM];
    real tmp[DIM*DIM];
    rvec eigval;
    int  d, m;

    fprintf(fp, "MSF     X         Y         Z\n");
    for (d = 0; d < DIM; d++)
    {
        fprintf(fp, " %c ", 'X'+d-XX);
        for (m = 0; m < DIM; m++)
        {
            fprintf(fp, " %9.2e", Uaver[3*m+d]);
        }
        fprintf(fp, "%s\n", m == DIM ? " (nm^2)" : "");
    }

    for (m = 0; m < DIM*DIM; m++)
    {
        tmp[m] = Uaver[m];
    }


    eigensolver(tmp, DIM, 0, DIM, eigval, eigvec);

    fprintf(fp, "\n             Eigenvectors\n\n");
    fprintf(fp, "Eigv  %-8.2e %-8.2e %-8.2e (nm^2)\n\n",
            eigval[2], eigval[1], eigval[0]);
    for (d = 0; d < DIM; d++)
    {
        fprintf(fp, "  %c   ", 'X'+d-XX);
        for (m = DIM-1; m >= 0; m--)
        {
            fprintf(fp, "%7.4f  ", eigvec[3*m+d]);
        }
        fprintf(fp, "\n");
    }
}

int gmx_rmsf(int argc, char *argv[])
{
    const char      *desc[] = {
        "[THISMODULE] computes the root mean square fluctuation (RMSF, i.e. standard ",
        "deviation) of atomic positions in the trajectory (supplied with [TT]-f[tt])",
        "after (optionally) fitting to a reference frame (supplied with [TT]-s[tt]).[PAR]",
        "With option [TT]-oq[tt] the RMSF values are converted to B-factor",
        "values, which are written to a [REF].pdb[ref] file with the coordinates, of the",
        "structure file, or of a [REF].pdb[ref] file when [TT]-q[tt] is specified.",
        "Option [TT]-ox[tt] writes the B-factors to a file with the average",
        "coordinates.[PAR]",
        "With the option [TT]-od[tt] the root mean square deviation with",
        "respect to the reference structure is calculated.[PAR]",
        "With the option [TT]-aniso[tt], [THISMODULE] will compute anisotropic",
        "temperature factors and then it will also output average coordinates",
        "and a [REF].pdb[ref] file with ANISOU records (corresonding to the [TT]-oq[tt]",
        "or [TT]-ox[tt] option). Please note that the U values",
        "are orientation-dependent, so before comparison with experimental data",
        "you should verify that you fit to the experimental coordinates.[PAR]",
        "When a [REF].pdb[ref] input file is passed to the program and the [TT]-aniso[tt]",
        "flag is set",
        "a correlation plot of the Uij will be created, if any anisotropic",
        "temperature factors are present in the [REF].pdb[ref] file.[PAR]",
        "With option [TT]-dir[tt] the average MSF (3x3) matrix is diagonalized.",
        "This shows the directions in which the atoms fluctuate the most and",
        "the least."
    };
    static gmx_bool  bRes    = FALSE, bAniso = FALSE, bdevX = FALSE, bFit = TRUE;
    t_pargs          pargs[] = {
        { "-res", FALSE, etBOOL, {&bRes},
          "Calculate averages for each residue" },
        { "-aniso", FALSE, etBOOL, {&bAniso},
          "Compute anisotropic termperature factors" },
        { "-fit", FALSE, etBOOL, {&bFit},
          "Do a least squares superposition before computing RMSF. Without this you must make sure that the reference structure and the trajectory match." }
    };
    int              natom;
    int              step, nre, natoms, i, g, m, teller = 0;
    real             t, lambda, *w_rls, *w_rms;

    t_tpxheader      header;
    t_inputrec       ir;
    t_topology       top;
    int              ePBC;
    t_atoms         *pdbatoms, *refatoms;
    gmx_bool         bCont;

    matrix           box, pdbbox;
    rvec            *x, *pdbx, *xref;
    t_trxstatus     *status;
    int              npdbatoms, res0;
    char             buf[256];
    const char      *label;
    char             title[STRLEN];

    FILE            *fp;          /* the graphics file */
    const char      *devfn, *dirfn;
    int              resind;

    gmx_bool         bReadPDB;
    atom_id         *index;
    int              isize;
    char            *grpnames;

    real             bfac, pdb_bfac, *Uaver;
    double         **U, *xav;
    atom_id          aid;
    rvec            *rmsd_x = NULL;
    double          *rmsf, invcount, totmass;
    int              d;
    real             count = 0;
    rvec             xcm;
    gmx_rmpbc_t      gpbc = NULL;

    output_env_t     oenv;

    const char      *leg[2] = { "MD", "X-Ray" };

    t_filenm         fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD  },
        { efTPS, NULL,  NULL,     ffREAD  },
        { efNDX, NULL,  NULL,     ffOPTRD },
        { efPDB, "-q",  NULL,     ffOPTRD },
        { efPDB, "-oq", "bfac",   ffOPTWR },
        { efPDB, "-ox", "xaver",  ffOPTWR },
        { efXVG, "-o",  "rmsf",   ffWRITE },
        { efXVG, "-od", "rmsdev", ffOPTWR },
        { efXVG, "-oc", "correl", ffOPTWR },
        { efLOG, "-dir", "rmsf",  ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, asize(pargs), pargs, asize(desc), desc, 0, NULL,
                           &oenv))
    {
        return 0;
    }

    bReadPDB = ftp2bSet(efPDB, NFILE, fnm);
    devfn    = opt2fn_null("-od", NFILE, fnm);
    dirfn    = opt2fn_null("-dir", NFILE, fnm);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xref, NULL, box, TRUE);
    snew(w_rls, top.atoms.nr);

    fprintf(stderr, "Select group(s) for root mean square calculation\n");
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpnames);

    /* Set the weight */
    for (i = 0; i < isize; i++)
    {
        w_rls[index[i]] = top.atoms.atom[index[i]].m;
    }

    /* Malloc the rmsf arrays */
    snew(xav, isize*DIM);
    snew(U, isize);
    for (i = 0; i < isize; i++)
    {
        snew(U[i], DIM*DIM);
    }
    snew(rmsf, isize);
    if (devfn)
    {
        snew(rmsd_x, isize);
    }

    if (bReadPDB)
    {
        get_stx_coordnum(opt2fn("-q", NFILE, fnm), &npdbatoms);
        snew(pdbatoms, 1);
        snew(refatoms, 1);
        init_t_atoms(pdbatoms, npdbatoms, TRUE);
        init_t_atoms(refatoms, npdbatoms, TRUE);
        snew(pdbx, npdbatoms);
        /* Read coordinates twice */
        read_stx_conf(opt2fn("-q", NFILE, fnm), title, pdbatoms, pdbx, NULL, NULL, pdbbox);
        read_stx_conf(opt2fn("-q", NFILE, fnm), title, refatoms, pdbx, NULL, NULL, pdbbox);
    }
    else
    {
        pdbatoms  = &top.atoms;
        refatoms  = &top.atoms;
        pdbx      = xref;
        npdbatoms = pdbatoms->nr;
        snew(pdbatoms->pdbinfo, npdbatoms);
        copy_mat(box, pdbbox);
    }

    if (bFit)
    {
        sub_xcm(xref, isize, index, top.atoms.atom, xcm, FALSE);
    }

    natom = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    if (bFit)
    {
        gpbc = gmx_rmpbc_init(&top.idef, ePBC, natom);
    }

    /* Now read the trj again to compute fluctuations */
    teller = 0;
    do
    {
        if (bFit)
        {
            /* Remove periodic boundary */
            gmx_rmpbc(gpbc, natom, box, x);

            /* Set center of mass to zero */
            sub_xcm(x, isize, index, top.atoms.atom, xcm, FALSE);

            /* Fit to reference structure */
            do_fit(natom, w_rls, xref, x);
        }

        /* Calculate Anisotropic U Tensor */
        for (i = 0; i < isize; i++)
        {
            aid = index[i];
            for (d = 0; d < DIM; d++)
            {
                xav[i*DIM + d] += x[aid][d];
                for (m = 0; m < DIM; m++)
                {
                    U[i][d*DIM + m] += x[aid][d]*x[aid][m];
                }
            }
        }

        if (devfn)
        {
            /* Calculate RMS Deviation */
            for (i = 0; (i < isize); i++)
            {
                aid = index[i];
                for (d = 0; (d < DIM); d++)
                {
                    rmsd_x[i][d] += sqr(x[aid][d]-xref[aid][d]);
                }
            }
        }
        count += 1.0;
        teller++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    close_trj(status);

    if (bFit)
    {
        gmx_rmpbc_done(gpbc);
    }


    invcount = 1.0/count;
    snew(Uaver, DIM*DIM);
    totmass = 0;
    for (i = 0; i < isize; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            xav[i*DIM + d] *= invcount;
        }
        for (d = 0; d < DIM; d++)
        {
            for (m = 0; m < DIM; m++)
            {
                U[i][d*DIM + m] = U[i][d*DIM + m]*invcount
                    - xav[i*DIM + d]*xav[i*DIM + m];
                Uaver[3*d+m] += top.atoms.atom[index[i]].m*U[i][d*DIM + m];
            }
        }
        totmass += top.atoms.atom[index[i]].m;
    }
    for (d = 0; d < DIM*DIM; d++)
    {
        Uaver[d] /= totmass;
    }

    if (bRes)
    {
        for (d = 0; d < DIM*DIM; d++)
        {
            average_residues(NULL, U, d, isize, index, w_rls, &top.atoms);
        }
    }

    if (bAniso)
    {
        for (i = 0; i < isize; i++)
        {
            aid = index[i];
            pdbatoms->pdbinfo[aid].bAnisotropic = TRUE;
            pdbatoms->pdbinfo[aid].uij[U11]     = 1e6*U[i][XX*DIM + XX];
            pdbatoms->pdbinfo[aid].uij[U22]     = 1e6*U[i][YY*DIM + YY];
            pdbatoms->pdbinfo[aid].uij[U33]     = 1e6*U[i][ZZ*DIM + ZZ];
            pdbatoms->pdbinfo[aid].uij[U12]     = 1e6*U[i][XX*DIM + YY];
            pdbatoms->pdbinfo[aid].uij[U13]     = 1e6*U[i][XX*DIM + ZZ];
            pdbatoms->pdbinfo[aid].uij[U23]     = 1e6*U[i][YY*DIM + ZZ];
        }
    }
    if (bRes)
    {
        label = "Residue";
    }
    else
    {
        label = "Atom";
    }

    for (i = 0; i < isize; i++)
    {
        rmsf[i] = U[i][XX*DIM + XX] + U[i][YY*DIM + YY] + U[i][ZZ*DIM + ZZ];
    }

    if (dirfn)
    {
        fprintf(stdout, "\n");
        print_dir(stdout, Uaver);
        fp = gmx_ffopen(dirfn, "w");
        print_dir(fp, Uaver);
        gmx_ffclose(fp);
    }

    for (i = 0; i < isize; i++)
    {
        sfree(U[i]);
    }
    sfree(U);

    /* Write RMSF output */
    if (bReadPDB)
    {
        bfac = 8.0*M_PI*M_PI/3.0*100;
        fp   = xvgropen(ftp2fn(efXVG, NFILE, fnm), "B-Factors",
                        label, "(A\\b\\S\\So\\N\\S2\\N)", oenv);
        xvgr_legend(fp, 2, leg, oenv);
        for (i = 0; (i < isize); i++)
        {
            if (!bRes || i+1 == isize ||
                top.atoms.atom[index[i]].resind != top.atoms.atom[index[i+1]].resind)
            {
                resind    = top.atoms.atom[index[i]].resind;
                pdb_bfac  = find_pdb_bfac(pdbatoms, &top.atoms.resinfo[resind],
                                          *(top.atoms.atomname[index[i]]));

                fprintf(fp, "%5d  %10.5f  %10.5f\n",
                        bRes ? top.atoms.resinfo[top.atoms.atom[index[i]].resind].nr : index[i]+1, rmsf[i]*bfac,
                        pdb_bfac);
            }
        }
        xvgrclose(fp);
    }
    else
    {
        fp = xvgropen(ftp2fn(efXVG, NFILE, fnm), "RMS fluctuation", label, "(nm)", oenv);
        for (i = 0; i < isize; i++)
        {
            if (!bRes || i+1 == isize ||
                top.atoms.atom[index[i]].resind != top.atoms.atom[index[i+1]].resind)
            {
                fprintf(fp, "%5d %8.4f\n",
                        bRes ? top.atoms.resinfo[top.atoms.atom[index[i]].resind].nr : index[i]+1, sqrt(rmsf[i]));
            }
        }
        xvgrclose(fp);
    }

    for (i = 0; i < isize; i++)
    {
        pdbatoms->pdbinfo[index[i]].bfac = 800*M_PI*M_PI/3.0*rmsf[i];
    }

    if (devfn)
    {
        for (i = 0; i < isize; i++)
        {
            rmsf[i] = (rmsd_x[i][XX]+rmsd_x[i][YY]+rmsd_x[i][ZZ])/count;
        }
        if (bRes)
        {
            average_residues(rmsf, NULL, 0, isize, index, w_rls, &top.atoms);
        }
        /* Write RMSD output */
        fp = xvgropen(devfn, "RMS Deviation", label, "(nm)", oenv);
        for (i = 0; i < isize; i++)
        {
            if (!bRes || i+1 == isize ||
                top.atoms.atom[index[i]].resind != top.atoms.atom[index[i+1]].resind)
            {
                fprintf(fp, "%5d %8.4f\n",
                        bRes ? top.atoms.resinfo[top.atoms.atom[index[i]].resind].nr : index[i]+1, sqrt(rmsf[i]));
            }
        }
        xvgrclose(fp);
    }

    if (opt2bSet("-oq", NFILE, fnm))
    {
        /* Write a .pdb file with B-factors and optionally anisou records */
        for (i = 0; i < isize; i++)
        {
            rvec_inc(xref[index[i]], xcm);
        }
        write_sto_conf_indexed(opt2fn("-oq", NFILE, fnm), title, pdbatoms, pdbx,
                               NULL, ePBC, pdbbox, isize, index);
    }
    if (opt2bSet("-ox", NFILE, fnm))
    {
        /* Misuse xref as a temporary array */
        for (i = 0; i < isize; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                xref[index[i]][d] = xcm[d] + xav[i*DIM + d];
            }
        }
        /* Write a .pdb file with B-factors and optionally anisou records */
        write_sto_conf_indexed(opt2fn("-ox", NFILE, fnm), title, pdbatoms, xref, NULL,
                               ePBC, pdbbox, isize, index);
    }
    if (bAniso)
    {
        correlate_aniso(opt2fn("-oc", NFILE, fnm), refatoms, pdbatoms, oenv);
        do_view(oenv, opt2fn("-oc", NFILE, fnm), "-nxy");
    }
    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");
    if (devfn)
    {
        do_view(oenv, opt2fn("-od", NFILE, fnm), "-nxy");
    }

    return 0;
}
