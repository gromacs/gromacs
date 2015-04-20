/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include <assert.h>
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define SQR(x) (pow(x, 2.0))
#define EPSI0 (EPSILON0*E_CHARGE*E_CHARGE*AVOGADRO/(KILO*NANO)) /* EPSILON0 in SI units */

static void index_atom2mol(int *n, int *index, t_block *mols)
{
    int nat, i, nmol, mol, j;

    nat  = *n;
    i    = 0;
    nmol = 0;
    mol  = 0;
    while (i < nat)
    {
        while (index[i] > mols->index[mol])
        {
            mol++;
            if (mol >= mols->nr)
            {
                gmx_fatal(FARGS, "Atom index out of range: %d", index[i]+1);
            }
        }
        for (j = mols->index[mol]; j < mols->index[mol+1]; j++)
        {
            if (i >= nat || index[i] != j)
            {
                gmx_fatal(FARGS, "The index group does not consist of whole molecules");
            }
            i++;
        }
        index[nmol++] = mol;
    }

    fprintf(stderr, "\nSplit group of %d atoms into %d molecules\n", nat, nmol);

    *n = nmol;
}

static gmx_bool precalc(t_topology top, real mass2[], real qmol[])
{

    real     mtot;
    real     qtot;
    real     qall;
    int      i, j, k, l;
    int      ai, ci;
    gmx_bool bNEU;
    ai   = 0;
    ci   = 0;
    qall = 0.0;



    for (i = 0; i < top.mols.nr; i++)
    {
        k    = top.mols.index[i];
        l    = top.mols.index[i+1];
        mtot = 0.0;
        qtot = 0.0;

        for (j = k; j < l; j++)
        {
            mtot += top.atoms.atom[j].m;
            qtot += top.atoms.atom[j].q;
        }

        for (j = k; j < l; j++)
        {
            top.atoms.atom[j].q -= top.atoms.atom[j].m*qtot/mtot;
            mass2[j]             = top.atoms.atom[j].m/mtot;
            qmol[j]              = qtot;
        }


        qall += qtot;

        if (qtot < 0.0)
        {
            ai++;
        }
        if (qtot > 0.0)
        {
            ci++;
        }
    }

    if (fabs(qall) > 0.01)
    {
        printf("\n\nSystem not neutral (q=%f) will not calculate translational part of the dipole moment.\n", qall);
        bNEU = FALSE;
    }
    else
    {
        bNEU = TRUE;
    }

    return bNEU;

}

static void remove_jump(matrix box, int natoms, rvec xp[], rvec x[])
{

    rvec hbox;
    int  d, i, m;

    for (d = 0; d < DIM; d++)
    {
        hbox[d] = 0.5*box[d][d];
    }
    for (i = 0; i < natoms; i++)
    {
        for (m = DIM-1; m >= 0; m--)
        {
            while (x[i][m]-xp[i][m] <= -hbox[m])
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] += box[m][d];
                }
            }
            while (x[i][m]-xp[i][m] > hbox[m])
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] -= box[m][d];
                }
            }
        }
    }
}

static void calc_mj(t_topology top, int ePBC, matrix box, gmx_bool bNoJump, int isize, int index0[], \
                    rvec fr[], rvec mj, real mass2[], real qmol[])
{

    int   i, j, k, l;
    rvec  tmp;
    rvec  center;
    rvec  mt1, mt2;
    t_pbc pbc;


    calc_box_center(ecenterRECT, box, center);

    if (!bNoJump)
    {
        set_pbc(&pbc, ePBC, box);
    }

    clear_rvec(tmp);


    for (i = 0; i < isize; i++)
    {
        clear_rvec(mt1);
        clear_rvec(mt2);
        k = top.mols.index[index0[i]];
        l = top.mols.index[index0[i+1]];
        for (j = k; j < l; j++)
        {
            svmul(mass2[j], fr[j], tmp);
            rvec_inc(mt1, tmp);
        }

        if (bNoJump)
        {
            svmul(qmol[k], mt1, mt1);
        }
        else
        {
            pbc_dx(&pbc, mt1, center, mt2);
            svmul(qmol[k], mt2, mt1);
        }

        rvec_inc(mj, mt1);

    }

}


static real calceps(real prefactor, real md2, real mj2, real cor, real eps_rf, gmx_bool bCOR)
{

    /* bCOR determines if the correlation is computed via
     * static properties (FALSE) or the correlation integral (TRUE)
     */

    real epsilon = 0.0;


    if (bCOR)
    {
        epsilon = md2-2.0*cor+mj2;
    }
    else
    {
        epsilon = md2+mj2+2.0*cor;
    }

    if (eps_rf == 0.0)
    {
        epsilon = 1.0+prefactor*epsilon;

    }
    else
    {
        epsilon  = 2.0*eps_rf+1.0+2.0*eps_rf*prefactor*epsilon;
        epsilon /= 2.0*eps_rf+1.0-prefactor*epsilon;
    }


    return epsilon;

}


static real calc_cacf(FILE *fcacf, real prefactor, real cacf[], real time[], int nfr, int vfr[], int ei, int nshift)
{

    int  i;
    real corint;
    real deltat = 0.0;
    real rfr;
    real sigma     = 0.0;
    real sigma_ret = 0.0;
    corint = 0.0;

    if (nfr > 1)
    {
        i = 0;

        while (i < nfr)
        {

            rfr      = (real) (nfr/nshift-i/nshift);
            cacf[i] /= rfr;

            if (time[vfr[i]] <= time[vfr[ei]])
            {
                sigma_ret = sigma;
            }

            fprintf(fcacf, "%.3f\t%10.6g\t%10.6g\n", time[vfr[i]], cacf[i], sigma);

            if ((i+1) < nfr)
            {
                deltat = (time[vfr[i+1]]-time[vfr[i]]);
            }
            corint = 2.0*deltat*cacf[i]*prefactor;
            if (i == 0 || (i+1) == nfr)
            {
                corint *= 0.5;
            }
            sigma += corint;

            i++;
        }

    }
    else
    {
        printf("Too less points.\n");
    }

    return sigma_ret;

}

static void calc_mjdsp(FILE *fmjdsp, real prefactor, real dsp2[], real time[], int nfr, real refr[])
{

    int     i;
    real    rtmp;
    real    rfr;


    fprintf(fmjdsp, "#Prefactor fit E-H: 1 / 6.0*V*k_B*T: %g\n", prefactor);



    for (i = 0; i < nfr; i++)
    {

        if (refr[i] != 0.0)
        {
            dsp2[i] *= prefactor/refr[i];
            fprintf(fmjdsp, "%.3f\t%10.6g\n", time[i], dsp2[i]);
        }


    }

}


static void dielectric(FILE *fmj, FILE *fmd, FILE *outf, FILE *fcur, FILE *mcor,
                       FILE *fmjdsp, gmx_bool bNoJump, gmx_bool bACF, gmx_bool bINT,
                       int ePBC, t_topology top, t_trxframe fr, real temp,
                       real trust, real bfit, real efit, real bvit, real evit,
                       t_trxstatus *status, int isize, int nmols, int nshift,
                       atom_id *index0, int indexm[], real mass2[],
                       real qmol[], real eps_rf, const output_env_t oenv)
{
    int       i, j, k, l, f;
    int       valloc, nalloc, nfr, nvfr, m, itrust = 0;
    int       vshfr;
    real     *xshfr       = NULL;
    int      *vfr         = NULL;
    real      refr        = 0.0;
    real      rfr         = 0.0;
    real     *cacf        = NULL;
    real     *time        = NULL;
    real     *djc         = NULL;
    real      corint      = 0.0;
    real      prefactorav = 0.0;
    real      prefactor   = 0.0;
    real      volume;
    real      volume_av = 0.0;
    real      dk_s, dk_d, dk_f;
    real      dm_s, dm_d;
    real      mj    = 0.0;
    real      mj2   = 0.0;
    real      mjd   = 0.0;
    real      mjdav = 0.0;
    real      md2   = 0.0;
    real      mdav2 = 0.0;
    real      sgk;
    rvec      mja_tmp;
    rvec      mjd_tmp;
    rvec      mdvec;
    rvec     *mu    = NULL;
    rvec     *xp    = NULL;
    rvec     *v0    = NULL;
    rvec     *mjdsp = NULL;
    real     *dsp2  = NULL;
    real      t0;
    real      rtmp;
    real      qtmp;



    rvec   tmp;
    rvec  *mtrans = NULL;

    /*
     * Variables for the least-squares fit for Einstein-Helfand and Green-Kubo
     */

    int          bi, ei, ie, ii;
    real         rest  = 0.0;
    real         sigma = 0.0;
    real         malt  = 0.0;
    real         err   = 0.0;
    real        *xfit;
    real        *yfit;
    gmx_rmpbc_t  gpbc = NULL;

    /*
     * indices for EH
     */

    ei = 0;
    bi = 0;

    /*
     * indices for GK
     */

    ii  = 0;
    ie  = 0;
    t0  = 0;
    sgk = 0.0;


    /* This is the main loop over frames */


    nfr = 0;


    nvfr   = 0;
    vshfr  = 0;
    nalloc = 0;
    valloc = 0;

    clear_rvec(mja_tmp);
    clear_rvec(mjd_tmp);
    clear_rvec(mdvec);
    clear_rvec(tmp);
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, fr.natoms);

    do
    {

        refr = (real)(nfr+1);

        if (nfr >= nalloc)
        {
            nalloc += 100;
            srenew(time, nalloc);
            srenew(mu, nalloc);
            srenew(mjdsp, nalloc);
            srenew(dsp2, nalloc);
            srenew(mtrans, nalloc);
            srenew(xshfr, nalloc);

            for (i = nfr; i < nalloc; i++)
            {
                clear_rvec(mjdsp[i]);
                clear_rvec(mu[i]);
                clear_rvec(mtrans[i]);
                dsp2[i]  = 0.0;
                xshfr[i] = 0.0;
            }
        }
        assert(time != NULL);


        if (nfr == 0)
        {
            t0 = fr.time;

        }

        time[nfr] = fr.time-t0;

        if (time[nfr] <= bfit)
        {
            bi = nfr;
        }
        if (time[nfr] <= efit)
        {
            ei = nfr;
        }

        if (bNoJump)
        {

            if (xp)
            {
                remove_jump(fr.box, fr.natoms, xp, fr.x);
            }
            else
            {
                snew(xp, fr.natoms);
            }

            for (i = 0; i < fr.natoms; i++)
            {
                copy_rvec(fr.x[i], xp[i]);
            }

        }

        gmx_rmpbc_trxfr(gpbc, &fr);

        calc_mj(top, ePBC, fr.box, bNoJump, nmols, indexm, fr.x, mtrans[nfr], mass2, qmol);

        for (i = 0; i < isize; i++)
        {
            j = index0[i];
            svmul(top.atoms.atom[j].q, fr.x[j], fr.x[j]);
            rvec_inc(mu[nfr], fr.x[j]);
        }

        /*if(mod(nfr,nshift)==0){*/
        if (nfr%nshift == 0)
        {
            for (j = nfr; j >= 0; j--)
            {
                rvec_sub(mtrans[nfr], mtrans[j], tmp);
                dsp2[nfr-j]  += norm2(tmp);
                xshfr[nfr-j] += 1.0;
            }
        }

        if (fr.bV)
        {
            if (nvfr >= valloc)
            {
                valloc += 100;
                srenew(vfr, valloc);
                if (bINT)
                {
                    srenew(djc, valloc);
                }
                srenew(v0, valloc);
                if (bACF)
                {
                    srenew(cacf, valloc);
                }
            }
            if (time[nfr] <= bvit)
            {
                ii = nvfr;
            }
            if (time[nfr] <= evit)
            {
                ie = nvfr;
            }
            vfr[nvfr] = nfr;
            clear_rvec(v0[nvfr]);
            if (bACF)
            {
                cacf[nvfr] = 0.0;
            }
            if (bINT)
            {
                djc[nvfr] = 0.0;
            }
            for (i = 0; i < isize; i++)
            {
                j = index0[i];
                svmul(mass2[j], fr.v[j], fr.v[j]);
                svmul(qmol[j], fr.v[j], fr.v[j]);
                rvec_inc(v0[nvfr], fr.v[j]);
            }

            fprintf(fcur, "%.3f\t%.6f\t%.6f\t%.6f\n", time[nfr], v0[nfr][XX], v0[nfr][YY], v0[nfr][ZZ]);
            if (bACF || bINT)
            {
                /*if(mod(nvfr,nshift)==0){*/
                if (nvfr%nshift == 0)
                {
                    for (j = nvfr; j >= 0; j--)
                    {
                        if (bACF)
                        {
                            cacf[nvfr-j] += iprod(v0[nvfr], v0[j]);
                        }
                        if (bINT)
                        {
                            djc[nvfr-j] += iprod(mu[vfr[j]], v0[nvfr]);
                        }
                    }
                    vshfr++;
                }
            }
            nvfr++;
        }

        volume     = det(fr.box);
        volume_av += volume;

        rvec_inc(mja_tmp, mtrans[nfr]);
        mjd += iprod(mu[nfr], mtrans[nfr]);
        rvec_inc(mdvec, mu[nfr]);

        mj2 += iprod(mtrans[nfr], mtrans[nfr]);
        md2 += iprod(mu[nfr], mu[nfr]);

        fprintf(fmj, "%.3f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n", time[nfr], mtrans[nfr][XX], mtrans[nfr][YY], mtrans[nfr][ZZ], mj2/refr, norm(mja_tmp)/refr);
        fprintf(fmd, "%.3f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n",    \
                time[nfr], mu[nfr][XX], mu[nfr][YY], mu[nfr][ZZ], md2/refr, norm(mdvec)/refr);

        nfr++;

    }
    while (read_next_frame(oenv, status, &fr));

    gmx_rmpbc_done(gpbc);

    volume_av /= refr;

    prefactor  = 1.0;
    prefactor /= 3.0*EPSILON0*volume_av*BOLTZ*temp;


    prefactorav  = E_CHARGE*E_CHARGE;
    prefactorav /= volume_av*BOLTZMANN*temp*NANO*6.0;

    fprintf(stderr, "Prefactor fit E-H: 1 / 6.0*V*k_B*T: %g\n", prefactorav);

    calc_mjdsp(fmjdsp, prefactorav, dsp2, time, nfr, xshfr);

    /*
     * Now we can average and calculate the correlation functions
     */


    mj2 /= refr;
    mjd /= refr;
    md2 /= refr;

    svmul(1.0/refr, mdvec, mdvec);
    svmul(1.0/refr, mja_tmp, mja_tmp);

    mdav2 = norm2(mdvec);
    mj    = norm2(mja_tmp);
    mjdav = iprod(mdvec, mja_tmp);


    printf("\n\nAverage translational dipole moment M_J [enm] after %d frames (|M|^2): %f %f %f (%f)\n", nfr, mja_tmp[XX], mja_tmp[YY], mja_tmp[ZZ], mj2);
    printf("\n\nAverage molecular dipole moment M_D [enm] after %d frames (|M|^2): %f %f %f (%f)\n", nfr, mdvec[XX], mdvec[YY], mdvec[ZZ], md2);

    if (v0 != NULL)
    {
        trust *= (real) nvfr;
        itrust = (int) trust;

        if (bINT)
        {

            printf("\nCalculating M_D - current correlation integral ... \n");
            corint = calc_cacf(mcor, prefactorav/EPSI0, djc, time, nvfr, vfr, ie, nshift);

        }

        if (bACF)
        {

            printf("\nCalculating current autocorrelation ... \n");
            sgk = calc_cacf(outf, prefactorav/PICO, cacf, time, nvfr, vfr, ie, nshift);

            if (ie > ii)
            {

                snew(xfit, ie-ii+1);
                snew(yfit, ie-ii+1);

                for (i = ii; i <= ie; i++)
                {

                    xfit[i-ii] = log(time[vfr[i]]);
                    rtmp       = fabs(cacf[i]);
                    yfit[i-ii] = log(rtmp);

                }

                lsq_y_ax_b(ie-ii, xfit, yfit, &sigma, &malt, &err, &rest);

                malt = exp(malt);

                sigma += 1.0;

                malt *= prefactorav*2.0e12/sigma;

                sfree(xfit);
                sfree(yfit);

            }
        }
    }


    /* Calculation of the dielectric constant */

    fprintf(stderr, "\n********************************************\n");
    dk_s = calceps(prefactor, md2, mj2, mjd, eps_rf, FALSE);
    fprintf(stderr, "\nAbsolute values:\n epsilon=%f\n", dk_s);
    fprintf(stderr, " <M_D^2> , <M_J^2>, <(M_J*M_D)^2>:  (%f, %f, %f)\n\n", md2, mj2, mjd);
    fprintf(stderr, "********************************************\n");


    dk_f = calceps(prefactor, md2-mdav2, mj2-mj, mjd-mjdav, eps_rf, FALSE);
    fprintf(stderr, "\n\nFluctuations:\n epsilon=%f\n\n", dk_f);
    fprintf(stderr, "\n deltaM_D , deltaM_J, deltaM_JD:  (%f, %f, %f)\n", md2-mdav2, mj2-mj, mjd-mjdav);
    fprintf(stderr, "\n********************************************\n");
    if (bINT)
    {
        dk_d = calceps(prefactor, md2-mdav2, mj2-mj, corint, eps_rf, TRUE);
        fprintf(stderr, "\nStatic dielectric constant using integral and fluctuations: %f\n", dk_d);
        fprintf(stderr, "\n < M_JM_D > via integral:  %.3f\n", -1.0*corint);
    }

    fprintf(stderr, "\n***************************************************");
    fprintf(stderr, "\n\nAverage volume V=%f nm^3 at T=%f K\n", volume_av, temp);
    fprintf(stderr, "and corresponding refactor 1.0 / 3.0*V*k_B*T*EPSILON_0: %f \n", prefactor);



    if (bACF && (ii < nvfr))
    {
        fprintf(stderr, "Integral and integrated fit to the current acf yields at t=%f:\n", time[vfr[ii]]);
        fprintf(stderr, "sigma=%8.3f (pure integral: %.3f)\n", sgk-malt*pow(time[vfr[ii]], sigma), sgk);
    }

    if (ei > bi)
    {
        fprintf(stderr, "\nStart fit at %f ps (%f).\n", time[bi], bfit);
        fprintf(stderr, "End fit at %f ps (%f).\n\n", time[ei], efit);

        snew(xfit, ei-bi+1);
        snew(yfit, ei-bi+1);

        for (i = bi; i <= ei; i++)
        {
            xfit[i-bi] = time[i];
            yfit[i-bi] = dsp2[i];
        }

        lsq_y_ax_b(ei-bi, xfit, yfit, &sigma, &malt, &err, &rest);

        sigma *= 1e12;
        dk_d   = calceps(prefactor, md2, 0.5*malt/prefactorav, corint, eps_rf, TRUE);

        fprintf(stderr, "Einstein-Helfand fit to the MSD of the translational dipole moment yields:\n\n");
        fprintf(stderr, "sigma=%.4f\n", sigma);
        fprintf(stderr, "translational fraction of M^2: %.4f\n", 0.5*malt/prefactorav);
        fprintf(stderr, "Dielectric constant using EH: %.4f\n", dk_d);

        sfree(xfit);
        sfree(yfit);
    }
    else
    {
        fprintf(stderr, "Too few points for a fit.\n");
    }


    if (v0 != NULL)
    {
        sfree(v0);
    }
    if (bACF)
    {
        sfree(cacf);
    }
    if (bINT)
    {
        sfree(djc);
    }

    sfree(time);


    sfree(mjdsp);
    sfree(mu);
}



int gmx_current(int argc, char *argv[])
{

    static int             nshift  = 1000;
    static real            temp    = 300.0;
    static real            trust   = 0.25;
    static real            eps_rf  = 0.0;
    static gmx_bool        bNoJump = TRUE;
    static real            bfit    = 100.0;
    static real            bvit    = 0.5;
    static real            efit    = 400.0;
    static real            evit    = 5.0;
    t_pargs                pa[]    = {
        { "-sh", FALSE, etINT, {&nshift},
          "Shift of the frames for averaging the correlation functions and the mean-square displacement."},
        { "-nojump", FALSE, etBOOL, {&bNoJump},
          "Removes jumps of atoms across the box."},
        { "-eps", FALSE, etREAL, {&eps_rf},
          "Dielectric constant of the surrounding medium. The value zero corresponds to infinity (tin-foil boundary conditions)."},
        { "-bfit", FALSE, etREAL, {&bfit},
          "Begin of the fit of the straight line to the MSD of the translational fraction of the dipole moment."},
        { "-efit", FALSE, etREAL, {&efit},
          "End of the fit of the straight line to the MSD of the translational fraction of the dipole moment."},
        { "-bvit", FALSE, etREAL, {&bvit},
          "Begin of the fit of the current autocorrelation function to a*t^b."},
        { "-evit", FALSE, etREAL, {&evit},
          "End of the fit of the current autocorrelation function to a*t^b."},
        { "-tr", FALSE, etREAL, {&trust},
          "Fraction of the trajectory taken into account for the integral."},
        { "-temp", FALSE, etREAL, {&temp},
          "Temperature for calculating epsilon."}
    };

    output_env_t           oenv;
    t_topology             top;
    char                   title[STRLEN];
    char                 **grpname = NULL;
    const char            *indexfn;
    t_trxframe             fr;
    real                  *mass2 = NULL;
    rvec                  *xtop, *vtop;
    matrix                 box;
    atom_id               *index0;
    int                   *indexm = NULL;
    int                    isize;
    t_trxstatus           *status;
    int                    flags = 0;
    gmx_bool               bTop;
    gmx_bool               bNEU;
    gmx_bool               bACF;
    gmx_bool               bINT;
    int                    ePBC = -1;
    int                    natoms;
    int                    nmols;
    int                    i, j, k = 0, l;
    int                    step;
    real                   t;
    real                   lambda;
    real                  *qmol;
    FILE                  *outf   = NULL;
    FILE                  *outi   = NULL;
    FILE                  *tfile  = NULL;
    FILE                  *mcor   = NULL;
    FILE                  *fmj    = NULL;
    FILE                  *fmd    = NULL;
    FILE                  *fmjdsp = NULL;
    FILE                  *fcur   = NULL;
    t_filenm               fnm[]  = {
        { efTPS,  NULL,  NULL, ffREAD }, /* this is for the topology */
        { efNDX, NULL, NULL, ffOPTRD },
        { efTRX, "-f", NULL, ffREAD },   /* and this for the trajectory */
        { efXVG, "-o",   "current", ffWRITE },
        { efXVG, "-caf", "caf",     ffOPTWR },
        { efXVG, "-dsp", "dsp",     ffWRITE },
        { efXVG, "-md",  "md",      ffWRITE },
        { efXVG, "-mj",  "mj",      ffWRITE },
        { efXVG, "-mc",  "mc",      ffOPTWR }
    };

#define NFILE asize(fnm)


    const char *desc[] = {
        "[THISMODULE] is a tool for calculating the current autocorrelation function, the correlation",
        "of the rotational and translational dipole moment of the system, and the resulting static",
        "dielectric constant. To obtain a reasonable result, the index group has to be neutral.",
        "Furthermore, the routine is capable of extracting the static conductivity from the current ",
        "autocorrelation function, if velocities are given. Additionally, an Einstein-Helfand fit ",
        "can be used to obtain the static conductivity."
        "[PAR]",
        "The flag [TT]-caf[tt] is for the output of the current autocorrelation function and [TT]-mc[tt] writes the",
        "correlation of the rotational and translational part of the dipole moment in the corresponding",
        "file. However, this option is only available for trajectories containing velocities.",
        "Options [TT]-sh[tt] and [TT]-tr[tt] are responsible for the averaging and integration of the",
        "autocorrelation functions. Since averaging proceeds by shifting the starting point",
        "through the trajectory, the shift can be modified with [TT]-sh[tt] to enable the choice of uncorrelated",
        "starting points. Towards the end, statistical inaccuracy grows and integrating the",
        "correlation function only yields reliable values until a certain point, depending on",
        "the number of frames. The option [TT]-tr[tt] controls the region of the integral taken into account",
        "for calculating the static dielectric constant.",
        "[PAR]",
        "Option [TT]-temp[tt] sets the temperature required for the computation of the static dielectric constant.",
        "[PAR]",
        "Option [TT]-eps[tt] controls the dielectric constant of the surrounding medium for simulations using",
        "a Reaction Field or dipole corrections of the Ewald summation ([TT]-eps[tt]\\=0 corresponds to",
        "tin-foil boundary conditions).",
        "[PAR]",
        "[TT]-[no]nojump[tt] unfolds the coordinates to allow free diffusion. This is required to get a continuous",
        "translational dipole moment, required for the Einstein-Helfand fit. The results from the fit allow",
        "the determination of the dielectric constant for system of charged molecules. However, it is also possible to extract",
        "the dielectric constant from the fluctuations of the total dipole moment in folded coordinates. But this",
        "option has to be used with care, since only very short time spans fulfill the approximation that the density",
        "of the molecules is approximately constant and the averages are already converged. To be on the safe side,",
        "the dielectric constant should be calculated with the help of the Einstein-Helfand method for",
        "the translational part of the dielectric constant."
    };


    /* At first the arguments will be parsed and the system information processed */
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    bACF = opt2bSet("-caf", NFILE, fnm);
    bINT = opt2bSet("-mc", NFILE, fnm);

    bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xtop, &vtop, box, TRUE);



    sfree(xtop);
    sfree(vtop);
    indexfn = ftp2fn_null(efNDX, NFILE, fnm);
    snew(grpname, 1);

    get_index(&(top.atoms), indexfn, 1, &isize, &index0, grpname);

    flags = flags | TRX_READ_X | TRX_READ_V;

    read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);

    snew(mass2, top.atoms.nr);
    snew(qmol, top.atoms.nr);

    bNEU = precalc(top, mass2, qmol);


    snew(indexm, isize);

    for (i = 0; i < isize; i++)
    {
        indexm[i] = index0[i];
    }

    nmols = isize;


    index_atom2mol(&nmols, indexm, &top.mols);

    if (fr.bV)
    {
        if (bACF)
        {
            outf = xvgropen(opt2fn("-caf", NFILE, fnm),
                            "Current autocorrelation function", output_env_get_xvgr_tlabel(oenv),
                            "ACF (e nm/ps)\\S2", oenv);
            fprintf(outf, "# time\t acf\t average \t std.dev\n");
        }
        fcur = xvgropen(opt2fn("-o", NFILE, fnm),
                        "Current", output_env_get_xvgr_tlabel(oenv), "J(t) (e nm/ps)", oenv);
        fprintf(fcur, "# time\t Jx\t Jy \t J_z \n");
        if (bINT)
        {
            mcor = xvgropen(opt2fn("-mc", NFILE, fnm),
                            "M\\sD\\N - current  autocorrelation function",
                            output_env_get_xvgr_tlabel(oenv),
                            "< M\\sD\\N (0)\\c7\\CJ(t) >  (e nm/ps)\\S2", oenv);
            fprintf(mcor, "# time\t M_D(0) J(t) acf \t Integral acf\n");
        }
    }

    fmj = xvgropen(opt2fn("-mj", NFILE, fnm),
                   "Averaged translational part of M", output_env_get_xvgr_tlabel(oenv),
                   "< M\\sJ\\N > (enm)", oenv);
    fprintf(fmj, "# time\t x\t y \t z \t average of M_J^2 \t std.dev\n");
    fmd = xvgropen(opt2fn("-md", NFILE, fnm),
                   "Averaged rotational part of M", output_env_get_xvgr_tlabel(oenv),
                   "< M\\sD\\N > (enm)", oenv);
    fprintf(fmd, "# time\t x\t y \t z \t average of M_D^2 \t std.dev\n");

    fmjdsp = xvgropen(opt2fn("-dsp", NFILE, fnm),
                      "MSD of the squared translational dipole moment M",
                      output_env_get_xvgr_tlabel(oenv),
                      "<|M\\sJ\\N(t)-M\\sJ\\N(0)|\\S2\\N > / 6.0*V*k\\sB\\N*T / Sm\\S-1\\Nps\\S-1\\N",
                      oenv);


    /* System information is read and prepared, dielectric() processes the frames
     * and calculates the requested quantities */

    dielectric(fmj, fmd, outf, fcur, mcor, fmjdsp, bNoJump, bACF, bINT, ePBC, top, fr,
               temp, trust, bfit, efit, bvit, evit, status, isize, nmols, nshift,
               index0, indexm, mass2, qmol, eps_rf, oenv);

    xvgrclose(fmj);
    xvgrclose(fmd);
    xvgrclose(fmjdsp);
    if (fr.bV)
    {
        if (bACF)
        {
            xvgrclose(outf);
        }
        xvgrclose(fcur);
        if (bINT)
        {
            xvgrclose(mcor);
        }
    }

    return 0;
}
