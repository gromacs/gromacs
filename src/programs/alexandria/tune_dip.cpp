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
#include <math.h>
#include "gromacs/fileio/confio.h"
#include "gromacs/utility/futil.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/bonded/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/shellfc.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/random/random.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/legacyheaders/nmsimplex.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "gmx_simple_comm.h"

// Alexandria stuff
#include "gentop_qgen.h"
#include "gentop_core.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "molprop_util.h"
#include "mymol.h"
#include "moldip.h"

static void print_stats(FILE *fp, const char *prop, gmx_stats_t lsq, gmx_bool bHeader,
                        char *xaxis, char *yaxis)
{
    real a, da, b, db, chi2, rmsd, Rfit;
    int  n;

    if (bHeader)
    {
        fprintf(fp, "Fitting data to y = ax+b, where x = %s and y = %s\n",
                xaxis, yaxis);
        fprintf(fp, "%-12s %5s %13s %13s %8s %8s\n",
                "Property", "N", "a", "b", "R", "RMSD");
        fprintf(fp, "---------------------------------------------------------------\n");
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_rmsd(lsq, &rmsd);
    gmx_stats_get_npoints(lsq, &n);
    fprintf(fp, "%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f%% %8.4f\n",
            prop, n, a, da, b, db, Rfit*100, rmsd);
}

static void print_lsq_set(FILE *fp, gmx_stats_t lsq)
{
    real   x, y;

    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, NULL, NULL, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
}

static void print_quad(FILE *fp, tensor Q_exp, tensor Q_calc, char *calc_name,
                       real q_toler)
{
    tensor dQ;
    real   delta;
    if (NULL != calc_name)
    {
        m_sub(Q_exp, Q_calc, dQ);
        delta = sqrt(sqr(dQ[XX][XX])+sqr(dQ[XX][YY])+sqr(dQ[XX][ZZ])+
                     sqr(dQ[YY][YY])+sqr(dQ[YY][ZZ]));
        fprintf(fp,
                "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n",
                calc_name,
                Q_calc[XX][XX], Q_calc[XX][YY], Q_calc[XX][ZZ],
                dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                "", Q_calc[YY][YY], Q_calc[YY][ZZ],
                "", dQ[YY][YY], dQ[YY][ZZ]);
    }
    else
    {
        fprintf(fp, "Quadrupole analysis (5 independent components only)\n");
        fprintf(fp,
                "Exp  (%6.2f %6.2f %6.2f)\n"
                "     (%6s %6.2f %6.2f)\n",
                Q_exp[XX][XX], Q_exp[XX][YY], Q_exp[XX][ZZ],
                "", Q_exp[YY][YY], Q_exp[YY][ZZ]);
    }
}

static void print_dip(FILE *fp, rvec mu_exp, rvec mu_calc, char *calc_name,
                      real toler)
{
    rvec dmu;
    real ndmu, cosa;
    char ebuf[32];

    if (NULL != calc_name)
    {
        rvec_sub(mu_exp, mu_calc, dmu);
        ndmu = norm(dmu);
        cosa = cos_angle(mu_exp, mu_calc);
        if (ndmu > toler)
        {
            sprintf(ebuf, "XXX");
        }
        else if (fabs(cosa) < 0.1)
        {
            sprintf(ebuf, "YYY");
        }
        else
        {
            ebuf[0] = '\0';
        }
        fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                calc_name,
                mu_calc[XX], mu_calc[YY], mu_calc[ZZ], norm(mu_calc),
                dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
    }
    else
    {
        fprintf(fp, "Dipole analysis\n");
        fprintf(fp, "Exp  (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f\n",
                mu_exp[XX], mu_exp[YY], mu_exp[ZZ], norm(mu_exp));
    }
}

static void xvgr_symbolize(FILE *xvgf, int nsym, const char *leg[],
                           const output_env_t oenv)
{
    int i;

    xvgr_legend(xvgf, nsym, leg, oenv);
    for (i = 0; (i < nsym); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", i, i+1);
    }
}

static void print_mols(FILE *fp, const char *xvgfn, const char *qhisto,
                       const char *cdiff, const char *mudiff, const char *Qdiff,
                       const char *espdiff,
                       std::vector<alexandria::MyMol> mol,
                       real hfac,
                       real dip_toler, real quad_toler, real q_toler, output_env_t oenv)
{
    FILE         *xvgf, *qdiff, *mud, *tdiff, *hh, *espd;
    double        d2 = 0;
    real          rms, sigma, aver, error, qq, chi2, espx, espy, espdx, espdy, wtot;
    int           j, k, n, nout, nlsqt = 0, mm, nn;
    char         *resnm, *atomnm;
    const  char **atomtypes = NULL;
    enum {
        eprEEM, eprESP, eprNR
    };
    gmx_stats_t             lsq_q, lsq_mu[eprNR], lsq_quad[eprNR], *lsqt = NULL, lsq_esp;
    const char             *eprnm[eprNR] = { "EEM", "ESP" };
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     at_global, resnr;


    xvgf  = xvgropen(xvgfn, "Correlation between dipoles",
                     "Experimental", "Predicted", oenv);
    xvgr_symbolize(xvgf, 2, eprnm, oenv);
    lsq_q       = gmx_stats_init();
    lsq_quad[0] = gmx_stats_init();
    lsq_quad[1] = gmx_stats_init();
    lsq_mu[0]   = gmx_stats_init();
    lsq_mu[1]   = gmx_stats_init();
    lsq_esp     = gmx_stats_init();
    n           = 0;
    for (std::vector<alexandria::MyMol>::iterator mi = mol.begin(); (mi < mol.end()); mi++)
    {
        if (mi->eSupp != eSupportNo)
        {
            fprintf(fp, "Molecule %d: %s. Qtot: %d, Multiplicity %d\n",
                    n+1,
                    mi->GetMolname().c_str(),
                    mi->GetCharge(),
                    mi->GetMultiplicity());

            print_dip(fp, mi->mu_exp, NULL, NULL, dip_toler);
            print_dip(fp, mi->mu_exp, mi->mu_calc, (char *)"EEM", dip_toler);
            print_dip(fp, mi->mu_exp, mi->mu_esp, (char *)"ESP", dip_toler);

            print_quad(fp, mi->Q_exp, NULL, NULL, quad_toler);
            print_quad(fp, mi->Q_exp, mi->Q_calc, (char *)"EEM", quad_toler);
            print_quad(fp, mi->Q_exp, mi->Q_esp, (char *)"ESP", quad_toler);
            chi2 = gmx_resp_get_rms(mi->gr_, &wtot);
            fprintf(fp, "ESP chi2 %g Hartree/e wtot = %g\n", chi2, wtot);
            gmx_resp_pot_lsq(mi->gr_, lsq_esp);
            while (estatsOK == gmx_stats_get_point(lsq_esp, &espx, &espy,
                                                   &espdx, &espdy, 5))
            {
                fprintf(fp, "ESP outlier: EEM = %g, should be %g\n", espy, espx);
            }

            fprintf(xvgf, "%10g  %10g\n", mi->dip_exp, mi->dip_calc);
            for (mm = 0; (mm < DIM); mm++)
            {
                gmx_stats_add_point(lsq_mu[0], mi->mu_exp[mm], mi->mu_calc[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[1], mi->mu_exp[mm], mi->mu_esp[mm], 0, 0);
                if (0)
                {
                    for (nn = mm; (nn < DIM); nn++)
                    {
                        if (mm < ZZ)
                        {
                            gmx_stats_add_point(lsq_quad[0], mi->Q_exp[mm][nn],
                                                mi->Q_calc[mm][nn], 0, 0);
                            gmx_stats_add_point(lsq_quad[1], mi->Q_exp[mm][nn],
                                                mi->Q_esp[mm][nn], 0, 0);
                        }
                    }
                }
                else
                {
                    /* Ignore off-diagonal components */
                    gmx_stats_add_point(lsq_quad[0], mi->Q_exp[mm][mm],
                                        mi->Q_calc[mm][mm], 0, 0);
                    gmx_stats_add_point(lsq_quad[1], mi->Q_exp[mm][mm],
                                        mi->Q_esp[mm][mm], 0, 0);

                }
            }

            d2 += sqr(mi->dip_exp-mi->dip_calc);
            fprintf(fp, "Atom   Type      q_EEM     q_ESP       x       y       z\n");
            aloop = gmx_mtop_atomloop_all_init(mi->mtop_);
            j     = 0;
            while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
            {
                gmx_mtop_atomloop_all_names(aloop, &atomnm, &resnr, &resnm);
                for (k = 0; (k < nlsqt); k++)
                {
                    if (strcmp(atomtypes[k], *(mi->topology_->atoms.atomtype[j])) == 0)
                    {
                        break;
                    }
                }
                if (k == nlsqt)
                {
                    srenew(lsqt, ++nlsqt);
                    srenew(atomtypes, nlsqt);
                    atomtypes[k] = strdup(*(mi->topology_->atoms.atomtype[j]));
                    lsqt[k]      = gmx_stats_init();
                }
                qq = atom->q;
                fprintf(fp, "%-2s%3d  %-5s  %8.4f  %8.4f%8.3f%8.3f%8.3f %s\n",
                        atomnm, j+1,
                        *(mi->topology_->atoms.atomtype[j]), qq, mi->qESP[j],
                        mi->x_[j][XX], mi->x_[j][YY], mi->x_[j][ZZ],
                        fabs(qq-mi->qESP[j]) > q_toler ? "ZZZ" : "");
                gmx_stats_add_point(lsqt[k], mi->qESP[j], atom->q, 0, 0);
                gmx_stats_add_point(lsq_q, mi->qESP[j], atom->q, 0, 0);
                j++;
            }
            gmx_assert(j, mi->topology_->atoms.nr);
            fprintf(fp, "\n");
            n++;
        }
    }
    fclose(xvgf);

    print_stats(fp, (char *)"dipoles", lsq_mu[0], TRUE, (char *)"Elec", (char *)"EEM");
    print_stats(fp, (char *)"quadrupoles", lsq_quad[0], FALSE, (char *)"Elec", (char *)"EEM");
    print_stats(fp, (char *)"charges", lsq_q, FALSE, (char *)"ESP", (char *)"EEM");
    print_stats(fp, (char *)"esp", lsq_esp, FALSE, (char *)"Elec", (char *)"EEM");
    fprintf(fp, "\n");

    print_stats(fp, (char *)"dipoles", lsq_mu[1], TRUE, (char *)"Elec", (char *)"ESP");
    print_stats(fp, (char *)"quadrupoles", lsq_quad[1], FALSE, (char *)"Elec", (char *)"ESP");

    mud = xvgropen(mudiff, "Correlation between Mu Elec and others",
                   "muElec", "mu", oenv);
    xvgr_symbolize(mud, 2, eprnm, oenv);
    print_lsq_set(mud, lsq_mu[0]);
    print_lsq_set(mud, lsq_mu[1]);
    fclose(mud);

    espd = xvgropen(espdiff, "Correlation between Esp Elec and others",
                    "ESP (Hartree/e)", "ESP (Hartree/e)", oenv);
    xvgr_symbolize(espd, 2, eprnm, oenv);
    print_lsq_set(espd, lsq_esp);
    fclose(espd);

    tdiff = xvgropen(Qdiff, "Correlation between Theta Elec and others",
                     "thetaElec", "theta", oenv);
    xvgr_symbolize(tdiff, 2, eprnm, oenv);
    print_lsq_set(tdiff, lsq_quad[0]);
    print_lsq_set(tdiff, lsq_quad[1]);
    fclose(tdiff);
    qdiff = xvgropen(cdiff, "Correlation between ESP and EEM", "qESP", "qEEM", oenv);
    xvgr_legend(qdiff, nlsqt, atomtypes, oenv);
    xvgr_symbolize(qdiff, nlsqt, atomtypes, oenv);
    hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgr_legend(hh, nlsqt, atomtypes, oenv);
    fprintf(fp, "\nDeviations of the charges separated per atomtype:\n");
    for (k = 0; (k < nlsqt); k++)
    {
        int   N;
        real *x, *y;

        print_stats(fp, atomtypes[k], lsqt[k], (k == 0), (char *)"ESP",
                    (char *)"EEM");
        print_lsq_set(qdiff, lsqt[k]);
        if (gmx_stats_get_npoints(lsqt[k], &N) == estatsOK)
        {
            N = N/4;
            if (gmx_stats_make_histogram(lsqt[k], 0, &N, ehistoY, 0, &x, &y) == estatsOK)
            {
                fprintf(hh, "@type xy\n");
                for (int i = 0; (i < N); i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                sfree(x);
                sfree(y);
            }
        }
    }
    fclose(qdiff);
    fclose(hh);

    rms = sqrt(d2/n);
    fprintf(fp, "RMSD = %.3f D\n", rms);
    fprintf(fp, "hfac = %g\n", hfac);
    gmx_stats_get_ase(lsq_mu[0], &aver, &sigma, &error);
    sigma = rms;
    nout  = 0;
    fprintf(fp, "Overview of outliers (> %.3f off)\n", 2*sigma);
    fprintf(fp, "----------------------------------\n");
    fprintf(fp, "%-20s  %12s  %12s  %12s\n",
            "Name", "Predicted", "Experimental", "Mu-Deviation");
    for (std::vector<alexandria::MyMol>::iterator mi = mol.begin(); (mi < mol.end()); mi++)
    {
        rvec dmu;
        rvec_sub(mi->mu_exp, mi->mu_calc, dmu);
        if ((mi->eSupp != eSupportNo) &&
            (mi->dip_exp > sigma) &&
            (norm(dmu) > 2*sigma))
        {
            fprintf(fp, "%-20s  %12.3f  %12.3f  %12.3f\n",
                    mi->GetMolname().c_str(),
                    mi->dip_calc, mi->dip_exp,
                    mi->dip_calc-mi->dip_exp);
            nout++;
        }
    }
    if (nout)
    {
        printf("There were %d outliers. See at the very bottom of the log file\n",
               nout);
    }
    else
    {
        printf("No outliers! Well done.\n");
    }
    do_view(oenv, xvgfn, NULL);

    gmx_stats_done(lsq_q);
    gmx_stats_done(lsq_quad[0]);
    gmx_stats_done(lsq_quad[1]);
    gmx_stats_done(lsq_mu[0]);
    gmx_stats_done(lsq_mu[1]);
    sfree(lsq_q);
    sfree(lsq_quad[0]);
    sfree(lsq_quad[1]);
    sfree(lsq_mu[0]);
    sfree(lsq_mu[1]);
}

static double dipole_function(void *params, double v[])
{
    alexandria::MolDip *md = (alexandria::MolDip *) params;
    int                 j, k, zz, nzeta;
    double              chi0, z, J0, bounds = 0;
    char               *name;
    char               *qstr, *rowstr;
    char                zstr[STRLEN], buf[STRLEN];

#define HARMONIC(x, xmin, xmax) (x < xmin) ? (sqr(x-xmin)) : ((x > xmax) ? (sqr(x-xmax)) : 0)

    /* Set parameters in eem record. There is a penalty if parameters
     * go out of bounds as well.
     */
    k = 0;
    while ((name = opt_index_count(md->_ic)) != NULL)
    {
        J0      = v[k++];
        bounds += HARMONIC(J0, md->_J0_0, md->_J0_1);
        if (strcasecmp(name, md->_fixchi) != 0)
        {
            chi0    = v[k++];
            bounds += HARMONIC(chi0, md->_Chi0_0, md->_Chi0_1);
        }
        else
        {
            chi0 = gmx_poldata_get_chi0(md->_pd, md->_iModel, name);
        }

        qstr    = gmx_poldata_get_qstr(md->_pd, md->_iModel, name);
        rowstr  = gmx_poldata_get_rowstr(md->_pd, md->_iModel, name);
        nzeta   = gmx_poldata_get_nzeta(md->_pd, md->_iModel, name);
        zstr[0] = '\0';
        for (zz = 0; (zz < nzeta); zz++)
        {
            z = gmx_poldata_get_zeta(md->_pd, md->_iModel, name, zz);
            if ((0 != z) && (md->_bFitZeta))
            {
                z       = v[k++];
                bounds += HARMONIC(z, md->_w_0, md->_w_1);
            }
            sprintf(buf, "  %g", z);
            strcat(zstr, buf);
        }
        gmx_poldata_set_eemprops(md->_pd, md->_iModel, name, J0, chi0,
                                 zstr, qstr, rowstr);
    }
    if (md->_bOptHfac)
    {
        md->_hfac = v[k++];
        if (md->_hfac >  md->_hfac0)
        {
            bounds += 100*sqr(md->_hfac - md->_hfac0);
        }
        else if (md->_hfac < -md->_hfac0)
        {
            bounds += 100*sqr(md->_hfac + md->_hfac0);
        }
    }
    for (j = 0; (j < ermsNR); j++)
    {
        md->_ener[j] = 0;
    }
    md->CalcDeviation();

    /* This contribution is not scaled with force constant because
     * it are essential to charge generation convergence and can hence
     * not be left out of the optimization.
     */
    md->_ener[ermsBOUNDS] += bounds;
    md->_ener[ermsTOT]    += bounds;

    return md->_ener[ermsTOT];
}

static real guess_new_param(real x, real step, real x0, real x1, gmx_rng_t rng,
                            gmx_bool bRandom)
{
    real r = gmx_rng_uniform_real(rng);

    if (bRandom)
    {
        x = x0+(x1-x0)*r;
    }
    else
    {
        x = x*(1-step+2*step*r);
    }

    if (x < x0)
    {
        return x0;
    }
    else if (x > x1)
    {
        return x1;
    }
    else
    {
        return x;
    }
}

static int guess_all_param(FILE *fplog, alexandria::MolDip *md,
                           int run, int iter, real stepsize,
                           gmx_bool bRandom, gmx_rng_t rng,
                           double orig_param[], double test_param[])
{
    double     J00, xxx, chi0, zeta;
    char      *name, *qstr, *rowstr;
    char       zstr[STRLEN], buf[STRLEN];
    gmx_bool   bStart = (/*(run == 0) &&*/ (iter == 0));
    gmx_bool   bRand  = bRandom && (iter == 0);
    int        zz, nzeta, k = 0;

    fprintf(fplog, "%-5s %10s %10s %10s Run/Iter %d/%d - %s randomization\n", "Name",
            "J00", "chi0", "zeta", run, iter,
            bRand ? "Complete" : (bStart ? "Initial" : "Limited"));
    while ((name = opt_index_count(md->_ic)) != NULL)
    {
        if (bStart)
        {
            J00 = gmx_poldata_get_j00(md->_pd, md->_iModel, name);
            xxx = guess_new_param(J00, stepsize, md->_J0_0, md->_J0_1, rng, bRand);
            if (bRand)
            {
                orig_param[k] = xxx;
            }
            else
            {
                orig_param[k] = J00;
            }
            J00 = xxx;
        }
        else
        {
            J00 = guess_new_param(orig_param[k], stepsize, md->_J0_0, md->_J0_1, rng, bRand);
        }
        test_param[k++] = J00;

        chi0 = gmx_poldata_get_chi0(md->_pd, md->_iModel, name);
        if (strcasecmp(name, md->_fixchi) != 0)
        {
            if (bStart)
            {
                xxx = guess_new_param(chi0, stepsize, md->_Chi0_0, md->_Chi0_1, rng, bRand);
                if (bRand)
                {
                    orig_param[k] = xxx;
                }
                else
                {
                    orig_param[k] = chi0;
                }
                chi0 = xxx;
            }
            else
            {
                chi0 = guess_new_param(orig_param[k], stepsize, md->_Chi0_0, md->_Chi0_1, rng, bRand);
            }
            test_param[k++] = chi0;
        }
        if ((qstr = gmx_poldata_get_qstr(md->_pd, md->_iModel, name)) == NULL)
        {
            gmx_fatal(FARGS, "No qstr for atom %s model %d\n", name, md->_iModel);
        }
        if ((rowstr = gmx_poldata_get_rowstr(md->_pd, md->_iModel, name)) == NULL)
        {
            gmx_fatal(FARGS, "No rowstr for atom %s model %d\n", name, md->_iModel);
        }
        nzeta   = gmx_poldata_get_nzeta(md->_pd, md->_iModel, name);
        zstr[0] = '\0';
        for (zz = 0; (zz < nzeta); zz++)
        {
            zeta = gmx_poldata_get_zeta(md->_pd, md->_iModel, name, zz);
            if ((md->_bFitZeta) && (0 != zeta))
            {
                if (bStart)
                {
                    xxx           = guess_new_param(zeta, stepsize, md->_w_0, md->_w_1, rng, bRand);
                    orig_param[k] = (bRand) ? xxx : zeta;
                    zeta          = xxx;
                }
                else
                {
                    zeta = guess_new_param(orig_param[k], stepsize, md->_w_0, md->_w_1, rng, bRand);
                }
                test_param[k++] = zeta;
            }
            sprintf(buf, "  %10g", zeta);
            strcat(zstr, buf);
        }
        gmx_poldata_set_eemprops(md->_pd, md->_iModel, name, J00, chi0,
                                 zstr, qstr, rowstr);
        fprintf(fplog, "%-5s %10g %10g %10s\n", name, J00, chi0, zstr);
    }
    fprintf(fplog, "\n");
    fflush(fplog);
    if (md->_bOptHfac)
    {
        test_param[k++] = md->_hfac;
    }
    return k;
}

static void optimize_moldip(FILE *fp, FILE *fplog, const char *convfn,
                            alexandria::MolDip *md, int maxiter, real tol,
                            int nrun, real stepsize, int seed,
                            gmx_bool bRandom, output_env_t oenv)
{
    FILE      *cfp = NULL;
    double     chi2, chi2_min;
    int        nzeta, zz;
    int        i, k, n, nparam;
    double    *test_param, *orig_param, *best_param, *start;
    gmx_bool   bMinimum = FALSE;
    double     J00, chi0, zeta;
    char      *name, *qstr, *rowstr;
    char       zstr[STRLEN], buf[STRLEN];
    gmx_rng_t  rng;

    if (MASTER(md->_cr))
    {
        rng = gmx_rng_init(seed);

        nparam = 0;
        while ((name = opt_index_count(md->_ic)) != NULL)
        {
            /* One parameter for J00 and one for chi0 */
            nparam++;
            if (strcasecmp(name, md->_fixchi) != 0)
            {
                nparam++;
            }
            if (md->_bFitZeta)
            {
                nzeta  = gmx_poldata_get_nzeta(md->_pd, md->_iModel, name);
                for (i = 0; (i < nzeta); i++)
                {
                    zeta = gmx_poldata_get_zeta(md->_pd, md->_iModel, name, i);
                    if (zeta > 0)
                    {
                        nparam++;
                    }
                }
            }
        }
        /* Check whether we have to optimize the fudge factor for J00 H */
        if (md->_hfac != 0)
        {
            nparam++;
        }
        snew(test_param, nparam+1);
        snew(orig_param, nparam+1);
        snew(best_param, nparam+1);

        /* Starting point */
        snew(start, nparam);

        /* Monitor convergence graphically */
        if (NULL != convfn)
        {
            cfp = xvgropen(convfn, "Convergence", "Value", "Iter", oenv);
        }
        chi2_min = GMX_REAL_MAX;
        for (n = 0; (n < nrun); n++)
        {
            if ((NULL != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n+1, nrun);
                fprintf(fp, "%5s %8s %8s %8s %8s %8s %8s %8s\n",
                        "Run", "d2", "Total", "Bounds",
                        "Dipole", "Quad.", "Charge", "ESP");
            }
            if (NULL != cfp)
            {
                fflush(cfp);
            }

            k = guess_all_param(fplog, md, n, 0, stepsize, bRandom, rng,
                                orig_param, test_param);
            if (k != nparam)
            {
                gmx_fatal(FARGS, "Inconsistency in guess_all_param: k = %d, should be %d", k, nparam);
            }
            for (k = 0; (k < nparam); k++)
            {
                start[k] = test_param[k];
            }

            nmsimplex(cfp, (void *)md, dipole_function, start, nparam,
                      tol, 1, maxiter, &chi2);

            if (chi2 < chi2_min)
            {
                bMinimum = TRUE;
                /* Print convergence if needed */
                if (NULL != cfp)
                {
                    fprintf(cfp, "%5d  ", n*maxiter);
                }
                for (k = 0; (k < nparam); k++)
                {
                    best_param[k] = start[k];
                    if (NULL != cfp)
                    {
                        fprintf(cfp, " %10g", best_param[k]);
                    }
                }
                if (NULL != cfp)
                {
                    fprintf(cfp, "\n");
                }
                chi2_min   = chi2;
            }

            if (fp)
            {
                fprintf(fp, "%5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                        n, chi2,
                        sqrt(md->_ener[ermsTOT]),
                        sqrt(md->_ener[ermsBOUNDS]),
                        sqrt(md->_ener[ermsMU]),
                        sqrt(md->_ener[ermsQUAD]),
                        sqrt(md->_ener[ermsCHARGE]),
                        sqrt(md->_ener[ermsESP]));
            }
        }

        if (bMinimum)
        {
            for (k = 0; (k < nparam); k++)
            {
                start[k] = best_param[k];
            }
            k = 0;
            while ((name = opt_index_count(md->_ic)) != NULL)
            {
                J00    = start[k++];
                chi0   = gmx_poldata_get_chi0(md->_pd, md->_iModel, name);
                if (strcasecmp(name, md->_fixchi) != 0)
                {
                    chi0 = start[k++];
                }
                qstr    = gmx_poldata_get_qstr(md->_pd, md->_iModel, name);
                rowstr  = gmx_poldata_get_rowstr(md->_pd, md->_iModel, name);
                nzeta   = gmx_poldata_get_nzeta(md->_pd, md->_iModel, name);
                zstr[0] = '\0';
                for (zz = 0; (zz < nzeta); zz++)
                {
                    zeta = gmx_poldata_get_zeta(md->_pd, md->_iModel, name, zz);
                    if ((0 != zeta) && md->_bFitZeta)
                    {
                        zeta = start[k++];
                    }
                    sprintf(buf, " %g", zeta);
                    strcat(zstr, buf);
                }
                gmx_poldata_set_eemprops(md->_pd, md->_iModel, name, J00, chi0,
                                         zstr, qstr, rowstr);
            }
            if (md->_bOptHfac)
            {
                md->_hfac = start[k++];
            }
            gmx_assert(k, nparam);

            md->CalcDeviation();
            chi2        = sqrt(md->_ener[ermsTOT]);
            md->_bFinal = TRUE;
            md->CalcDeviation();
            if (fplog)
            {
                fprintf(fplog, "\nMinimum value for RMSD during optimization: %.3f.\n", chi2);
            }
        }
        md->_bDone = TRUE;
        md->CalcDeviation();

        if (NULL != cfp)
        {
            fclose(cfp);
        }
    }
    else
    {
        /* Slave calculators */
        do
        {
            md->CalcDeviation();
        }
        while (!md->_bDone);
    }
}

int alex_tune_dip(int argc, char *argv[])
{
    static const char    *desc[] = {
        "tune_dip read a series of molecules and corresponding experimental",
        "dipole moments from a file, and tunes parameters in an algorithm",
        "until the experimental dipole moments are reproduced by the",
        "charge generating algorithm AX as implemented in the gentop program.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "One of the electronegativities (chi) is redundant in the optimization,",
        "only the relative values are meaningful.",
        "Therefore by default we fix the value for hydrogen to what is written",
        "in the eemprops.dat file (or whatever is given with the [tt]-d[TT] flag).",
        "A suitable value would be 2.3, the original, value due to Pauling,",
        "this can by overridden by setting the [tt]-fixchi[TT] flag to something else (e.g. a non-existing atom).[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm              fnm[] = {
        { efDAT, "-f", "allmols",    ffREAD  },
        { efDAT, "-d", "gentop",     ffOPTRD },
        { efDAT, "-o", "tunedip",    ffWRITE },
        { efDAT, "-sel", "molselect", ffREAD },
        { efLOG, "-g", "charges",    ffWRITE },
        { efXVG, "-x", "dipcorr",    ffWRITE },
        { efXVG, "-qhisto", "q_histo",    ffWRITE },
        { efXVG, "-qdiff", "q_diff",    ffWRITE },
        { efXVG, "-mudiff", "mu_diff",    ffWRITE },
        { efXVG, "-thetadiff", "theta_diff",    ffWRITE },
        { efXVG, "-espdiff", "esp_diff",    ffWRITE },
        { efXVG, "-conv", "convergence", ffOPTWR }
    };
#define NFILE asize(fnm)
    static int            nrun         = 1, maxiter = 100, reinit = 0, seed = 0;
    static int            minimum_data = 3, compress = 1;
    static real           tol          = 1e-3, stol = 1e-6, watoms = 1;
    static gmx_bool       bRandom      = FALSE, bZero = TRUE, bWeighted = TRUE, bOptHfac = FALSE, bQM = FALSE, bCharged = TRUE, bGaussianBug = TRUE, bPol = FALSE, bFitZeta = TRUE;
    static real           J0_0         = 5, Chi0_0 = 1, w_0 = 5, step = 0.01, hfac = 0, rDecrZeta = -1;
    static real           J0_1         = 30, Chi0_1 = 30, w_1 = 50, epsr = 1;
    static real           fc_mu        = 1, fc_bound = 1, fc_quad = 1, fc_charge = 0, fc_esp = 0;
    static real           th_toler     = 170, ph_toler = 5, dip_toler = 0.5, quad_toler = 5, q_toler = 0.25;
    static char          *opt_elem     = NULL, *const_elem = NULL, *fixchi = (char *)"H";
    static char          *lot          = (char *)"B3LYP/aug-cc-pVTZ";
    static char          *qgen[]       = { NULL, (char *)"AXp", (char *)"AXs", (char *)"AXg", NULL };
    t_pargs               pa[]         = {
        { "-tol",   FALSE, etREAL, {&tol},
          "Tolerance for convergence in optimization" },
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations for optimization" },
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-stol",   FALSE, etREAL, {&stol},
          "If reinit is -1 then a reinit will be done as soon as the simplex size is below this treshold." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-charged", FALSE, etBOOL, {&bCharged},
          "Use charged molecules in the parameter tuning as well" },
        { "-qgen",   FALSE, etENUM, {qgen},
          "Algorithm used for charge generation" },
        { "-fixchi", FALSE, etSTR,  {&fixchi},
          "Electronegativity for this element is fixed. Set to FALSE if you want this variable as well, but read the help text above." },
        { "-seed",   FALSE, etINT,  {&seed},
          "Random number seed. If zero, a seed will be generated." },
        { "-j0",    FALSE, etREAL, {&J0_0},
          "Minimum value that J0 (eV) can obtain in fitting" },
        { "-chi0",    FALSE, etREAL, {&Chi0_0},
          "Minimum value that Chi0 (eV) can obtain in fitting" },
        { "-z0",    FALSE, etREAL, {&w_0},
          "Minimum value that inverse radius (1/nm) can obtain in fitting" },
        { "-j1",    FALSE, etREAL, {&J0_1},
          "Maximum value that J0 (eV) can obtain in fitting" },
        { "-chi1",    FALSE, etREAL, {&Chi0_1},
          "Maximum value that Chi0 (eV) can obtain in fitting" },
        { "-z1",    FALSE, etREAL, {&w_1},
          "Maximum value that inverse radius (1/nm) can obtain in fitting" },
        { "-decrzeta", FALSE, etREAL, {&rDecrZeta},
          "Generate decreasing zeta with increasing row numbers for atoms that have multiple distributed charges. In this manner the 1S electrons are closer to the nucleus than 2S electrons and so on. If this number is < 0, nothing is done, otherwise a penalty is imposed in fitting if the Z2-Z1 < this number." },
        { "-fc_bound",    FALSE, etREAL, {&fc_bound},
          "Force constant in the penalty function for going outside the borders given with the above six options." },
        { "-fc_mu",    FALSE, etREAL, {&fc_mu},
          "Force constant in the penalty function for the magnitude of the dipole components." },
        { "-fc_quad",  FALSE, etREAL, {&fc_quad},
          "Force constant in the penalty function for the magnitude of the quadrupole components." },
        { "-fc_esp",   FALSE, etREAL, {&fc_esp},
          "Force constant in the penalty function for the magnitude of the electrostatic potential." },
        { "-fc_charge",  FALSE, etREAL, {&fc_charge},
          "Force constant in the penalty function for the magnitude of the charges with respect to the ESP charges." },
        { "-step",  FALSE, etREAL, {&step},
          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." },
        { "-min_data",  FALSE, etINT, {&minimum_data},
          "Minimum number of data points in order to be able to optimize the parameters for a given atomtype" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of elements to optimize, e.g. \"H C Br\". The other available elements in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of elements to include but keep constant, e.g. \"O N\". These elements from gentop.dat are left unmodified" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use zero. For point+smeared charges 1 is recommended (the default)." },
        { "-pol",    FALSE, etBOOL, {&bPol},
          "Add polarizabilities to the topology and coordinate file" },
        { "-fitzeta", FALSE, etBOOL, {&bFitZeta},
          "Controls whether or not the Gaussian/Slater widths are optimized." },
        { "-zero", FALSE, etBOOL, {&bZero},
          "Use molecules with zero dipole in the fit as well" },
        { "-weight", FALSE, etBOOL, {&bWeighted},
          "Perform a weighted fit, by using the errors in the dipoles presented in the input file. This may or may not improve convergence." },
        { "-hfac",  FALSE, etREAL, {&hfac},
          "Fudge factor to scale the J00 of hydrogen by (1 + hfac * qH). Default hfac is 0, means no fudging." },
        { "-epsr", FALSE, etREAL, {&epsr},
          "Relative dielectric constant to account for intramolecular polarization. Should be >= 1." },
        { "-opthfac",  FALSE, etBOOL, {&bOptHfac},
          "[HIDDEN]Optimize the fudge factor to scale the J00 of hydrogen (see above). If set, then [TT]-hfac[tt] set the absolute value of the largest hfac. Above this, a penalty is incurred." },
        { "-dip_toler", FALSE, etREAL, {&dip_toler},
          "Tolerance (Debye) for marking dipole as an outlier in the log file" },
        { "-quad_toler", FALSE, etREAL, {&quad_toler},
          "Tolerance (Buckingham) for marking quadrupole as an outlier in the log file" },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "Minimum angle to be considered a linear A-B-C bond" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-bgaussquad", FALSE, etBOOL, {&bGaussianBug},
          "[HIDDEN]Work around a bug in the off-diagonal quadrupole components in Gaussian" }
    };
    alexandria::MolDip    md;
    FILE                 *fp;
    ChargeGenerationModel iModel;
    t_commrec            *cr;
    output_env_t          oenv;
    gmx_molselect_t       gms;
    time_t                my_t;
    char                  pukestr[STRLEN];

    cr = init_commrec();

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    if (qgen[0])
    {
        iModel = name2eemtype(qgen[0]);
    }
    else
    {
        iModel = eqgNone;
    }

    if (MASTER(cr))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# %s is part of G R O M A C S:\n#\n", ShortProgram());
        bromacs(pukestr, 99);
        fprintf(fp, "# %s\n#\n", pukestr);
    }
    else
    {
        fp = NULL;
    }
    if (MASTER(cr))
    {
        gms = gmx_molselect_init(opt2fn("-sel", NFILE, fnm));
    }
    else
    {
        gms = NULL;
    }
    md.Init(cr, bQM, bGaussianBug, iModel, rDecrZeta, epsr,
            J0_0, Chi0_0, w_0, J0_1, Chi0_1, w_1,
            fc_bound, fc_mu, fc_quad, fc_charge,
            fc_esp, 1, 1, fixchi, bOptHfac, hfac, bPol, bFitZeta);
    if (0 == seed)
    {
        seed = gmx_rng_make_seed();
    }
    md.Read(fp ? fp : (debug ? debug : NULL),
            opt2fn("-f", NFILE, fnm),
            opt2fn_null("-d", NFILE, fnm),
            minimum_data, bZero,
            opt_elem, const_elem,
            lot, oenv, gms, watoms, TRUE, seed);
    printf("Read %d molecules\n", (int)md._mymol.size());

    optimize_moldip(MASTER(cr) ? stderr : NULL, fp, opt2fn_null("-conv", NFILE, fnm),
                    &md, maxiter, tol, nrun, step, seed,
                    bRandom, oenv);
    if (MASTER(cr))
    {
        print_mols(fp, opt2fn("-x", NFILE, fnm), opt2fn("-qhisto", NFILE, fnm),
                   opt2fn("-qdiff", NFILE, fnm), opt2fn("-mudiff", NFILE, fnm),
                   opt2fn("-thetadiff", NFILE, fnm), opt2fn("-espdiff", NFILE, fnm),
                   md._mymol, md._hfac,
                   dip_toler, quad_toler, q_toler, oenv);

        gmx_ffclose(fp);

        gmx_poldata_write(opt2fn("-o", NFILE, fnm), md._pd, compress);
    }

    return 0;
}
