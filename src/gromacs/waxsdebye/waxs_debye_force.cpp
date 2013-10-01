/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
#include <fstream>
#include <stdio.h>
#include <string.h>
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/exceptions.h"
#include "waxs_debye_force_c.h"
#include "waxs_debye_force.h"
#include "scattering_factors.h"

/*! \internal \file
 * \brief
 * Implements gmx::WaxsDebyeForce
 *
 * \libinternal
 * \author Alexander Bjorling <alexander.bjorling@chem.gu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

using namespace std;

//! Compute distance using periodic boundary conditions
static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

namespace gmx
{

void WaxsDebyeForce::readData(const char *waxs_ref,
                              const char *waxs_diff,
                              FILE       *fplog)
{
    std::vector<double> ref_error, diff_error;

    /* Read the experimental data */
    for (unsigned int i = 0; (i <= 1); i++)
    {
        const char         *waxs_data;
        bool                bDiff = (bool) i;
        if (bDiff)
        {
            waxs_data = waxs_diff;
        }
        else
        {
            waxs_data = waxs_ref;
        }
        if (NULL == waxs_data)
        {
            continue;
        }
        gmx::File   ifs(waxs_data, "r");
        std::string line;
        int         n = 0;
        while (true == ifs.readLine(&line))
        {
            double      q, S, dS;
            char        buf[STRLEN];

            strncpy(buf, line.c_str(), STRLEN);
            strip_comment(buf);

            /* if experimental errors are provided: */
            if (sscanf(buf, "%15lf%15lf%15lf", &q, &S, &dS) == 3)
            {
                if (!bDiff)
                {
                    q_.push_back(q);
                    Sref_.push_back(S);
                    ref_error.push_back(dS);
                    /* Initiate internal variables */
                    sigma_.push_back(1);
                    A_.push_back(0);
                    Scalc_.push_back(0);
                    dSexp_.push_back(0);
                }
                else
                {
                    // If the q are not exactly equal we can not
                    // do the comparison between experiment and calculation
                    if (q != q_[n])
                    {
                        char errbuf[STRLEN];
                        sprintf(errbuf, "Inconsistency between reference data and difference data. Expected q[%d] to be %f but found %f", n, q_[n], q);
                        GMX_THROW(InvalidInputError(errbuf));
                    }

                    dSexp_[n] = S;
                    diff_error.push_back(dS);
                }
                n++;
            }
            // If experimental errors are not provided:
            else if (sscanf(buf, "%15lf%15lf", &q, &S) == 2)
            {
                if (!bDiff)
                {
                    q_.push_back(q);
                    Sref_.push_back(S);
                    /* Initiate internal variables */
                    sigma_.push_back(1);
                    A_.push_back(0);
                    Scalc_.push_back(0);
                    dSexp_.push_back(0);
                }
                else
                {
                    // See comment above
                    if (q != q_[n])
                    {
                        char errbuf[STRLEN];
                        sprintf(errbuf, "Inconsistency between reference data and difference data. Expected q[%d] to be %f but found %f", n, q_[n], q);
                        GMX_THROW(InvalidInputError(errbuf));
                    }

                    dSexp_[n]  = S;
                }
                n++;
            }
        }
        ifs.close();
        if (NULL != fplog)
        {
            fprintf(fplog, "WAXS: Read %d data points from %s.\n", n, waxs_data);
        }
    }

    /* Fill in sigma_: errors from waxs_diff take precedence */
    if (diff_error.size() == q_.size())
    {
        for (unsigned int i = 0; i < q_.size(); i++)
        {
            sigma_[i] = diff_error[i]/dSexp_[i] + 1.0;
        }
        if (NULL != fplog)
        {
            fprintf(fplog, "WAXS: Using experimental errors from %s.\n", waxs_diff);
        }
    }
    else if (ref_error.size() == q_.size())
    {
        for (unsigned int i = 0; i < q_.size(); i++)
        {
            sigma_[i] = ref_error[i]/Sref_[i] + 1.0;
        }
        if (NULL != fplog)
        {
            fprintf(fplog, "WAXS: Using experimental errors from %s.\n", waxs_ref);
        }
    }
    else if ( (ref_error.empty()) && (diff_error.empty()) )
    {
        if (NULL != fplog)
        {
            fprintf(fplog, "WAXS: No experimental errors read, using zero.\n");
        }
    }
    else
    {
        char errbuf[STRLEN];
        sprintf(errbuf, "Inconsistent experimental WAXS error, expected 0 or %d values", (int) q_.size());

        GMX_THROW(InvalidInputError(errbuf));
    }

}

WaxsDebyeForce::WaxsDebyeForce(FILE             *fplog,
                               const char       *sfactor,
                               const char       *waxs_ref,
                               const char       *waxs_diff,
                               const char       *waxs_out,
                               const char       *waxs_alpha,
                               const t_commrec  *cr,
                               const t_inputrec *ir)
{
    std::vector<double> tmp;
    unsigned int        m, n;
    int                 atomic_number;

    if (NULL != fplog)
    {
        fprintf(fplog, "WAXS: Refinement code initiating.\n");
    }

    /* Initialize step count */
    step_       = 0;

    // Parameters from mdp file.
    rmin2_  = sqr(ir->waxs.debye_r_min);
    rmax2_  = sqr(ir->waxs.debye_r_max);
    kWaxs_  = ir->waxs.kwaxs;

    /* Starting values for the optimization of alpha and other run parameters */
    alpha_mode_ = ir->userint1;
    alpha_min_  = ir->waxs.debye_alpha_min;
    alpha_max_  = ir->waxs.debye_alpha_max;
    nstout_     = ir->waxs.nstout;

    /* Starting values for alpha-dynamics */
    alpha_      = 0.5*(alpha_min_ + alpha_max_);
    alpha_dot_  = 0;
    Falpha_     = 0;
    Malpha_     = 1000;  // 1000
    alpha_step_ = 5e-4;

    /* Open data files */
    readData(waxs_ref, waxs_diff, fplog);

    // Output file for Scalc(q,t)
    if ((nstout_ > 0) && MASTER(cr) && (NULL != waxs_out))
    {
        fp_out_.reset(new File(waxs_out, "w"));
    }

    if ((alpha_step_ != 0.0) && MASTER(cr) && (NULL != waxs_alpha))
    {
        fp_alpha_.reset(new File(waxs_alpha, "w"));
    }

    /* Read in the table for scattering factor data */
    ScatteringFactorTable cmt;
    if (!cmt.read(sfactor))
    {
        GMX_THROW(InvalidInputError(sfactor));
    }
    if (!cmt.write("dump.out"))
    {
        GMX_THROW(InvalidInputError("Can not write to dump.out"));
    }

    /* Fill the F_ array for all types
       allocate space first: atomic_number is one-based so add one to get F_'s size.
     */
    F_.resize(cmt.maxType()+1);
    F0_.resize(cmt.maxType()+1);
    for (m = 0; (m < cmt.size()); m++)
    {
        // Fetch the right atomic number from cmt
        atomic_number = cmt.type(m);

        // Scattering factor at q = 0
        F0_[atomic_number] = cmt.computeScatteringFactor(atomic_number, 0);

        // loop over q-points for this specific scattering factor & store in tmp
        for (n = 0; n < Scalc_.size(); n++)
        {
            tmp.push_back(cmt.computeScatteringFactor(atomic_number, q_[n]));
        }

        // add tmp to the F_ array
        F_[atomic_number] = tmp;
        tmp.clear();
    }

    if (NULL != fplog)
    {
        if (NULL != waxs_diff)
        {
            fprintf(fplog, "WAXS: Difference data file %s read.\n", waxs_diff);
        }
        fprintf(fplog, "WAXS: There were %d (atomic/residue) scattering factors in %s.\n",
                (int)cmt.size(), sfactor);
        if (alpha_mode_ == 1)
        {
            fprintf(fplog, "WAXS: Doing golden ratio minimization of alpha between %5.3f and %5.3f\n", alpha_min_, alpha_max_);
        }
        else if (alpha_mode_ == 2)
        {
            fprintf(fplog, "WAXS: Doing alpha-dynamics between %5.3f and %5.3f\n", alpha_min_, alpha_max_);
        }
    }
}

void WaxsDebyeForce::finalize()
{
    if (fp_out_)
    {
        fp_out_->close();
    }
    if (fp_alpha_)
    {
        fp_alpha_->close();
    }
}

real WaxsDebyeForce::calcS0(FILE            *fplog,
                            const int        nbonds,
                            const t_iatom    iatoms[],
                            const t_iparams  forceparams[],
                            const t_commrec *cr)
{
    double       S0;

    /* Compute S0 */
    S0 = 0;
    for (int i = 0; (i < nbonds); )
    {
        int type = iatoms[i++];
        int ai   = iatoms[i++];
        int aj   = iatoms[i++];
        int gi   = forceparams[type].waxs_debye.tpi;
        int gj   = forceparams[type].waxs_debye.tpj;

        // Should probably check whether gi and gj are less
        if (((0 < gi) && (gi < (int)F0_.size())) &&
            ((0 < gj) && (gj < (int)F0_.size())))
        {
            double dS = F0_[gi]*F0_[gj];
            if (ai == aj)
            {
                S0 += dS;
            }
            else
            {
                // Since we typically only store half the matrix of pairs
                S0 += 2*dS;
            }
        }
    }

    /* Global sum of Scalc_ */
    if (PAR(cr))
    {
        gmx_sumd(1, &S0, cr);
    }
    if (NULL != fplog)
    {
        fprintf(fplog, "WAXS: S(0) = %g.\n", S0);
    }

    return S0;
}


real WaxsDebyeForce::calc(FILE            *fplog,
                          const int        nbonds,
                          const t_iatom    iatoms[],
                          const t_iparams  forceparams[],
                          rvec             x[],
                          rvec             f[],
                          const t_pbc     *pbc,
                          const t_commrec *cr,
                          t_nrnb          *nrnb)
{
    int          i, m, ai, aj, gi, gj, type;
    unsigned int n;
    real         dr, dr2, fscal, fff, vtot, qdr, dr_1, dS;
    rvec         dx;

    if (0 == step_)
    {
        (void) calcS0(fplog, nbonds, iatoms, forceparams, cr);
    }
    /* Reset Scalc */
    for (n = 0; n < Scalc_.size(); n++)
    {
        Scalc_[n] = 0;
    }

    /* Compute Scalc_ for this node only */
    for (i = 0; (i < nbonds); )
    {
        type = iatoms[i++];
        ai   = iatoms[i++];
        aj   = iatoms[i++];
        gi   = forceparams[type].waxs_debye.tpi;
        gj   = forceparams[type].waxs_debye.tpj;

        pbc_rvec_sub(pbc, x[ai], x[aj], dx);

        dr2 = iprod(dx, dx);
        // Special case where the distance is zero
        if (ai == aj)
        {
            for (n = 0; (n < Scalc_.size()); n++)
            {
                double f2 = sqr(F_[gi][n]);
                Scalc_[n] += f2;
            }
        }
        // We only consider particles within a range of distances
        else if ((dr2 > rmin2_) && (dr2 < rmax2_))
        {
            dr = dr2*gmx_invsqrt(dr2);
            for (n = 0; (n < Scalc_.size()); n++)
            {
                qdr        = q_[n]*dr;
                double fgigj = 2.0*F_[gi][n]*F_[gj][n];
                Scalc_[n] += fgigj*sin(qdr)/(qdr);
            }
        }
    }

    /* Global sum of Scalc_ */
    if (PAR(cr))
    {
        gmx_sumd(Scalc_.size(), &(Scalc_[0]), cr);
    }

    /* Update alpha_ based on Scalc_ before calculating forces */
    updateAlpha();

    /* Construction to export Scalc(q,t) to a textfile */
    dump();

    /* Now we calculate A_, the factor
     * A_ = k_chi/(2 SigmaQ^2) * (dSexp - alpha*dScalc)
     * that we will need repeatedly in the force calculation loop.
     * In addition we compute d chi^2/d alpha
     */
    Falpha_     = 0.0;
    vtot        = 0.0;
    for (n = 0; n < Scalc_.size(); n++)
    {
        dS       = dSexp_[n] - alpha_*(Scalc_[n] - Sref_[n]);
        A_[n]    = kWaxs_/(2*sqr(sigma_[n])) * dS;
        vtot    += A_[n]*dS;
        Falpha_ += kWaxs_ * dS/sqr(sigma_[n]) * (Scalc_[n] - Sref_[n]);
    }

    step_++;

    /* loop over the scattering pairs again */
    for (i = 0; (i < nbonds); )
    {
        type = iatoms[i++];
        ai   = iatoms[i++];
        aj   = iatoms[i++];
        gi   = forceparams[type].waxs_debye.tpi;
        gj   = forceparams[type].waxs_debye.tpj;

        pbc_rvec_sub(pbc, x[ai], x[aj], dx);

        // inner product = dx.dx = |dx|^2
        dr2 = iprod(dx, dx);
        if ((dr2 > rmin2_) && (ai != aj))
        {
            /* Compute scalar force */
            fscal = 0;
            dr_1  = gmx_invsqrt(dr2);
            dr    = dr2*dr_1;
            for (n = 0; (n < Scalc_.size()); n++)
            {
                qdr    = q_[n]*dr;
                if (qdr > 0)
                {
                    fscal += (4.0*A_[n]*alpha_)*(F_[gi][n]*F_[gj][n])*(cos(qdr)-(sin(qdr)/qdr))*dr_1;
                }
            }

            /* loop over the xyz dimensions */
            for (m = 0; m < DIM; m++)
            {
                fff       = fscal*dx[m]*dr_1;
                f[ai][m] += fff;
                f[aj][m] -= fff;
            }
        }
    }
    if (NULL != nrnb)
    {
        inc_nrnb(nrnb, eNR_WAXS_DEBYE, nbonds*Scalc_.size());
    }
    return vtot;
}

void WaxsDebyeForce::dump()
{
    unsigned int n;

    if (step_ % nstout_ == 0)
    {
        if (fp_out_)
        {
            fprintf(fp_out_->handle(), "@type xy\n");
            for (n = 0; n < Scalc_.size(); n++)
            {
                fprintf(fp_out_->handle(), "%8.3f  %8.3f  %8.3f\n",
                        q_[n], Scalc_[n], alpha_*(Scalc_[n]-Sref_[n]));
            }
            fprintf(fp_out_->handle(), "&\n");
            fflush(fp_out_->handle());
        }
        if (fp_alpha_)
        {
            fprintf(fp_alpha_->handle(), "%8d  %8.5f  %8.3f\n", step_, alpha_, alpha_dot_);
            fflush(fp_alpha_->handle());
        }
    }
}

void WaxsDebyeForce::updateAlpha()
{
    switch (alpha_mode_)
    {
        case 0:
        {
            /* no optimization */
            alpha_ = 0.5 * (alpha_min_ + alpha_max_);
        }
        break;

        case 1:
        {
            /* golden section minimization */
            double       a = alpha_min_, c = alpha_max_, b;
            double       x, fb, fx;
            unsigned int n, i;

            alpha_ = 0.5 * (alpha_min_ + alpha_max_);

            /* best convergence with golden ratio */
            b = (c+1.61803*a)/(1.61803+1);

            /* convergence parameters, can be hard-coded I think */
            unsigned int max_iterations = 10;
            float        tolerance      = .01;

            /* initialize a,b,c and fb=chi_sq(b). fa and fc are not needed! */
            fb = 0;
            for (n = 0; n < q_.size(); n++)
            {
                fb += sqr(dSexp_[n] - b*(Scalc_[n] - Sref_[n]));
            }

            fprintf(stderr, "\nalpha=");
            for (i = 0; i < max_iterations; i++)
            {
                fprintf(stderr, " %f", alpha_);

                /* x also comes from the golden ratio, and is in the biggest interval */
                x = c - b + a;

                /* evaluate fx = chi_sq(alpha=x) */
                fx = 0;
                for (n = 0; n < q_.size(); n++)
                {
                    fx += sqr(dSexp_[n] - x*(Scalc_[n] - Sref_[n]));
                }

                /* four possibilies for new triplet: b<x or b>x on left and right */
                bool left = (b-a > c-b);

                if (fb <= fx)
                {
                    if (left)
                    {
                        /* (x,b,c) */
                        a      = x;
                        alpha_ = 0.5*(x+c);
                    }
                    else
                    {
                        /* (a,b,x) */
                        c      = x;
                        alpha_ = 0.5*(a+x);
                    }
                }
                else
                {
                    if (left)
                    {
                        /* (a,x,b) */
                        c      = b;
                        b      = x;
                        fb     = fx;
                        alpha_ = 0.5*(a+b);
                    }
                    else
                    {
                        /* (b,x,c) */
                        a      = b;
                        b      = x;
                        fb     = fx;
                        alpha_ = 0.5*(b+c);
                    }
                }

                if (c-a < tolerance)
                {
                    break;
                }
            }
        }
        break;

        case 2:
        {
            /* use alpha-dynamics */
            double tmp;
            alpha_dot_ = alpha_dot_ + Falpha_/Malpha_*alpha_step_;
            tmp        = alpha_ + alpha_dot_*alpha_step_;

            if (tmp < alpha_min_)
            {
                alpha_     = alpha_min_;
                alpha_dot_ = 0.0;
            }
            else if (tmp > alpha_max_)
            {
                alpha_     = alpha_max_;
                alpha_dot_ = 0.0;
            }
            else
            {
                alpha_ = tmp;
            }
        }
        default:
            GMX_THROW(InvalidInputError("alpha_mode_ should be 0, 1, or 2"));
    }
}

}

//! C interface to the C++ core routines for WAXS calculations.
typedef struct waxs_debye_force_c
{
    //! The C++ class
    gmx::WaxsDebyeForce *wdf;
} t_waxs_debye_force;

waxs_debye_force_t waxs_debye_force_init(FILE             *fplog,
                                         const char       *sfactor,
                                         const char       *waxs_ref,
                                         const char       *waxs_diff,
                                         const char       *waxs_out,
                                         const char       *waxs_alpha,
                                         const t_commrec  *cr,
                                         const t_inputrec *ir)
{
    waxs_debye_force_t wdf = NULL;

    if (ir->waxs.waxs_type == eWaxsDebye)
    {
        snew(wdf, 1);
        try
        {
            wdf->wdf = new gmx::WaxsDebyeForce(fplog, sfactor, waxs_ref, waxs_diff,
                                               waxs_out, waxs_alpha, cr, ir);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    return wdf;
}

void waxs_debye_force_done(waxs_debye_force_t wdf)
{
    if (NULL != wdf)
    {
        try
        {
            wdf->wdf->finalize();
            delete wdf->wdf;
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    sfree(wdf);
}

void waxs_debye_force_update_alpha(waxs_debye_force_t wdf)
{
    if (NULL != wdf)
    {
        try
        {
            wdf->wdf->updateAlpha();
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

real waxs_debye_force_calc(FILE              *fplog,
                           waxs_debye_force_t wdf,
                           const int          nbonds,
                           const t_iatom      iatoms[],
                           const t_iparams    forceparams[],
                           rvec               x[],
                           rvec               f[],
                           const t_pbc       *pbc,
                           const t_commrec   *cr,
                           t_nrnb            *nrnb)
{
    real ener = 0.0;

    if (NULL != wdf)
    {
        try
        {
            ener = wdf->wdf->calc(fplog, nbonds, iatoms, forceparams, x, f,
                                  pbc, cr, nrnb);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    return ener;
}

/* End of C interface */
