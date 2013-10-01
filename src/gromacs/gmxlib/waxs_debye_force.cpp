/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include <iostream>
#include <string.h>
#include "string2.h"
#include "maths.h"
#include "smalloc.h"
#include "vec.h"
#include "pbc.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "nrnb.h"
#include "waxs_debye_force_c.h"
#include "waxs_debye_force.h"
#include "scattering_factors.h"

using namespace std;

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

double integrate_vector(std::vector<double> q,std::vector<double> S)
{
    double sq = 0;
    
    for(unsigned int n=1; (n<q.size()); n++)
    {
        double dsq = (q[n]-q[n-1])*(S[n]+S[n-1])/2;
        sq += dsq;
    }
    return sq;
}

namespace gmx
{
/*! \brief
 * Implements gmx::WaxsDebyeForce
 *
 * \inpublicapi
 * \author Alexander Bjorling <alexander.bjorling@chem.gu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

void WaxsDebyeForce::ReadData(const char *waxs_ref, 
                              const char *waxs_diff)
{
    std::ifstream ifs;
    
    /* Read the experimental data */
    for(int i = 0; (i<=1); i++)
    {
        const char *waxs_data;
        bool bDiff = (bool) i;
        if (bDiff)
        {
            waxs_data = waxs_diff;
        }
        else
        {
            waxs_data = waxs_ref;
        }
        ifs.open(waxs_data, ios_base::in);
        if (ifs.is_open())
        {
            int n = 0;
            while (!ifs.eof())
            {
                double      q, S, dS;
                char        buf[STRLEN];
                std::string line;
                
                std::getline(ifs, line);
                strncpy(buf,line.c_str(),STRLEN);
                strip_comment(buf);
                
                if (sscanf(buf, "%15lf%15lf%15lf", &q, &S, &dS) == 3)
                    
                {
                    if (!bDiff)
                    {
                        _q.push_back(q);
                        _Sref.push_back(S);
                        _dSref.push_back(dS);
                        /* Initiate internal variables to zero */
                        _A.push_back(0);
                        _Scalc.push_back(0);
                        _dSexp.push_back(0);
                        _ddSexp.push_back(0);
                    }
                    else
                    {
                        if (_q[n] != q)
                        {
                            gmx_fatal(FARGS, "Inconsistency between reference data and difference data. Expected q[%d] to be %f but found %f", n, _q[n], q);
                        }
                        _dSexp[n]  = S;
                        _ddSexp[n] = dS;
                    }
                    n++;
                }
            }
        }
        else if (!bDiff)
        {
            gmx_fatal(FARGS, "Error opening input file %s!", waxs_data);
        }
        ifs.close();
    }
}

WaxsDebyeForce::WaxsDebyeForce(FILE       *fplog,
                               const char *sfactor, 
                               const char *waxs_ref, 
                               const char *waxs_diff, 
                               const char *waxs_out, 
                               const char *waxs_alpha, 
                               const t_commrec *cr,
                               const t_inputrec *ir)
{
    /* read the file */
    std::vector<double> tmp;
    unsigned int        m, n;
    int                 atomic_number;

    //! Parameters from mdp file.
    _rmin2  = sqr(ir->waxs.debye_r_min);
    _rmax2  = sqr(ir->waxs.debye_r_max);
    _k_chi  = ir->waxs.kwaxs;

    /* Starting values for the optimization of alpha and other run parameters */
    _alpha_min  = ir->waxs.debye_alpha_min;
    _alpha_max  = ir->waxs.debye_alpha_max;
    _nstout     = ir->waxs.nstout;

    /* Open data files */
    ReadData(waxs_ref, waxs_diff);

    //! Output file for Scalc(q,t)    
    _dump = NULL;
    if ((_nstout > 0) && MASTER(cr))
    {
        _dump = fopen(waxs_out, "w");
    }

    /* Starting values for the optimization of alpha */
    _alpha      = 0.5*(_alpha_min + _alpha_max);
    _alpha_dot  = 0;
    _Falpha     = 0;
    _Malpha     = 1000;  // 1000
    _alpha_step = 5e-4;
    _step       = 0;

    _aa = NULL;
    if ((_alpha_step != 0.0) && MASTER(cr) && (NULL != waxs_alpha))
    {
        _aa = fopen(waxs_alpha, "w");
    }

    /* Read in the table for Cromer-Mann data */
    ScatteringFactorTable *cmt = new ScatteringFactorTable(sfactor);

    /* Fill the _F array for all types 
       allocate space first: atomic_number is one-based so add one to get _F's size. 
    */
    _F.resize(cmt->max_atomic_number()+1);
    _F0.resize(cmt->max_atomic_number()+1);
    for (m = 0; (m < cmt->size()); m++)
    {
        // Fetch the right atomic number from cmt
        atomic_number = cmt->atomic_number(m);

        // Scattering factor at q = 0
        _F0[atomic_number] = cmt->calc(atomic_number, 0);
        
        // loop over q-points for this specific scattering factor & store in tmp
        for (n = 0; n < _Scalc.size(); n++)
        {
            tmp.push_back(cmt->calc(atomic_number, _q[n]));
        }

        // add tmp to the _F array
        _F[atomic_number] = tmp;
        tmp.clear();
    }

    if (NULL != fplog)
    {
        fprintf(fplog, "WAXS: Refinement code initiating.\n");
        fprintf(fplog, "WAXS: Found %d structure factors in %s.\n", 
                (int)_Scalc.size(), waxs_ref);
        if (NULL != waxs_diff)
        {
            fprintf(fplog, "WAXS: Difference data file %s read.\n", waxs_diff);
        }
        fprintf(fplog, "WAXS: There were %d (atomic/residue) scattering factors.\n",
                (int)cmt->size());
    }
    delete cmt;
}

WaxsDebyeForce::~WaxsDebyeForce()
{
    if (NULL != _dump)
    {
        fclose(_dump);
    }
    if (NULL != _aa)
    {
        fclose(_aa);
    }
}

real WaxsDebyeForce::CalcS0(FILE *fplog,
                            const int nbonds, 
                            const t_iatom iatoms[],
                            const t_iparams forceparams[],
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

        //! Should probably check whether gi and gj are less
        if (((0 < gi) && (gi < (int)_F0.size())) &&
            ((0 < gj) && (gj < (int)_F0.size())))
        {
            double dS = _F0[gi]*_F0[gj];
            if (ai == aj)
            {
                S0 += dS;
            }
            else
            {
                //! Since we typically only store half the matrix of pairs
                S0 += 2*dS;
            }
        }
    }

    /* Global sum of _Scalc */
    if (PAR(cr))
    {
        gmx_sumd(1, &S0, cr);
    }
    if (NULL != fplog)
    {
        fprintf(fplog, "WAXS: S(0) = %g.\n",S0);
    }
        
    return S0;
}


real WaxsDebyeForce::Calc(FILE *fplog,
                          const int nbonds, 
                          const t_iatom iatoms[],
                          const t_iparams forceparams[],
                          rvec x[], 
                          rvec f[],
                          const t_pbc *pbc, 
                          const t_commrec *cr,
                          t_nrnb *nrnb)
{
    int          i, m, ai, aj, gi, gj, type;
    unsigned int n;
    real         dr, dr2, fscal, fff, vtot, qdr, dr_1, dS;
    rvec         dx;
    ivec         dt;

    if (0 == _step)
    {
        (void) CalcS0(fplog, nbonds, iatoms, forceparams, cr);
    }
    /* Reset Scalc */
    for (n = 0; n < _Scalc.size(); n++)
    {
        _Scalc[n] = 0;
    }

    /* Compute _Scalc for this node only */
    for (i = 0; (i < nbonds); )
    {
        type = iatoms[i++];
        ai   = iatoms[i++];
        aj   = iatoms[i++];
        gi   = forceparams[type].waxs_debye.tpi;
        gj   = forceparams[type].waxs_debye.tpj;

        pbc_rvec_sub(pbc, x[ai], x[aj], dx);

        dr2 = iprod(dx, dx);
        //! Special case where the distance is zero
        if (ai == aj)
        {
            for (n = 0; (n < _Scalc.size()); n++)
            {
                double f2 = sqr(_F[gi][n]);
                _Scalc[n] += f2;
            }
        }
        //! We only consider particles within a range of distances
        else if ((dr2 > _rmin2) && (dr2 < _rmax2))
        {
            dr = dr2*gmx_invsqrt(dr2);
            for (n = 0; (n < _Scalc.size()); n++)
            {
                qdr        = _q[n]*dr;
                double fgigj = 2.0*_F[gi][n]*_F[gj][n];
                _Scalc[n] += fgigj*sin(qdr)/(qdr);
            }
        }
    }

    /* Global sum of _Scalc */
    if (PAR(cr))
    {
        gmx_sumd(_Scalc.size(), &(_Scalc[0]), cr);
    }
        
    /* Construction to export Scalc(q,t) to a textfile */
    Dump();

    /* Now we calculate _A, the factor 
     * _A = k_chi/(2 SigmaQ^2) * (dSexp - alpha*dScalc)
     * that we will need repeatedly in the force calculation loop.
     * In addition we compute d chi^2/d alpha
     */
    _Falpha     = 0.0;
    vtot        = 0.0;
    for (n = 0; n < _Scalc.size(); n++)
    {
        dS       = _dSexp[n] - _alpha*(_Scalc[n] - _Sref[n]);
        _A[n]    = _k_chi/(2*sqr(_dSref[n])) * dS;
        vtot    += _A[n]*dS;
        _Falpha += _k_chi * dS/sqr(_dSref[n]) * (_Scalc[n] - _Sref[n]);
    }

    _step++;

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
        if ((dr2 > _rmin2) && (ai != aj))
        {
            /* Compute scalar force */
            fscal = 0;
            dr_1  = gmx_invsqrt(dr2);
            dr    = dr2*dr_1;
            for (n = 0; (n < _Scalc.size()); n++)
            {
                qdr    = _q[n]*dr;
                fscal += (4.0*_A[n]*_alpha)*(_F[gi][n]*_F[gj][n])*(cos(qdr)-(sin(qdr)/qdr))*dr_1;
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
    inc_nrnb(nrnb, eNR_WAXS_DEBYE, nbonds*_Scalc.size());
    
    return vtot;
}

void WaxsDebyeForce::Dump()
{
    unsigned int n;
    FILE *fp;

    if (_step % _nstout == 0)
    {
        if (NULL != _dump)
        {
            fprintf(_dump, "@type xy\n");
            for (n = 0; n < _Scalc.size(); n++)
            {
                fprintf(_dump, "%8.3f  %8.3f  %8.3f\n",
                        _q[n], _alpha*(_Scalc[n]-_Sref[n]), _Scalc[n]);
            }
            fprintf(_dump, "&\n");
            fflush(_dump);
        }
        if (NULL != _aa)
        {
            fprintf(_aa, "%8d  %8.5f  %8.3f\n", _step, _alpha, _alpha_dot);
            fflush(_aa);
        }
    }
}

void WaxsDebyeForce::UpdateAlpha()
{
    double tmp;

    _alpha_dot = _alpha_dot + _Falpha/_Malpha*_alpha_step;
    tmp        = _alpha + _alpha_dot*_alpha_step;

    if (tmp < _alpha_min)
    {
        _alpha     = _alpha_min;
        _alpha_dot = 0.0;
    }
    else if (tmp > _alpha_max)
    {
        _alpha     = _alpha_max;
        _alpha_dot = 0.0;
    }
    else
    {
        _alpha = tmp;
    }

}

}

/* C interface */
typedef struct waxs_debye_force_c 
{
    gmx::WaxsDebyeForce *wdf;
} t_waxs_debye_force;

waxs_debye_force_t waxs_debye_force_init(FILE *fp,
                                         const char *sfactor, 
                                         const char *waxs_ref, 
                                         const char *waxs_diff, 
                                         const char *waxs_out, 
                                         const char *waxs_alpha, 
                                         const t_commrec *cr,
                                         const t_inputrec *ir)
{
    waxs_debye_force_t wdf = NULL;

    if (ir->waxs.waxs_type == eWaxsDebye)
    {
        snew(wdf, 1);
        wdf->wdf = new gmx::WaxsDebyeForce(fp, sfactor, waxs_ref, waxs_diff,
                                           waxs_out, waxs_alpha, cr, ir);
    }
    return wdf;
}

void waxs_debye_force_done(waxs_debye_force_t wdf)
{
    if (NULL != wdf)
    {
        delete wdf->wdf;
    }
    sfree(wdf);
}

void waxs_debye_force_update_alpha(waxs_debye_force_t wdf)
{
    if (NULL != wdf)
    {
        wdf->wdf->UpdateAlpha();
    }
}

real waxs_debye_force_calc(FILE *fplog,
                           waxs_debye_force_t wdf,
                           const int nbonds, 
                           const t_iatom iatoms[],
                           const t_iparams forceparams[],
                           rvec x[], 
                           rvec f[], 
                           const t_pbc *pbc, 
                           const t_commrec *cr,
                           t_nrnb *nrnb)
{
    if (NULL == wdf)
    {
        return 0.0;
    }
    else
    {
        return wdf->wdf->Calc(fplog,nbonds, iatoms, forceparams, x, f,
                              pbc, cr, nrnb);
    }
}

/* End of C interface */
