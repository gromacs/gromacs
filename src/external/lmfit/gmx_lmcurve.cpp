/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.c
 *
 * Contents:  Implements lmcurve, a simplified API for curve fitting
 *            using the generic Levenberg-Marquardt routine lmmin.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 *
 * Note to programmers: Don't patch and fork, but copy and variate!
 *   If you need to compute residues differently, then please do not patch
 * lmcurve.c, but copy it to a differently named file, and rename lmcurve(),
 * lmcurve_evaluate() and lmcurve_data_struct before adapting them to your
 * needs, like we have done in lmcurve_tyd.c.
 */
#include "gmx_lmcurve.h"
#include "gmx_lmmin.h"

typedef struct {
    const double* t;
    const double* y;
    const double* dy;
    double (*f)(const double t, const double* par);
} lmcurve_data_struct;

static void lmcurve_evaluate(
    const double* par, const int m_dat, const void* data, double* fvec,
    int* info)
{
    lmcurve_data_struct* D = (lmcurve_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
    {
        double dy = D->dy[i];
        if (dy == 0)
        {
            dy = 1;
        }
        fvec[i] = (D->y[i] - D->f(D->t[i], par))/dy;
    }
    *info = 0;
}

void gmx_lmcurve(
    const int n_par, double* par, const int m_dat,
    const double* t, const double* y, const double *dy,
    double (*f)(const double t, const double* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_data_struct data = { t, y, dy, f };

    gmx_lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_evaluate,
          control, status);
}
