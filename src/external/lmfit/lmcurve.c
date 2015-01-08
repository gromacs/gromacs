/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.c
 *
 * Contents:  Levenberg-Marquardt curve-fitting
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#include "lmmin.h"
#include "gromacs/utility/basedefinitions.h"

typedef struct {
    const double *t;
    const double *y;
    const double *dy;
    double (*f)(double t, const double *par);
} lmcurve_data_struct;


void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                       double *fvec, gmx_unused int *info )
{
    int    i;
    double fy;
    lmcurve_data_struct *d = (lmcurve_data_struct*) data;
    for (i = 0; i < m_dat; i++)
    {
        fy      = d->f(d->t[i], par );
        fvec[i] = (d->y[i] - fy)/d->dy[i];
    }
}


void lmcurve( int n_par, double *par, int m_dat,
              const double *t, const double *y, const double *dy,
              double (*f)( double t, const double *par ),
              const lm_control_struct *control,
              lm_status_struct *status )
{
    lmcurve_data_struct data;
    data.t  = t;
    data.y  = y;
    data.dy = dy;
    data.f  = f;

    lmmin( n_par, par, m_dat, (const void*) &data,
           lmcurve_evaluate, control, status );
}
