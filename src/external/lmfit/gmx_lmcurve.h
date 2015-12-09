/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.h
 *
 * Contents:  Declarations for Levenberg-Marquardt curve fitting.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#ifndef GMX_LMCURVE_H
#define GMX_LMCURVE_H

#include "gmx_lmstruct.h"

void gmx_lmcurve( const int n_par, double *par, const int m_dat,
                  const double *t, const double *y, const double *dy,
                  double (*f)(const double t, const double *par ),
                  const lm_control_struct *control,
                  lm_status_struct *status );

#endif /* GMX_LMCURVE_H */
