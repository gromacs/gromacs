/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#include <math.h>
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void
F77_FUNC(dlartg,DLARTG)(double *f,
	double *g,
	double *cs,
	double *sn,
	double *r)
{
  double minval,safemin, safemin2, safemx2, eps;
  double f1,g1,f1a,g1a,scale;
  int i,n,count;

  eps = GMX_DOUBLE_EPS;
  minval = GMX_DOUBLE_MIN;
  safemin = minval*(1.0+eps);
  n = 0.5*log( safemin/eps ) / log(2);
  safemin2 = pow(2,n);

  safemx2 = 1.0 / safemin2;

  if(fabs(*g)<GMX_DOUBLE_MIN) {
    *cs = 1.0;
    *sn = 0.0;
    *r = *f;
  } else if (fabs(*f)<GMX_DOUBLE_MIN) {
    *cs = 0.0;
    *sn = 1.0;
    *r = *g;
  } else {
    f1 = *f;
    g1 = *g;
    f1a = fabs(f1);
    g1a = fabs(g1);
    scale = (f1a > g1a) ? f1a : g1a;
    if(scale >= safemx2) {
      count = 0;
      while(scale >= safemx2) {
	count++;
	f1 *= safemin2;
	g1 *= safemin2;
	f1a = fabs(f1);
	g1a = fabs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r = sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemx2;
    } else if (scale<=safemin2) {
      count = 0;
      while(scale <= safemin2) {
	count++;
	f1 *= safemx2;
	g1 *= safemx2;
	f1a = fabs(f1);
	g1a = fabs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r = sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemin2;
    } else {
      *r = sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
    }
    if(fabs(*f)>fabs(*g) && *cs<0.0) {
      *cs *= -1.0;
      *sn *= -1.0;
      *r  *= -1.0;
    }
  }
  return;
}
      
