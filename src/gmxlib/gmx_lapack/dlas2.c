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
#include <types/simple.h>

#include "gmx_lapack.h"

void
F77_FUNC(dlas2,DLAS2)(double *f,
       double *g,
       double *h,
       double *ssmin,
       double *ssmax)
{
  double fa = fabs(*f);
  double ga = fabs(*g);
  double ha = fabs(*h);
  double fhmin,fhmax,tmax,tmin,tmp1,tmp2;
  double as,at,au,c;

  fhmin = (fa<ha) ? fa : ha;
  fhmax = (fa>ha) ? fa : ha;
  
  if(fabs(fhmin)<GMX_DOUBLE_MIN) {
    *ssmin = 0.0;
    if(fabs(fhmax)<GMX_DOUBLE_MIN) 
      *ssmax = ga;
    else {
      tmax = (fhmax>ga) ? fhmax : ga;
      tmin = (fhmax<ga) ? fhmax : ga;
      tmp1 = tmin / tmax;
      tmp1 = tmp1 * tmp1;
      *ssmax = tmax*sqrt(1.0 + tmp1);
    }
  } else {
    if(ga<fhmax) {
      as = 1.0 + fhmin / fhmax;
      at = (fhmax-fhmin) / fhmax;
      au = (ga/fhmax);
      au = au * au;
      c = 2.0 / ( sqrt(as*as+au) + sqrt(at*at+au) );
      *ssmin = fhmin * c;
      *ssmax = fhmax / c;
    } else {
      au = fhmax / ga;
      if(fabs(au)<GMX_DOUBLE_MIN) {
	*ssmin = (fhmin*fhmax)/ga;
	*ssmax = ga;
      } else {
	as = 1.0 + fhmin / fhmax;
	at = (fhmax-fhmin)/fhmax;
	tmp1 = as*au;
	tmp2 = at*au;
	c = 1.0 / ( sqrt(1.0+tmp1*tmp1) + sqrt(1.0+tmp2*tmp2));
	*ssmin = (fhmin*c)*au;
	*ssmin = *ssmin + *ssmin;
	*ssmax = ga / (c+c);
      }
    }
  }
  return;
}
