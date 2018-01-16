/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <m.ghahremanpour@hotmail.com>
 */

#ifndef TUNNING_UTILITY_H
#define TUNNING_UTILITY_H

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_simple_comm.h"
#include "molgen.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"

namespace alexandria
{

void print_stats(FILE        *fp,
                 const char  *prop,
                 gmx_stats_t  lsq,
                 gmx_bool     bHeader,
                 char        *xaxis,
                 char        *yaxis);

void print_lsq_set(FILE *fp, gmx_stats_t lsq);

void xvgr_symbolize(FILE                   *xvgf,
                    int                     nsym,
                    const char             *leg[],
                    const gmx_output_env_t *oenv);

void print_polarizability(FILE              *fp,
                          alexandria::MyMol *mol,
                          char              *calc_name,
                          real               q_toler);

void print_dipole(FILE                      *fp,
                  alexandria::MyMol         *mol,
                  char                      *calc_name,
                  real                       toler);

void print_quadrapole(FILE                  *fp,
                      alexandria::MyMol     *mol,
                      char                  *calc_name,
                      real                   toler);

void print_electric_props(FILE                           *fp,
                          std::vector<alexandria::MyMol>  mymol,
                          const Poldata                  &pd,
                          const char                     *qhisto,
                          const char                     *dipcorr,
                          const char                     *mucorr,
                          const char                     *Qcorr,
                          const char                     *espcorr,
                          const char                     *alphacorr,
                          const char                     *isopolCorr,
                          const char                     *anisopolCorr,
                          real                            dip_toler,
                          real                            quad_toler,
                          real                            alpha_toler,
                          const gmx_output_env_t         *oenv,
                          bool                            bPolar,
                          bool                            bDipole,
                          bool                            bQuadrupole,
                          bool                            bfullTensor,
                          IndexCount                     *indexCount,
                          t_commrec                      *cr,
                          real                            efield);


}

#endif
