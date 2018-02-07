/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef FCA_ENTROPY_BASIC_HISTO_H
#define FCA_ENTROPY_BASIC_HISTO_H

/*
 * $Id: g_project.c,v 1.5 2004/11/12 17:17:36 olange Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to  consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */
//static char *SRCID_g_anaeig_c = "$Id: g_project.c,v 1.5 2004/11/12 17:17:36 olange Exp $";
/*! \internal \file
 * \brief
 * Declares FcaEntropyHist.
 */
#include "gromacs/utility/real.h"

#include "fca_minimizing.h"

namespace gmx
{

struct FcaBasicEntropy
{
    real log_kappa1; /* proportionality for binwidth */
    real log_kappa2; /* prop. for 2D bandwith */
    int  ndata;
};

class FcaEntropyHisto
{

    /*! \brief
     * results:
     * S -- resulting 1D entropy
     * input :
     * s -- data-series of length entr->ndata
     * entr -- entropy object
     */
    static real entropy_1D_basic_histo(real* pS, const real s[], const FcaBasicEntropy &entr);

    static real entropy_2D_basic_histo(const real sx[], const real sy[], const FcaBasicEntropy &entr);

    public:

        static void Basic_compute_MI_matrix(FcaMaster* fca, const real log_kappa1, const real log_kappa2);

};

} //gmx namespace

#endif
