/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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

#ifndef _hxprops_h
#define _hxprops_h

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"

#define PHI_AHX (-55.0)
#define PSI_AHX (-45.0)
/* Canonical values of the helix phi/psi angles */


typedef struct {
    real     phi, psi, pprms2;
    real     jcaha;
    real     d3, d4, d5, rmsa;
    gmx_bool bHelix;
    int      nhx;
    int      nrms, resno;
    int      Cprev, N, H, CA, C, O, Nnext;
    char     label[32];
} t_bb;

enum {
    efhRAD,  efhTWIST, efhRISE, efhLEN,
    efhDIP,  efhRMS,   efhRMSA, efhCD222,
    efhPPRMS, efhCPHI,  efhPHI,  efhPSI,
    efhHB3,  efhHB4,   efhHB5,  efhJCA,
    efhAHX,  efhNR
};

extern real ahx_len(int gnx, atom_id index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern real ellipticity(int nres, t_bb bb[]);

extern real radius(FILE *fp, int nca, atom_id ca_index[], rvec x[]);
/* Assume we have calphas */

extern real twist(int nca, atom_id caindex[], rvec x[]);
/* Calculate the twist of the helix */

extern real pprms(FILE *fp, int nbb, t_bb bb[]);
/* Calculate the average RMS from canonical phi/psi values
 * and the distance per residue
 */

extern real ca_phi(int gnx, atom_id index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern real dip(int nbb, atom_id bbind[], rvec x[], t_atom atom[]);

extern real rise(int gnx, atom_id index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern void av_hblen(FILE *fp3, FILE *fp3a,
                     FILE *fp4, FILE *fp4a,
                     FILE *fp5, FILE *fp5a,
                     real t, int nres, t_bb bb[]);

extern void av_phipsi(FILE *fphi, FILE *fpsi, FILE *fphi2, FILE *fpsi2,
                      real t, int nres, t_bb bb[]);

extern t_bb *mkbbind(const char *fn, int *nres, int *nbb, int res0,
                     int *nall, atom_id **index,
                     char ***atomname, t_atom atom[],
                     t_resinfo *resinfo);

extern void do_start_end(int nres, t_bb bb[], int *nbb,
                         atom_id bbindex[], int *nca, atom_id caindex[],
                         gmx_bool bRange, int rStart, int rEnd);

extern void calc_hxprops(int nres, t_bb bb[], rvec x[]);

extern void pr_bb(FILE *fp, int nres, t_bb bb[]);

#endif
