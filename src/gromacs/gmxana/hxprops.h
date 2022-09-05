/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_GMXANA_HXPROPS_H
#define GMX_GMXANA_HXPROPS_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_atom;
struct t_resinfo;

#define PHI_AHX (-55.0)
#define PSI_AHX (-45.0)
/* Canonical values of the helix phi/psi angles */


/*! \internal \brief Struct containing properties of a residue in a protein backbone. */
struct t_bb
{
    //! Protein backbone phi angle.
    real phi;
    //! Protein backbone psi angle.
    real psi;
    //! RMS distance of phi and psi angles from ideal helix
    real pprms2;
    //! Estimated J-coupling value
    real jcaha;
    //! Value of 3 turn helix?
    real d3;
    //! Value of 4 turn helix?
    real d4;
    //! Value of 5 turn?
    real d5;
    //! Average of RMS for analysis.
    real rmsa;
    //! If the structure is helical.
    gmx_bool bHelix;
    //! Number of elliptical elements
    int nhx;
    //! Average RMS Deviation when atoms of this residue are fitted to ideal helix
    int nrms;
    //! Residue index for output, relative to gmx_helix -r0 value
    int resno;
    //! Index for previous carbon.
    int Cprev;
    //! Index for backbone nitrogen.
    int N;
    //! Index for backbone NH hydrogen.
    int H;
    //! Index for alpha carbon.
    int CA;
    //! Index for carbonyl carbon.
    int C;
    //! Index for carbonyl oxygen.
    int O;
    //! Index for next backbone nitrogen.
    int Nnext;
    //! Name for this residue.
    char label[32];
};

enum
{
    efhRAD,
    efhTWIST,
    efhRISE,
    efhLEN,
    efhDIP,
    efhRMS,
    efhRMSA,
    efhCD222,
    efhPPRMS,
    efhCPHI,
    efhPHI,
    efhPSI,
    efhHB3,
    efhHB4,
    efhHB5,
    efhJCA,
    efhAHX,
    efhNR
};

extern real ahx_len(int gnx, const int index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern real ellipticity(int nres, t_bb bb[]);

extern real radius(FILE* fp, int nca, const int ca_index[], rvec x[]);
/* Assume we have calphas */

extern real twist(int nca, const int caindex[], rvec x[]);
/* Calculate the twist of the helix */

extern real pprms(FILE* fp, int nbb, t_bb bb[]);
/* Calculate the average RMS from canonical phi/psi values
 * and the distance per residue
 */

extern real ca_phi(int gnx, const int index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern real dip(int nbb, const int bbind[], const rvec x[], const t_atom atom[]);

extern real rise(int gnx, const int index[], rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern void
av_hblen(FILE* fp3, FILE* fp3a, FILE* fp4, FILE* fp4a, FILE* fp5, FILE* fp5a, real t, int nres, t_bb bb[]);

extern void av_phipsi(FILE* fphi, FILE* fpsi, FILE* fphi2, FILE* fpsi2, real t, int nres, t_bb bb[]);

/*! \brief Allocate and fill an array of information about residues in a protein backbone.
 *
 * The user is propted for an index group of protein residues (little
 * error checking occurs). For the number of residues found in the
 * selected group, nbb entries are made in the returned array.  Each
 * entry contains the atom indices of the N, H, CA, C and O atoms (for
 * PRO, H means CD), as well as the C of the previous residue and the
 * N of the next (-1 if not found).
 *
 * In the output array, the first residue will be numbered starting
 * from res0. */
extern t_bb* mkbbind(const char* fn,
                     int*        nres,
                     int*        nbb,
                     int         res0,
                     int*        nall,
                     int**       index,
                     char***     atomname,
                     t_atom      atom[],
                     t_resinfo*  resinfo);

extern void do_start_end(int      nres,
                         t_bb     bb[],
                         int*     nbb,
                         int      bbindex[],
                         int*     nca,
                         int      caindex[],
                         gmx_bool bRange,
                         int      rStart,
                         int      rEnd);

extern void calc_hxprops(int nres, t_bb bb[], const rvec x[]);

extern void pr_bb(FILE* fp, int nres, t_bb bb[]);

#endif
