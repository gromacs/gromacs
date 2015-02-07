/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/3dtransforms.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static void rot_conf(t_atoms *atoms, rvec x[], rvec v[], real trans, real angle,
                     rvec head, rvec tail, int isize, atom_id index[],
                     rvec xout[], rvec vout[])
{
    rvec     arrow, xcm;
    real     theta, phi, arrow_len;
    mat4     Rx, Ry, Rz, Rinvy, Rinvz, Mtot;
    mat4     temp1, temp2, temp3;
    vec4     xv;
    int      i, j, ai;

    rvec_sub(tail, head, arrow);
    arrow_len = norm(arrow);
    if (debug)
    {
        fprintf(debug, "Arrow vector:   %10.4f  %10.4f  %10.4f\n",
                arrow[XX], arrow[YY], arrow[ZZ]);
        fprintf(debug, "Effective translation %g nm\n", trans);
    }
    if (arrow_len == 0.0)
    {
        gmx_fatal(FARGS, "Arrow vector not given");
    }

    /* Copy all aoms to output */
    for (i = 0; (i < atoms->nr); i++)
    {
        copy_rvec(x[i], xout[i]);
        copy_rvec(v[i], vout[i]);
    }

    /* Compute center of mass and move atoms there */
    clear_rvec(xcm);
    for (i = 0; (i < isize); i++)
    {
        rvec_inc(xcm, x[index[i]]);
    }
    for (i = 0; (i < DIM); i++)
    {
        xcm[i] /= isize;
    }
    if (debug)
    {
        fprintf(debug, "Center of mass: %10.4f  %10.4f  %10.4f\n",
                xcm[XX], xcm[YY], xcm[ZZ]);
    }
    for (i = 0; (i < isize); i++)
    {
        rvec_sub(x[index[i]], xcm, xout[index[i]]);
    }

    /* Compute theta and phi that describe the arrow */
    theta = acos(arrow[ZZ]/arrow_len);
    phi   = atan2(arrow[YY]/arrow_len, arrow[XX]/arrow_len);
    if (debug)
    {
        fprintf(debug, "Phi = %.1f, Theta = %.1f\n", RAD2DEG*phi, RAD2DEG*theta);
    }

    /* Now the total rotation matrix: */
    /* Rotate a couple of times */
    gmx_mat4_init_rotation(ZZ, -phi, Rz);
    gmx_mat4_init_rotation(YY, M_PI/2-theta, Ry);
    gmx_mat4_init_rotation(XX, angle*DEG2RAD, Rx);
    Rx[WW][XX] = trans;
    gmx_mat4_init_rotation(YY, theta-M_PI/2, Rinvy);
    gmx_mat4_init_rotation(ZZ, phi, Rinvz);

    gmx_mat4_mmul(temp1, Ry, Rz);
    gmx_mat4_mmul(temp2, Rinvy, Rx);
    gmx_mat4_mmul(temp3, temp2, temp1);
    gmx_mat4_mmul(Mtot, Rinvz, temp3);

    if (debug)
    {
        gmx_mat4_print(debug, "Rz", Rz);
        gmx_mat4_print(debug, "Ry", Ry);
        gmx_mat4_print(debug, "Rx", Rx);
        gmx_mat4_print(debug, "Rinvy", Rinvy);
        gmx_mat4_print(debug, "Rinvz", Rinvz);
        gmx_mat4_print(debug, "Mtot", Mtot);
    }

    for (i = 0; (i < isize); i++)
    {
        ai = index[i];
        gmx_mat4_transform_point(Mtot, xout[ai], xv);
        rvec_add(xv, xcm, xout[ai]);
        gmx_mat4_transform_point(Mtot, v[ai], xv);
        copy_rvec(xv, vout[ai]);
    }
}

int gmx_dyndom(int argc, char *argv[])
{
    const char  *desc[] = {
        "[THISMODULE] reads a [REF].pdb[ref] file output from DynDom",
        "(http://www.cmp.uea.ac.uk/dyndom/).",
        "It reads the coordinates, the coordinates of the rotation axis,",
        "and an index file containing the domains.",
        "Furthermore, it takes the first and last atom of the arrow file",
        "as command line arguments (head and tail) and",
        "finally it takes the translation vector (given in DynDom info file)",
        "and the angle of rotation (also as command line arguments). If the angle",
        "determined by DynDom is given, one should be able to recover the",
        "second structure used for generating the DynDom output.",
        "Because of limited numerical accuracy this should be verified by",
        "computing an all-atom RMSD (using [gmx-confrms]) rather than by file",
        "comparison (using diff).[PAR]",
        "The purpose of this program is to interpolate and extrapolate the",
        "rotation as found by DynDom. As a result unphysical structures with",
        "long or short bonds, or overlapping atoms may be produced. Visual",
        "inspection, and energy minimization may be necessary to",
        "validate the structure."
    };
    static real  trans0 = 0;
    static rvec  head   = { 0, 0, 0 };
    static rvec  tail   = { 0, 0, 0 };
    static real  angle0 = 0, angle1 = 0, maxangle = 0;
    static int   label  = 0, nframes = 11;
    t_pargs      pa[]   = {
        { "-firstangle",    FALSE, etREAL, {&angle0},
          "Angle of rotation about rotation vector" },
        { "-lastangle",    FALSE, etREAL, {&angle1},
          "Angle of rotation about rotation vector" },
        { "-nframe",   FALSE, etINT,  {&nframes},
          "Number of steps on the pathway" },
        { "-maxangle", FALSE, etREAL, {&maxangle},
          "DymDom dtermined angle of rotation about rotation vector" },
        { "-trans",    FALSE, etREAL, {&trans0},
          "Translation (Angstrom) along rotation vector (see DynDom info file)" },
        { "-head",     FALSE, etRVEC, {head},
          "First atom of the arrow vector" },
        { "-tail",     FALSE, etRVEC, {tail},
          "Last atom of the arrow vector" }
    };
    int          i, j, natoms, isize;
    t_trxstatus *status;
    atom_id     *index = NULL, *index_all;
    char         title[256], *grpname;
    t_atoms      atoms;
    real         angle, trans;
    rvec        *x, *v, *xout, *vout;
    matrix       box;
    output_env_t oenv;

    t_filenm     fnm[] = {
        { efPDB, "-f", "dyndom",  ffREAD },
        { efTRO, "-o", "rotated", ffWRITE },
        { efNDX, "-n", "domains", ffREAD }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    if (maxangle == 0)
    {
        gmx_fatal(FARGS, "maxangle not given");
    }

    get_stx_coordnum (opt2fn("-f", NFILE, fnm), &natoms);
    init_t_atoms(&atoms, natoms, TRUE);
    snew(x, natoms);
    snew(v, natoms);
    read_stx_conf(opt2fn("-f", NFILE, fnm), title, &atoms, x, v, NULL, box);
    snew(xout, natoms);
    snew(vout, natoms);

    printf("Select group to rotate:\n");
    rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
    printf("Going to rotate %s containg %d atoms\n", grpname, isize);

    snew(index_all, atoms.nr);
    for (i = 0; (i < atoms.nr); i++)
    {
        index_all[i] = i;
    }

    status = open_trx(opt2fn("-o", NFILE, fnm), "w");

    label = 'A';
    for (i = 0; (i < nframes); i++, label++)
    {
        angle = angle0 + (i*(angle1-angle0))/(nframes-1);
        trans = trans0*0.1*angle/maxangle;
        printf("Frame: %2d (label %c), angle: %8.3f deg., trans: %8.3f nm\n",
               i, label, angle, trans);
        rot_conf(&atoms, x, v, trans, angle, head, tail, isize, index, xout, vout);

        if (label > 'Z')
        {
            label -= 26;
        }
        for (j = 0; (j < atoms.nr); j++)
        {
            atoms.resinfo[atoms.atom[j].resind].chainid = label;
        }

        write_trx(status, atoms.nr, index_all, &atoms, i, angle, box, xout, vout, NULL);
    }
    close_trx(status);

    return 0;
}
