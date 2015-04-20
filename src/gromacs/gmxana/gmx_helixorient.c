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

#include <math.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

int gmx_helixorient(int argc, char *argv[])
{
    const char      *desc[] = {
        "[THISMODULE] calculates the coordinates and direction of the average",
        "axis inside an alpha helix, and the direction/vectors of both the",
        "C[GRK]alpha[grk] and (optionally) a sidechain atom relative to the axis.[PAR]",
        "As input, you need to specify an index group with C[GRK]alpha[grk] atoms",
        "corresponding to an [GRK]alpha[grk]-helix of continuous residues. Sidechain",
        "directions require a second index group of the same size, containing",
        "the heavy atom in each residue that should represent the sidechain.[PAR]",
        "[BB]Note[bb] that this program does not do any fitting of structures.[PAR]",
        "We need four C[GRK]alpha[grk] coordinates to define the local direction of the helix",
        "axis.[PAR]",
        "The tilt/rotation is calculated from Euler rotations, where we define",
        "the helix axis as the local [IT]x[it]-axis, the residues/C[GRK]alpha[grk] vector as [IT]y[it], and the",
        "[IT]z[it]-axis from their cross product. We use the Euler Y-Z-X rotation, meaning",
        "we first tilt the helix axis (1) around and (2) orthogonal to the residues",
        "vector, and finally apply the (3) rotation around it. For debugging or other",
        "purposes, we also write out the actual Euler rotation angles as [TT]theta[1-3].xvg[tt]"
    };

    t_topology      *top = NULL;
    real             t;
    rvec            *x = NULL, dx;
    matrix           box;
    t_trxstatus     *status;
    int              natoms;
    real             theta1, theta2, theta3;

    int              d, i, j, teller = 0;
    int              iCA, iSC;
    atom_id         *ind_CA;
    atom_id         *ind_SC;
    char            *gn_CA;
    char            *gn_SC;
    rvec             averageaxis;
    rvec             v1, v2, p1, p2, vtmp, vproj;
    rvec            *x_CA, *x_SC;
    rvec            *r12;
    rvec            *r23;
    rvec            *r34;
    rvec            *diff13;
    rvec            *diff24;
    rvec            *helixaxis;
    rvec            *residuehelixaxis;
    rvec            *residueorigin;
    rvec            *residuevector;
    rvec            *sidechainvector;

    rvec             axes_t0[3];
    rvec             axes[3];
    rvec            *residuehelixaxis_t0;
    rvec            *residuevector_t0;
    rvec            *axis3_t0;
    rvec            *residuehelixaxis_tlast;
    rvec            *residuevector_tlast;
    rvec            *axis3_tlast;
    rvec             refaxes[3], newaxes[3];
    rvec             unitaxes[3];
    rvec             rot_refaxes[3], rot_newaxes[3];

    real             tilt, rotation;
    rvec            *axis3;
    real            *twist, *residuetwist;
    real            *radius, *residueradius;
    real            *rise, *residuerise;
    real            *residuebending;

    real             tmp, rotangle;
    real             weight[3];
    t_pbc            pbc;
    matrix           A;

    FILE            *fpaxis, *fpcenter, *fptilt, *fprotation;
    FILE            *fpradius, *fprise, *fptwist;
    FILE            *fptheta1, *fptheta2, *fptheta3;
    FILE            *fpbending;
    int              ePBC;

    output_env_t     oenv;
    gmx_rmpbc_t      gpbc = NULL;

    static  gmx_bool bSC          = FALSE;
    static gmx_bool  bIncremental = FALSE;

    static t_pargs   pa[] = {
        { "-sidechain",      FALSE, etBOOL, {&bSC},
          "Calculate sidechain directions relative to helix axis too." },
        { "-incremental",        FALSE, etBOOL, {&bIncremental},
          "Calculate incremental rather than total rotation/tilt." },
    };
#define NPA asize(pa)

    t_filenm fnm[] = {
        { efTPR, NULL, NULL, ffREAD },
        { efTRX, "-f", NULL, ffREAD },
        { efNDX, NULL, NULL, ffOPTRD },
        { efDAT, "-oaxis",    "helixaxis", ffWRITE },
        { efDAT, "-ocenter",  "center", ffWRITE },
        { efXVG, "-orise",    "rise", ffWRITE },
        { efXVG, "-oradius",  "radius", ffWRITE },
        { efXVG, "-otwist",   "twist", ffWRITE },
        { efXVG, "-obending", "bending", ffWRITE },
        { efXVG, "-otilt",    "tilt", ffWRITE },
        { efXVG, "-orot",     "rotation", ffWRITE }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC);

    for (i = 0; i < 3; i++)
    {
        weight[i] = 1.0;
    }

    /* read index files */
    printf("Select a group of Calpha atoms corresponding to a single continuous helix:\n");
    get_index(&(top->atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &iCA, &ind_CA, &gn_CA);
    snew(x_CA, iCA);
    snew(x_SC, iCA); /* sic! */

    snew(r12, iCA-3);
    snew(r23, iCA-3);
    snew(r34, iCA-3);
    snew(diff13, iCA-3);
    snew(diff24, iCA-3);
    snew(helixaxis, iCA-3);
    snew(twist, iCA);
    snew(residuetwist, iCA);
    snew(radius, iCA);
    snew(residueradius, iCA);
    snew(rise, iCA);
    snew(residuerise, iCA);
    snew(residueorigin, iCA);
    snew(residuehelixaxis, iCA);
    snew(residuevector, iCA);
    snew(sidechainvector, iCA);
    snew(residuebending, iCA);
    snew(residuehelixaxis_t0, iCA);
    snew(residuevector_t0, iCA);
    snew(axis3_t0, iCA);
    snew(residuehelixaxis_tlast, iCA);
    snew(residuevector_tlast, iCA);
    snew(axis3_tlast, iCA);
    snew(axis3, iCA);

    if (bSC)
    {
        printf("Select a group of atoms defining the sidechain direction (1/residue):\n");
        get_index(&(top->atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &iSC, &ind_SC, &gn_SC);
        if (iSC != iCA)
        {
            gmx_fatal(FARGS, "Number of sidechain atoms (%d) != number of CA atoms (%d)", iSC, iCA);
        }

    }

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);

    fpaxis    = gmx_ffopen(opt2fn("-oaxis", NFILE, fnm), "w");
    fpcenter  = gmx_ffopen(opt2fn("-ocenter", NFILE, fnm), "w");
    fprise    = gmx_ffopen(opt2fn("-orise", NFILE, fnm), "w");
    fpradius  = gmx_ffopen(opt2fn("-oradius", NFILE, fnm), "w");
    fptwist   = gmx_ffopen(opt2fn("-otwist", NFILE, fnm), "w");
    fpbending = gmx_ffopen(opt2fn("-obending", NFILE, fnm), "w");

    fptheta1 = gmx_ffopen("theta1.xvg", "w");
    fptheta2 = gmx_ffopen("theta2.xvg", "w");
    fptheta3 = gmx_ffopen("theta3.xvg", "w");

    if (bIncremental)
    {
        fptilt = xvgropen(opt2fn("-otilt", NFILE, fnm),
                          "Incremental local helix tilt", "Time(ps)", "Tilt (degrees)",
                          oenv);
        fprotation = xvgropen(opt2fn("-orot", NFILE, fnm),
                              "Incremental local helix rotation", "Time(ps)",
                              "Rotation (degrees)", oenv);
    }
    else
    {
        fptilt = xvgropen(opt2fn("-otilt", NFILE, fnm),
                          "Cumulative local helix tilt", "Time(ps)", "Tilt (degrees)", oenv);
        fprotation = xvgropen(opt2fn("-orot", NFILE, fnm),
                              "Cumulative local helix rotation", "Time(ps)",
                              "Rotation (degrees)", oenv);
    }

    clear_rvecs(3, unitaxes);
    unitaxes[0][0] = 1;
    unitaxes[1][1] = 1;
    unitaxes[2][2] = 1;

    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);

    do
    {
        /* initialisation for correct distance calculations */
        set_pbc(&pbc, ePBC, box);
        /* make molecules whole again */
        gmx_rmpbc(gpbc, natoms, box, x);

        /* copy coords to our smaller arrays */
        for (i = 0; i < iCA; i++)
        {
            copy_rvec(x[ind_CA[i]], x_CA[i]);
            if (bSC)
            {
                copy_rvec(x[ind_SC[i]], x_SC[i]);
            }
        }

        for (i = 0; i < iCA-3; i++)
        {
            rvec_sub(x_CA[i+1], x_CA[i], r12[i]);
            rvec_sub(x_CA[i+2], x_CA[i+1], r23[i]);
            rvec_sub(x_CA[i+3], x_CA[i+2], r34[i]);
            rvec_sub(r12[i], r23[i], diff13[i]);
            rvec_sub(r23[i], r34[i], diff24[i]);
            /* calculate helix axis */
            cprod(diff13[i], diff24[i], helixaxis[i]);
            svmul(1.0/norm(helixaxis[i]), helixaxis[i], helixaxis[i]);

            tmp       = cos_angle(diff13[i], diff24[i]);
            twist[i]  = 180.0/M_PI * acos( tmp );
            radius[i] = sqrt( norm(diff13[i])*norm(diff24[i]) ) / (2.0* (1.0-tmp) );
            rise[i]   = fabs(iprod(r23[i], helixaxis[i]));

            svmul(radius[i]/norm(diff13[i]), diff13[i], v1);
            svmul(radius[i]/norm(diff24[i]), diff24[i], v2);

            rvec_sub(x_CA[i+1], v1, residueorigin[i+1]);
            rvec_sub(x_CA[i+2], v2, residueorigin[i+2]);
        }
        residueradius[0] = residuetwist[0] = residuerise[0] = 0;

        residueradius[1] = radius[0];
        residuetwist[1]  = twist[0];
        residuerise[1]   = rise[0];

        residuebending[0] = residuebending[1] = 0;
        for (i = 2; i < iCA-2; i++)
        {
            residueradius[i]  = 0.5*(radius[i-2]+radius[i-1]);
            residuetwist[i]   = 0.5*(twist[i-2]+twist[i-1]);
            residuerise[i]    = 0.5*(rise[i-2]+rise[i-1]);
            residuebending[i] = 180.0/M_PI*acos( cos_angle(helixaxis[i-2], helixaxis[i-1]) );
        }
        residueradius[iCA-2]  = radius[iCA-4];
        residuetwist[iCA-2]   = twist[iCA-4];
        residuerise[iCA-2]    = rise[iCA-4];
        residueradius[iCA-1]  = residuetwist[iCA-1] = residuerise[iCA-1] = 0;
        residuebending[iCA-2] = residuebending[iCA-1] = 0;

        clear_rvec(residueorigin[0]);
        clear_rvec(residueorigin[iCA-1]);

        /* average helix axes to define them on the residues.
         * Just extrapolate second first/list atom.
         */
        copy_rvec(helixaxis[0], residuehelixaxis[0]);
        copy_rvec(helixaxis[0], residuehelixaxis[1]);

        for (i = 2; i < iCA-2; i++)
        {
            rvec_add(helixaxis[i-2], helixaxis[i-1], residuehelixaxis[i]);
            svmul(0.5, residuehelixaxis[i], residuehelixaxis[i]);
        }
        copy_rvec(helixaxis[iCA-4], residuehelixaxis[iCA-2]);
        copy_rvec(helixaxis[iCA-4], residuehelixaxis[iCA-1]);

        /* Normalize the axis */
        for (i = 0; i < iCA; i++)
        {
            svmul(1.0/norm(residuehelixaxis[i]), residuehelixaxis[i], residuehelixaxis[i]);
        }

        /* calculate vector from origin to residue CA */
        fprintf(fpaxis, "%15.12g  ", t);
        fprintf(fpcenter, "%15.12g  ", t);
        fprintf(fprise, "%15.12g  ", t);
        fprintf(fpradius, "%15.12g  ", t);
        fprintf(fptwist, "%15.12g  ", t);
        fprintf(fpbending, "%15.12g  ", t);

        for (i = 0; i < iCA; i++)
        {
            if (i == 0 || i == iCA-1)
            {
                fprintf(fpaxis, "%15.12g %15.12g %15.12g       ", 0.0, 0.0, 0.0);
                fprintf(fpcenter, "%15.12g %15.12g %15.12g       ", 0.0, 0.0, 0.0);
                fprintf(fprise, "%15.12g  ", 0.0);
                fprintf(fpradius, "%15.12g  ", 0.0);
                fprintf(fptwist, "%15.12g  ", 0.0);
                fprintf(fpbending, "%15.12g  ", 0.0);
            }
            else
            {
                rvec_sub( bSC ? x_SC[i] : x_CA[i], residueorigin[i], residuevector[i]);
                svmul(1.0/norm(residuevector[i]), residuevector[i], residuevector[i]);
                cprod(residuehelixaxis[i], residuevector[i], axis3[i]);
                fprintf(fpaxis, "%15.12g %15.12g %15.12g       ", residuehelixaxis[i][0], residuehelixaxis[i][1], residuehelixaxis[i][2]);
                fprintf(fpcenter, "%15.12g %15.12g %15.12g       ", residueorigin[i][0], residueorigin[i][1], residueorigin[i][2]);

                fprintf(fprise, "%15.12g  ", residuerise[i]);
                fprintf(fpradius, "%15.12g  ", residueradius[i]);
                fprintf(fptwist, "%15.12g  ", residuetwist[i]);
                fprintf(fpbending, "%15.12g  ", residuebending[i]);

                /* angle with local vector? */
                /*
                   printf("res[%2d]:  axis: %g %g %g    origin: %g %g %g   vector: %g %g %g   angle: %g\n",i,
                       residuehelixaxis[i][0],
                       residuehelixaxis[i][1],
                       residuehelixaxis[i][2],
                       residueorigin[i][0],
                       residueorigin[i][1],
                       residueorigin[i][2],
                       residuevector[i][0],
                       residuevector[i][1],
                       residuevector[i][2],
                       180.0/M_PI*acos( cos_angle(residuevector[i],residuehelixaxis[i]) ));
                 */
                /*      fprintf(fp,"%15.12g %15.12g %15.12g   %15.12g %15.12g %15.12g\n",
                              residuehelixaxis[i][0],
                              residuehelixaxis[i][1],
                              residuehelixaxis[i][2],
                              residuevector[i][0],
                              residuevector[i][1],
                              residuevector[i][2]);
                 */
            }
        }
        fprintf(fprise, "\n");
        fprintf(fpradius, "\n");
        fprintf(fpaxis, "\n");
        fprintf(fpcenter, "\n");
        fprintf(fptwist, "\n");
        fprintf(fpbending, "\n");

        if (teller == 0)
        {
            for (i = 0; i < iCA; i++)
            {
                copy_rvec(residuehelixaxis[i], residuehelixaxis_t0[i]);
                copy_rvec(residuevector[i], residuevector_t0[i]);
                copy_rvec(axis3[i], axis3_t0[i]);
            }
        }
        else
        {
            fprintf(fptilt, "%15.12g       ", t);
            fprintf(fprotation, "%15.12g       ", t);
            fprintf(fptheta1, "%15.12g      ", t);
            fprintf(fptheta2, "%15.12g      ", t);
            fprintf(fptheta3, "%15.12g      ", t);

            for (i = 0; i < iCA; i++)
            {
                if (i == 0 || i == iCA-1)
                {
                    tilt = rotation = 0;
                }
                else
                {
                    if (!bIncremental)
                    {
                        /* Total rotation & tilt */
                        copy_rvec(residuehelixaxis_t0[i], refaxes[0]);
                        copy_rvec(residuevector_t0[i], refaxes[1]);
                        copy_rvec(axis3_t0[i], refaxes[2]);
                    }
                    else
                    {
                        /* Rotation/tilt since last step */
                        copy_rvec(residuehelixaxis_tlast[i], refaxes[0]);
                        copy_rvec(residuevector_tlast[i], refaxes[1]);
                        copy_rvec(axis3_tlast[i], refaxes[2]);
                    }
                    copy_rvec(residuehelixaxis[i], newaxes[0]);
                    copy_rvec(residuevector[i], newaxes[1]);
                    copy_rvec(axis3[i], newaxes[2]);

                    /*
                       printf("frame %d, i=%d:\n  old: %g %g %g , %g %g %g , %g %g %g\n  new:  %g %g %g , %g %g %g , %g %g %g\n",
                       teller,i,
                       refaxes[0][0],refaxes[0][1],refaxes[0][2],
                               refaxes[1][0],refaxes[1][1],refaxes[1][2],
                               refaxes[2][0],refaxes[2][1],refaxes[2][2],
                               newaxes[0][0],newaxes[0][1],newaxes[0][2],
                               newaxes[1][0],newaxes[1][1],newaxes[1][2],
                               newaxes[2][0],newaxes[2][1],newaxes[2][2]);
                     */

                    /* rotate reference frame onto unit axes */
                    calc_fit_R(3, 3, weight, unitaxes, refaxes, A);
                    for (j = 0; j < 3; j++)
                    {
                        mvmul(A, refaxes[j], rot_refaxes[j]);
                        mvmul(A, newaxes[j], rot_newaxes[j]);
                    }

                    /* Determine local rotation matrix A */
                    calc_fit_R(3, 3, weight, rot_newaxes, rot_refaxes, A);
                    /* Calculate euler angles, from rotation order y-z-x, where
                     * x is helixaxis, y residuevector, and z axis3.
                     *
                     * A contains rotation column vectors.
                     */

                    /*
                       printf("frame %d, i=%d, A: %g %g %g  , %g %g %g , %g %g %g\n",
                       teller,i,A[0][0],A[0][1],A[0][2],A[1][0],A[1][1],A[1][2],A[2][0],A[2][1],A[2][2]);
                     */

                    theta1 = 180.0/M_PI*atan2(A[0][2], A[0][0]);
                    theta2 = 180.0/M_PI*asin(-A[0][1]);
                    theta3 = 180.0/M_PI*atan2(A[2][1], A[1][1]);

                    tilt     = sqrt(theta1*theta1+theta2*theta2);
                    rotation = theta3;
                    fprintf(fptheta1, "%15.12g  ", theta1);
                    fprintf(fptheta2, "%15.12g  ", theta2);
                    fprintf(fptheta3, "%15.12g  ", theta3);

                }
                fprintf(fptilt, "%15.12g  ", tilt);
                fprintf(fprotation, "%15.12g  ", rotation);
            }
            fprintf(fptilt, "\n");
            fprintf(fprotation, "\n");
            fprintf(fptheta1, "\n");
            fprintf(fptheta2, "\n");
            fprintf(fptheta3, "\n");
        }

        for (i = 0; i < iCA; i++)
        {
            copy_rvec(residuehelixaxis[i], residuehelixaxis_tlast[i]);
            copy_rvec(residuevector[i], residuevector_tlast[i]);
            copy_rvec(axis3[i], axis3_tlast[i]);
        }

        teller++;
    }
    while (read_next_x(oenv, status, &t, x, box));

    gmx_rmpbc_done(gpbc);

    gmx_ffclose(fpaxis);
    gmx_ffclose(fpcenter);
    xvgrclose(fptilt);
    xvgrclose(fprotation);
    gmx_ffclose(fprise);
    gmx_ffclose(fpradius);
    gmx_ffclose(fptwist);
    gmx_ffclose(fpbending);
    gmx_ffclose(fptheta1);
    gmx_ffclose(fptheta2);
    gmx_ffclose(fptheta3);

    close_trj(status);

    return 0;
}
