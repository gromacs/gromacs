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
#include "gmxpre.h"

#include "gromacs/fileio/trnio.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static char *RotStr = {"Enforced rotation:"};


static char s_vec[STRLEN];


static void string2dvec(char buf[], dvec nums)
{
    if (sscanf(buf, "%lf%lf%lf", &nums[0], &nums[1], &nums[2]) != 3)
    {
        gmx_fatal(FARGS, "Expected three numbers at input line %s", buf);
    }
}


extern char **read_rotparams(int *ninp_p, t_inpfile **inp_p, t_rot *rot,
                             warninp_t wi)
{
    int         ninp, g, m;
    t_inpfile  *inp;
    const char *tmp;
    char      **grpbuf;
    char        buf[STRLEN];
    char        warn_buf[STRLEN];
    dvec        vec;
    t_rotgrp   *rotg;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read rotation parameters */
    CTYPE("Output frequency for angle, torque and rotation potential energy for the whole group");
    ITYPE("rot-nstrout",     rot->nstrout, 100);
    CTYPE("Output frequency for per-slab data (angles, torques and slab centers)");
    ITYPE("rot-nstsout",     rot->nstsout, 1000);
    CTYPE("Number of rotation groups");
    ITYPE("rot-ngroups",     rot->ngrp, 1);

    if (rot->ngrp < 1)
    {
        gmx_fatal(FARGS, "rot-ngroups should be >= 1");
    }

    snew(rot->grp, rot->ngrp);

    /* Read the rotation groups */
    snew(grpbuf, rot->ngrp);
    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        snew(grpbuf[g], STRLEN);
        CTYPE("Rotation group name");
        sprintf(buf, "rot-group%d", g);
        STYPE(buf, grpbuf[g], "");

        CTYPE("Rotation potential. Can be iso, iso-pf, pm, pm-pf, rm, rm-pf, rm2, rm2-pf, flex, flex-t, flex2, flex2-t");
        sprintf(buf, "rot-type%d", g);
        ETYPE(buf, rotg->eType, erotg_names);

        CTYPE("Use mass-weighting of the rotation group positions");
        sprintf(buf, "rot-massw%d", g);
        ETYPE(buf, rotg->bMassW, yesno_names);

        CTYPE("Rotation vector, will get normalized");
        sprintf(buf, "rot-vec%d", g);
        STYPE(buf, s_vec, "1.0 0.0 0.0");
        string2dvec(s_vec, vec);
        /* Normalize the rotation vector */
        if (dnorm(vec) != 0)
        {
            dsvmul(1.0/dnorm(vec), vec, vec);
        }
        else
        {
            sprintf(warn_buf, "rot-vec%d = 0", g);
            warning_error(wi, warn_buf);
        }
        fprintf(stderr, "%s Group %d (%s) normalized rot. vector: %f %f %f\n",
                RotStr, g, erotg_names[rotg->eType], vec[0], vec[1], vec[2]);
        for (m = 0; m < DIM; m++)
        {
            rotg->vec[m] = vec[m];
        }

        CTYPE("Pivot point for the potentials iso, pm, rm, and rm2 (nm)");
        sprintf(buf, "rot-pivot%d", g);
        STYPE(buf, s_vec, "0.0 0.0 0.0");
        clear_dvec(vec);
        if ( (rotg->eType == erotgISO) || (rotg->eType == erotgPM) || (rotg->eType == erotgRM) || (rotg->eType == erotgRM2) )
        {
            string2dvec(s_vec, vec);
        }
        for (m = 0; m < DIM; m++)
        {
            rotg->pivot[m] = vec[m];
        }

        CTYPE("Rotation rate (degree/ps) and force constant (kJ/(mol*nm^2))");
        sprintf(buf, "rot-rate%d", g);
        RTYPE(buf, rotg->rate, 0.0);

        sprintf(buf, "rot-k%d", g);
        RTYPE(buf, rotg->k, 0.0);
        if (rotg->k <= 0.0)
        {
            sprintf(warn_buf, "rot-k%d <= 0", g);
            warning_note(wi, warn_buf);
        }

        CTYPE("Slab distance for flexible axis rotation (nm)");
        sprintf(buf, "rot-slab-dist%d", g);
        RTYPE(buf, rotg->slab_dist, 1.5);
        if (rotg->slab_dist <= 0.0)
        {
            sprintf(warn_buf, "rot-slab-dist%d <= 0", g);
            warning_error(wi, warn_buf);
        }

        CTYPE("Minimum value of Gaussian function for the force to be evaluated (for flex* potentials)");
        sprintf(buf, "rot-min-gauss%d", g);
        RTYPE(buf, rotg->min_gaussian, 1e-3);
        if (rotg->min_gaussian <= 0.0)
        {
            sprintf(warn_buf, "rot-min-gauss%d <= 0", g);
            warning_error(wi, warn_buf);
        }

        CTYPE("Value of additive constant epsilon' (nm^2) for rm2* and flex2* potentials");
        sprintf(buf, "rot-eps%d", g);
        RTYPE(buf, rotg->eps, 1e-4);
        if ( (rotg->eps <= 0.0) && (rotg->eType == erotgRM2 || rotg->eType == erotgFLEX2) )
        {
            sprintf(warn_buf, "rot-eps%d <= 0", g);
            warning_error(wi, warn_buf);
        }

        CTYPE("Fitting method to determine angle of rotation group (rmsd, norm, or potential)");
        sprintf(buf, "rot-fit-method%d", g);
        ETYPE(buf, rotg->eFittype, erotg_fitnames);
        CTYPE("For fit type 'potential', nr. of angles around the reference for which the pot. is evaluated");
        sprintf(buf, "rot-potfit-nsteps%d", g);
        ITYPE(buf, rotg->PotAngle_nstep, 21);
        if ( (rotg->eFittype == erotgFitPOT) && (rotg->PotAngle_nstep < 1) )
        {
            sprintf(warn_buf, "rot-potfit-nsteps%d < 1", g);
            warning_error(wi, warn_buf);
        }
        CTYPE("For fit type 'potential', distance in degrees between two consecutive angles");
        sprintf(buf, "rot-potfit-step%d", g);
        RTYPE(buf, rotg->PotAngle_step, 0.25);
    }

    *ninp_p   = ninp;
    *inp_p    = inp;

    return grpbuf;
}


/* Check whether the box is unchanged */
static void check_box_unchanged(matrix f_box, matrix box, char fn[], warninp_t wi)
{
    int      i, ii;
    gmx_bool bSame = TRUE;
    char     warn_buf[STRLEN];


    for (i = 0; i < DIM; i++)
    {
        for (ii = 0; ii < DIM; ii++)
        {
            if (f_box[i][ii] != box[i][ii])
            {
                bSame = FALSE;
            }
        }
    }
    if (!bSame)
    {
        sprintf(warn_buf, "%s Box size in reference file %s differs from actual box size!",
                RotStr, fn);
        warning(wi, warn_buf);
        pr_rvecs(stderr, 0, "Your box is:", box, 3);
        pr_rvecs(stderr, 0, "Box in file:", f_box, 3);
    }
}


/* Extract the reference positions for the rotation group(s) */
extern void set_reference_positions(
        t_rot *rot, rvec *x, matrix box,
        const char *fn, gmx_bool bSet, warninp_t wi)
{
    int         g, i, ii;
    t_rotgrp   *rotg;
    t_trnheader header;    /* Header information of reference file */
    char        base[STRLEN], extension[STRLEN], reffile[STRLEN];
    char       *extpos;
    rvec        f_box[3];  /* Box from reference file */


    /* Base name and extension of the reference file: */
    strncpy(base, fn, STRLEN - 1);
    base[STRLEN-1] = '\0';
    extpos         = strrchr(base, '.');
    strcpy(extension, extpos+1);
    *extpos = '\0';


    for (g = 0; g < rot->ngrp; g++)
    {
        rotg = &rot->grp[g];
        fprintf(stderr, "%s group %d has %d reference positions.\n", RotStr, g, rotg->nat);
        snew(rotg->x_ref, rotg->nat);

        /* Construct the name for the file containing the reference positions for this group: */
        sprintf(reffile, "%s.%d.%s", base, g, extension);

        /* If the base filename for the reference position files was explicitly set by
         * the user, we issue a fatal error if the group file can not be found */
        if (bSet && !gmx_fexist(reffile))
        {
            gmx_fatal(FARGS, "%s The file containing the reference positions was not found.\n"
                      "Expected the file '%s' for group %d.\n",
                      RotStr, reffile, g);
        }

        if (gmx_fexist(reffile))
        {
            fprintf(stderr, "  Reading them from %s.\n", reffile);
            read_trnheader(reffile, &header);
            if (rotg->nat != header.natoms)
            {
                gmx_fatal(FARGS, "Number of atoms in file %s (%d) does not match the number of atoms in rotation group (%d)!\n",
                          reffile, header.natoms, rotg->nat);
            }
            read_trn(reffile, &header.step, &header.t, &header.lambda, f_box, &header.natoms, rotg->x_ref, NULL, NULL);

            /* Check whether the box is unchanged and output a warning if not: */
            check_box_unchanged(f_box, box, reffile, wi);
        }
        else
        {
            fprintf(stderr, " Saving them to %s.\n", reffile);
            for (i = 0; i < rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], rotg->x_ref[i]);
            }
            write_trn(reffile, g, 0.0, 0.0, box, rotg->nat, rotg->x_ref, NULL, NULL);
        }
    }
}


extern void make_rotation_groups(t_rot *rot, char **rotgnames, t_blocka *grps, char **gnames)
{
    int       g, ig = -1, i;
    t_rotgrp *rotg;


    for (g = 0; g < rot->ngrp; g++)
    {
        rotg      = &rot->grp[g];
        ig        = search_string(rotgnames[g], grps->nr, gnames);
        rotg->nat = grps->index[ig+1] - grps->index[ig];

        if (rotg->nat > 0)
        {
            fprintf(stderr, "Rotation group %d '%s' has %d atoms\n", g, rotgnames[g], rotg->nat);
            snew(rotg->ind, rotg->nat);
            for (i = 0; i < rotg->nat; i++)
            {
                rotg->ind[i] = grps->a[grps->index[ig]+i];
            }
        }
        else
        {
            gmx_fatal(FARGS, "Rotation group %d '%s' is empty", g, rotgnames[g]);
        }
    }
}
