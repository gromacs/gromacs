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
#include "gmxpre.h"

#include <cstdio>

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

static const std::string RotStr = "Enforced rotation:";

static void string2dvec(char buf[], dvec nums)
{
    if (sscanf(buf, "%lf%lf%lf", &nums[0], &nums[1], &nums[2]) != 3)
    {
        gmx_fatal(FARGS, "Expected three numbers at input line %s", buf);
    }
}


extern std::vector<std::string> read_rotparams(std::vector<t_inpfile>* inp, t_rot* rot, WarningHandler* wi)
{
    int       g, m;
    char      buf[STRLEN];
    char      warn_buf[STRLEN];
    dvec      vec;
    t_rotgrp* rotg;

    /* read rotation parameters */
    printStringNoNewline(
            inp,
            "Output frequency for angle, torque and rotation potential energy for the whole group");
    rot->nstrout = get_eint(inp, "rot-nstrout", 100, wi);
    printStringNoNewline(inp,
                         "Output frequency for per-slab data (angles, torques and slab centers)");
    rot->nstsout = get_eint(inp, "rot-nstsout", 1000, wi);
    printStringNoNewline(inp, "Number of rotation groups");
    int numGroups = get_eint(inp, "rot-ngroups", 1, wi);

    if (numGroups < 1)
    {
        gmx_fatal(FARGS, "rot-ngroups should be >= 1");
    }

    rot->grp.resize(numGroups);

    /* Read the rotation groups */
    std::vector<std::string> rotateGroups(numGroups);
    char                     readBuffer[STRLEN];
    char                     s_vec[STRLEN];
    for (g = 0; g < numGroups; g++)
    {
        rotg = &rot->grp[g];
        printStringNoNewline(inp, "Rotation group name");
        sprintf(buf, "rot-group%d", g);
        setStringEntry(inp, buf, readBuffer, "");
        rotateGroups[g] = readBuffer;
        printStringNoNewline(inp,
                             "Rotation potential. Can be iso, iso-pf, pm, pm-pf, rm, rm-pf, rm2, "
                             "rm2-pf, flex, flex-t, flex2, flex2-t");
        sprintf(buf, "rot-type%d", g);
        rotg->eType = getEnum<EnforcedRotationGroupType>(inp, buf, wi);

        printStringNoNewline(inp, "Use mass-weighting of the rotation group positions");
        sprintf(buf, "rot-massw%d", g);
        rotg->bMassW = getEnum<Boolean>(inp, buf, wi) != Boolean::No;

        printStringNoNewline(inp, "Rotation vector, will get normalized");
        sprintf(buf, "rot-vec%d", g);
        setStringEntry(inp, buf, s_vec, "1.0 0.0 0.0");
        string2dvec(s_vec, vec);
        /* Normalize the rotation vector */
        if (dnorm(vec) != 0)
        {
            dsvmul(1.0 / dnorm(vec), vec, vec);
        }
        else
        {
            sprintf(warn_buf, "rot-vec%d = 0", g);
            wi->addError(warn_buf);
        }
        fprintf(stderr,
                "%s Group %d (%s) normalized rot. vector: %f %f %f\n",
                RotStr.c_str(),
                g,
                enumValueToString(rotg->eType),
                vec[0],
                vec[1],
                vec[2]);
        for (m = 0; m < DIM; m++)
        {
            rotg->inputVec[m] = vec[m];
        }

        printStringNoNewline(inp, "Pivot point for the potentials iso, pm, rm, and rm2 (nm)");
        sprintf(buf, "rot-pivot%d", g);
        setStringEntry(inp, buf, s_vec, "0.0 0.0 0.0");
        clear_dvec(vec);
        if ((rotg->eType == EnforcedRotationGroupType::Iso)
            || (rotg->eType == EnforcedRotationGroupType::Pm)
            || (rotg->eType == EnforcedRotationGroupType::Rm)
            || (rotg->eType == EnforcedRotationGroupType::Rm2))
        {
            string2dvec(s_vec, vec);
        }
        for (m = 0; m < DIM; m++)
        {
            rotg->pivot[m] = vec[m];
        }

        printStringNoNewline(inp, "Rotation rate (degree/ps) and force constant (kJ/(mol*nm^2))");
        sprintf(buf, "rot-rate%d", g);
        rotg->rate = get_ereal(inp, buf, 0.0, wi);

        sprintf(buf, "rot-k%d", g);
        rotg->k = get_ereal(inp, buf, 0.0, wi);
        if (rotg->k <= 0.0)
        {
            sprintf(warn_buf, "rot-k%d <= 0", g);
            wi->addNote(warn_buf);
        }

        printStringNoNewline(inp, "Slab distance for flexible axis rotation (nm)");
        sprintf(buf, "rot-slab-dist%d", g);
        rotg->slab_dist = get_ereal(inp, buf, 1.5, wi);
        if (rotg->slab_dist <= 0.0)
        {
            sprintf(warn_buf, "rot-slab-dist%d <= 0", g);
            wi->addError(warn_buf);
        }

        printStringNoNewline(inp,
                             "Minimum value of Gaussian function for the force to be evaluated "
                             "(for flex* potentials)");
        sprintf(buf, "rot-min-gauss%d", g);
        rotg->min_gaussian = get_ereal(inp, buf, 1e-3, wi);
        if (rotg->min_gaussian <= 0.0)
        {
            sprintf(warn_buf, "rot-min-gauss%d <= 0", g);
            wi->addError(warn_buf);
        }

        printStringNoNewline(
                inp, "Value of additive constant epsilon' (nm^2) for rm2* and flex2* potentials");
        sprintf(buf, "rot-eps%d", g);
        rotg->eps = get_ereal(inp, buf, 1e-4, wi);
        if ((rotg->eps <= 0.0)
            && (rotg->eType == EnforcedRotationGroupType::Rm2
                || rotg->eType == EnforcedRotationGroupType::Flex2))
        {
            sprintf(warn_buf, "rot-eps%d <= 0", g);
            wi->addError(warn_buf);
        }

        printStringNoNewline(
                inp,
                "Fitting method to determine angle of rotation group (rmsd, norm, or potential)");
        sprintf(buf, "rot-fit-method%d", g);
        rotg->eFittype = getEnum<RotationGroupFitting>(inp, buf, wi);
        printStringNoNewline(inp,
                             "For fit type 'potential', nr. of angles around the reference for "
                             "which the pot. is evaluated");
        sprintf(buf, "rot-potfit-nsteps%d", g);
        rotg->PotAngle_nstep = get_eint(inp, buf, 21, wi);
        if ((rotg->eFittype == RotationGroupFitting::Pot) && (rotg->PotAngle_nstep < 1))
        {
            sprintf(warn_buf, "rot-potfit-nsteps%d < 1", g);
            wi->addError(warn_buf);
        }
        printStringNoNewline(
                inp, "For fit type 'potential', distance in degrees between two consecutive angles");
        sprintf(buf, "rot-potfit-step%d", g);
        rotg->PotAngle_step = get_ereal(inp, buf, 0.25, wi);
    }

    return rotateGroups;
}


/* Check whether the box is unchanged */
static void check_box_unchanged(matrix f_box, matrix box, const char fn[], WarningHandler* wi)
{
    int  i, ii;
    bool bSame = TRUE;
    char warn_buf[STRLEN];


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
        sprintf(warn_buf, "%s Box size in reference file %s differs from actual box size!", RotStr.c_str(), fn);
        wi->addWarning(warn_buf);
        pr_rvecs(stderr, 0, "Your box is:", box, 3);
        pr_rvecs(stderr, 0, "Box in file:", f_box, 3);
    }
}


/* Extract the reference positions for the rotation group(s) */
extern void set_reference_positions(t_rot* rot, rvec* x, matrix box, const char* fn, bool bSet, WarningHandler* wi)
{
    int              i, ii;
    t_rotgrp*        rotg;
    gmx_trr_header_t header;   /* Header information of reference file */
    rvec             f_box[3]; /* Box from reference file */

    for (int g = 0; g < gmx::ssize(rot->grp); g++)
    {
        rotg = &rot->grp[g];
        fprintf(stderr, "%s group %d has %d reference positions.\n", RotStr.c_str(), g, rotg->nat);
        rotg->x_ref_original.resize(rotg->nat);

        /* Construct the name for the file containing the reference positions for this group: */
        const std::filesystem::path reffile =
                gmx::concatenateBeforeExtension(fn, gmx::formatString(".%d", g));
        const std::string reffileString = reffile.string();

        /* If the base filename for the reference position files was explicitly set by
         * the user, we issue a fatal error if the group file can not be found */
        if (bSet && !gmx_fexist(reffile))
        {
            gmx_fatal(FARGS,
                      "%s The file containing the reference positions was not found.\n"
                      "Expected the file '%s' for group %d.\n",
                      RotStr.c_str(),
                      reffileString.c_str(),
                      g);
        }

        if (gmx_fexist(reffile))
        {
            fprintf(stderr, "  Reading them from %s.\n", reffileString.c_str());
            gmx_trr_read_single_header(reffile, &header);
            if (rotg->nat != header.natoms)
            {
                gmx_fatal(FARGS,
                          "Number of atoms in file %s (%d) does not match the number of atoms in "
                          "rotation group (%d)!\n",
                          reffileString.c_str(),
                          header.natoms,
                          rotg->nat);
            }
            gmx_trr_read_single_frame(reffile,
                                      &header.step,
                                      &header.t,
                                      &header.lambda,
                                      f_box,
                                      &header.natoms,
                                      as_rvec_array(rotg->x_ref_original.data()),
                                      nullptr,
                                      nullptr);

            /* Check whether the box is unchanged and output a warning if not: */
            check_box_unchanged(f_box, box, reffileString.c_str(), wi);
        }
        else
        {
            fprintf(stderr, " Saving them to %s.\n", reffileString.c_str());
            for (i = 0; i < rotg->nat; i++)
            {
                ii = rotg->ind[i];
                copy_rvec(x[ii], rotg->x_ref_original[i]);
            }
            gmx_trr_write_single_frame(
                    reffile, g, 0.0, 0.0, box, rotg->nat, as_rvec_array(rotg->x_ref_original.data()), nullptr, nullptr);
        }
    }
}


extern void make_rotation_groups(t_rot*                           rot,
                                 gmx::ArrayRef<const std::string> rotateGroupNames,
                                 gmx::ArrayRef<const IndexGroup>  indexGroups)
{
    for (int g = 0; g < gmx::ssize(rot->grp); g++)
    {
        t_rotgrp* rotg = &rot->grp[g];
        int       ig   = getGroupIndex(rotateGroupNames[g], indexGroups);
        rotg->nat      = gmx::ssize(indexGroups[ig].particleIndices);

        if (rotg->nat > 0)
        {
            fprintf(stderr, "Rotation group %d '%s' has %d atoms\n", g, rotateGroupNames[g].c_str(), rotg->nat);
            snew(rotg->ind, rotg->nat);
            for (int i = 0; i < rotg->nat; i++)
            {
                rotg->ind[i] = indexGroups[ig].particleIndices[i];
            }
        }
        else
        {
            gmx_fatal(FARGS, "Rotation group %d '%s' is empty", g, rotateGroupNames[g].c_str());
        }
    }
}
