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

#include "genconf.h"

#include <cmath>
#include <cstdio>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/3dtransforms.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

static void
rand_rot(int natoms, rvec x[], rvec v[], vec4 xrot[], vec4 vrot[], gmx::DefaultRandomEngine* rng, const rvec max_rot)
{
    mat4                               mt1, mt2, mr[DIM], mtemp1, mtemp2, mtemp3, mxtot, mvtot;
    rvec                               xcm;
    real                               phi;
    int                                i, m;
    gmx::UniformRealDistribution<real> dist(-1.0, 1.0);

    clear_rvec(xcm);
    for (i = 0; (i < natoms); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            xcm[m] += x[i][m] / natoms; /* get center of mass of one molecule  */
        }
    }
    fprintf(stderr, "center of geometry: %f, %f, %f\n", xcm[0], xcm[1], xcm[2]);

    /* move c.o.ma to origin */
    gmx_mat4_init_translation(-xcm[XX], -xcm[YY], -xcm[ZZ], mt1);
    for (m = 0; (m < DIM); m++)
    {
        phi = M_PI * max_rot[m] * dist(*rng) / 180;
        gmx_mat4_init_rotation(m, phi, mr[m]);
    }
    gmx_mat4_init_translation(xcm[XX], xcm[YY], xcm[ZZ], mt2);

    /* For gmx_mat4_mmul() we need to multiply in the opposite order
     * compared to normal mathematical notation.
     */
    gmx_mat4_mmul(mtemp1, mt1, mr[XX]);
    gmx_mat4_mmul(mtemp2, mr[YY], mr[ZZ]);
    gmx_mat4_mmul(mtemp3, mtemp1, mtemp2);
    gmx_mat4_mmul(mxtot, mtemp3, mt2);
    gmx_mat4_mmul(mvtot, mr[XX], mtemp2);

    for (i = 0; (i < natoms); i++)
    {
        gmx_mat4_transform_point(mxtot, x[i], xrot[i]);
        gmx_mat4_transform_point(mvtot, v[i], vrot[i]);
    }
}

int gmx_genconf(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] multiplies a given coordinate file by simply stacking them",
        "on top of each other, like a small child playing with wooden blocks.",
        "The program makes a grid of [IT]user-defined[it]",
        "proportions ([TT]-nbox[tt]), ",
        "and interspaces the grid point with an extra space [TT]-dist[tt].[PAR]",
        "When option [TT]-rot[tt] is used the program does not check for overlap",
        "between molecules on grid points. It is recommended to make the box in",
        "the input file at least as big as the coordinates + ",
        "van der Waals radius.[PAR]",
        "If the optional trajectory file is given, conformations are not",
        "generated, but read from this file and translated appropriately to",
        "build the grid."

    };
    const char* bugs[] = { "The program should allow for random displacement of lattice points." };

    int               vol;
    rvec *            x, *xx, *v; /* coordinates? */
    real              t;
    vec4 *            xrot, *vrot;
    PbcType           pbcType;
    matrix            box, boxx; /* box length matrix */
    rvec              shift;
    int               natoms; /* number of atoms in one molecule  */
    int               nres;   /* number of molecules? */
    int               i, j, k, l, m, ndx, nrdx, nx, ny, nz;
    t_trxstatus*      status;
    bool              bTRX;
    gmx_output_env_t* oenv;

    t_filenm fnm[] = { { efSTX, "-f", "conf", ffREAD },
                       { efSTO, "-o", "out", ffWRITE },
                       { efTRX, "-trj", nullptr, ffOPTRD } };
#define NFILE asize(fnm)
    rvec     nrbox   = { 1, 1, 1 };
    int      seed    = 0;                 /* seed for random number generator */
    gmx_bool bRandom = FALSE;             /* False: no random rotations */
    gmx_bool bRenum  = TRUE;              /* renumber residues */
    rvec     dist    = { 0, 0, 0 };       /* space added between molecules ? */
    rvec     max_rot = { 180, 180, 180 }; /* maximum rotation */
    t_pargs  pa[]    = {
        { "-nbox", FALSE, etRVEC, { nrbox }, "Number of boxes" },
        { "-dist", FALSE, etRVEC, { dist }, "Distance between boxes" },
        { "-seed", FALSE, etINT, { &seed }, "Random generator seed (0 means generate)" },
        { "-rot", FALSE, etBOOL, { &bRandom }, "Randomly rotate conformations" },
        { "-maxrot", FALSE, etRVEC, { max_rot }, "Maximum random rotation" },
        { "-renumber", FALSE, etBOOL, { &bRenum }, "Renumber residues" }
    };

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
    }
    gmx::DefaultRandomEngine rng(seed);

    bTRX = ftp2bSet(efTRX, NFILE, fnm);
    nx   = gmx::roundToInt(nrbox[XX]);
    ny   = gmx::roundToInt(nrbox[YY]);
    nz   = gmx::roundToInt(nrbox[ZZ]);

    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    {
        gmx_fatal(FARGS, "Number of boxes (-nbox) should be larger than zero");
    }

    vol = nx * ny * nz; /* calculate volume in grid points (= nr. molecules) */

    gmx_mtop_t mtop;
    bool       haveTop = false;
    readConfAndTopology(opt2fn("-f", NFILE, fnm), &haveTop, &mtop, &pbcType, &x, &v, box);
    t_atoms atoms = gmx_mtop_global_atoms(mtop);
    natoms        = atoms.nr;
    nres          = atoms.nres; /* nr of residues in one element? */
    /* make space for all the atoms */
    add_t_atoms(&atoms, natoms * (vol - 1), nres * (vol - 1));
    srenew(x, natoms * vol); /* get space for coordinates of all atoms */
    srenew(v, natoms * vol); /* velocities. not really needed? */
    snew(xrot, natoms);      /* get space for rotation matrix? */
    snew(vrot, natoms);

    if (bTRX)
    {
        if (read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &xx, boxx) == 0)
        {
            gmx_fatal(FARGS, "No atoms in trajectory %s", ftp2fn(efTRX, NFILE, fnm));
        }
    }
    else
    {
        snew(xx, natoms);
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(x[i], xx[i]);
        }
    }


    for (k = 0; (k < nz); k++) /* loop over all gridpositions    */
    {
        shift[ZZ] = k * (dist[ZZ] + box[ZZ][ZZ]);

        for (j = 0; (j < ny); j++)
        {
            shift[YY] = j * (dist[YY] + box[YY][YY]) + k * box[ZZ][YY];

            for (i = 0; (i < nx); i++)
            {
                shift[XX] = i * (dist[XX] + box[XX][XX]) + j * box[YY][XX] + k * box[ZZ][XX];

                ndx  = (i * ny * nz + j * nz + k) * natoms;
                nrdx = (i * ny * nz + j * nz + k) * nres;

                /* Random rotation on input coords */
                if (bRandom)
                {
                    rand_rot(natoms, xx, v, xrot, vrot, &rng, max_rot);
                }

                for (l = 0; (l < natoms); l++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        if (bRandom)
                        {
                            x[ndx + l][m] = xrot[l][m];
                            v[ndx + l][m] = vrot[l][m];
                        }
                        else
                        {
                            x[ndx + l][m] = xx[l][m];
                            v[ndx + l][m] = v[l][m];
                        }
                    }
                    if (pbcType == PbcType::Screw && i % 2 == 1)
                    {
                        /* Rotate around x axis */
                        for (m = YY; m <= ZZ; m++)
                        {
                            x[ndx + l][m] = box[YY][m] + box[ZZ][m] - x[ndx + l][m];
                            v[ndx + l][m] = -v[ndx + l][m];
                        }
                    }
                    for (m = 0; (m < DIM); m++)
                    {
                        x[ndx + l][m] += shift[m];
                    }
                    atoms.atom[ndx + l].resind = nrdx + atoms.atom[l].resind;
                    atoms.atomname[ndx + l]    = atoms.atomname[l];
                }

                for (l = 0; (l < nres); l++)
                {
                    atoms.resinfo[nrdx + l] = atoms.resinfo[l];
                    if (bRenum)
                    {
                        atoms.resinfo[nrdx + l].nr += nrdx;
                    }
                }
                if (bTRX)
                {
                    if (!read_next_x(oenv, status, &t, xx, boxx) && ((i + 1) * (j + 1) * (k + 1) < vol))
                    {
                        gmx_fatal(FARGS, "Not enough frames in trajectory");
                    }
                }
            }
        }
    }
    if (bTRX)
    {
        close_trx(status);
    }

    /* make box bigger */
    for (m = 0; (m < DIM); m++)
    {
        box[m][m] += dist[m];
    }
    svmul(nx, box[XX], box[XX]);
    svmul(ny, box[YY], box[YY]);
    svmul(nz, box[ZZ], box[ZZ]);
    if (pbcType == PbcType::Screw && nx % 2 == 0)
    {
        /* With an even number of boxes in x we can forgot about the screw */
        pbcType = PbcType::Xyz;
    }

    /*depending on how you look at it, this is either a nasty hack or the way it should work*/
    if (bRenum)
    {
        for (i = 0; i < atoms.nres; i++)
        {
            atoms.resinfo[i].nr = i + 1;
        }
    }

    write_sto_conf(opt2fn("-o", NFILE, fnm), *mtop.name, &atoms, x, v, pbcType, box);

    sfree(x);
    sfree(v);
    sfree(xrot);
    sfree(vrot);
    sfree(xx);
    done_atom(&atoms);
    output_env_done(oenv);

    return 0;
}
