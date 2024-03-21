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

#include "editconf.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/princ.h"
#include "gromacs/gmxlib/conformation_utilities.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

static real calc_mass(t_atoms* atoms, gmx_bool bGetMass, AtomProperties* aps)
{
    real tmass;
    int  i;

    tmass = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        if (bGetMass)
        {
            aps->setAtomProperty(epropMass,
                                 std::string(*atoms->resinfo[atoms->atom[i].resind].name),
                                 std::string(*atoms->atomname[i]),
                                 &(atoms->atom[i].m));
        }
        tmass += atoms->atom[i].m;
    }

    return tmass;
}

static real calc_geom(int isize, const int* index, rvec* x, rvec geom_center, rvec minval, rvec maxval, gmx_bool bDiam)
{
    real diam2, d;
    int  ii, i, j;

    clear_rvec(geom_center);
    diam2 = 0;
    if (isize == 0)
    {
        clear_rvec(minval);
        clear_rvec(maxval);
    }
    else
    {
        if (index)
        {
            ii = index[0];
        }
        else
        {
            ii = 0;
        }
        for (j = 0; j < DIM; j++)
        {
            minval[j] = maxval[j] = x[ii][j];
        }
        for (i = 0; i < isize; i++)
        {
            if (index)
            {
                ii = index[i];
            }
            else
            {
                ii = i;
            }
            rvec_inc(geom_center, x[ii]);
            for (j = 0; j < DIM; j++)
            {
                if (x[ii][j] < minval[j])
                {
                    minval[j] = x[ii][j];
                }
                if (x[ii][j] > maxval[j])
                {
                    maxval[j] = x[ii][j];
                }
            }
            if (bDiam)
            {
                if (index)
                {
                    for (j = i + 1; j < isize; j++)
                    {
                        d     = distance2(x[ii], x[index[j]]);
                        diam2 = std::max(d, diam2);
                    }
                }
                else
                {
                    for (j = i + 1; j < isize; j++)
                    {
                        d     = distance2(x[i], x[j]);
                        diam2 = std::max(d, diam2);
                    }
                }
            }
        }
        svmul(1.0 / isize, geom_center, geom_center);
    }

    return std::sqrt(diam2);
}

static void center_conf(int natom, rvec* x, rvec center, rvec geom_cent)
{
    int  i;
    rvec shift;

    rvec_sub(center, geom_cent, shift);

    printf("    shift       :%7.3f%7.3f%7.3f (nm)\n", shift[XX], shift[YY], shift[ZZ]);

    for (i = 0; (i < natom); i++)
    {
        rvec_inc(x[i], shift);
    }
}

static void scale_conf(int natom, rvec x[], matrix box, const rvec scale)
{
    int i, j;

    for (i = 0; i < natom; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            x[i][j] *= scale[j];
        }
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            box[i][j] *= scale[j];
        }
    }
}

static void read_bfac(const char* fn, int* n_bfac, double** bfac_val, int** bfac_nr)
{
    int    i;
    char** bfac_lines;

    *n_bfac = get_lines(fn, &bfac_lines);
    snew(*bfac_val, *n_bfac);
    snew(*bfac_nr, *n_bfac);
    fprintf(stderr, "Reading %d B-factors from %s\n", *n_bfac, fn);
    for (i = 0; (i < *n_bfac); i++)
    {
        sscanf(bfac_lines[i], "%d %lf", &(*bfac_nr)[i], &(*bfac_val)[i]);
    }
}

static void set_pdb_conf_bfac(int natoms, int nres, t_atoms* atoms, int n_bfac, double* bfac, int* bfac_nr, gmx_bool peratom)
{
    real     bfac_min, bfac_max;
    int      i, n;
    gmx_bool found;

    if (n_bfac > atoms->nres)
    {
        peratom = TRUE;
    }

    bfac_max = -1e10;
    bfac_min = 1e10;
    for (i = 0; (i < n_bfac); i++)
    {
        /*    if ((bfac_nr[i]-1<0) || (bfac_nr[i]-1>=atoms->nr))
           gmx_fatal(FARGS,"Index of B-Factor %d is out of range: %d (%g)",
           i+1,bfac_nr[i],bfac[i]); */
        if (bfac[i] > bfac_max)
        {
            bfac_max = bfac[i];
        }
        if (bfac[i] < bfac_min)
        {
            bfac_min = bfac[i];
        }
    }
    while ((bfac_max > 99.99) || (bfac_min < -99.99))
    {
        fprintf(stderr,
                "Range of values for B-factors too large (min %g, max %g) "
                "will scale down a factor 10\n",
                bfac_min,
                bfac_max);
        for (i = 0; (i < n_bfac); i++)
        {
            bfac[i] /= 10;
        }
        bfac_max /= 10;
        bfac_min /= 10;
    }
    while ((std::abs(bfac_max) < 0.5) && (std::abs(bfac_min) < 0.5))
    {
        fprintf(stderr,
                "Range of values for B-factors too small (min %g, max %g) "
                "will scale up a factor 10\n",
                bfac_min,
                bfac_max);
        for (i = 0; (i < n_bfac); i++)
        {
            bfac[i] *= 10;
        }
        bfac_max *= 10;
        bfac_min *= 10;
    }

    for (i = 0; (i < natoms); i++)
    {
        atoms->pdbinfo[i].bfac = 0;
    }

    if (!peratom)
    {
        fprintf(stderr, "Will attach %d B-factors to %d residues\n", n_bfac, nres);
        for (i = 0; (i < n_bfac); i++)
        {
            found = FALSE;
            for (n = 0; (n < natoms); n++)
            {
                if (bfac_nr[i] == atoms->resinfo[atoms->atom[n].resind].nr)
                {
                    atoms->pdbinfo[n].bfac = bfac[i];
                    found                  = TRUE;
                }
            }
            if (!found)
            {
                gmx_warning("Residue nr %d not found\n", bfac_nr[i]);
            }
        }
    }
    else
    {
        fprintf(stderr, "Will attach %d B-factors to %d atoms\n", n_bfac, natoms);
        for (i = 0; (i < n_bfac); i++)
        {
            atoms->pdbinfo[bfac_nr[i] - 1].bfac = bfac[i];
        }
    }
}

static void pdb_legend(FILE* out, int natoms, int nres, t_atoms* atoms, rvec x[])
{
    real bfac_min, bfac_max, xmin, ymin, zmin;
    int  i;
    int  space = ' ';

    bfac_max = -1e10;
    bfac_min = 1e10;
    xmin     = 1e10;
    ymin     = 1e10;
    zmin     = 1e10;
    for (i = 0; (i < natoms); i++)
    {
        xmin     = std::min(xmin, x[i][XX]);
        ymin     = std::min(ymin, x[i][YY]);
        zmin     = std::min(zmin, x[i][ZZ]);
        bfac_min = std::min(bfac_min, atoms->pdbinfo[i].bfac);
        bfac_max = std::max(bfac_max, atoms->pdbinfo[i].bfac);
    }
    fprintf(stderr, "B-factors range from %g to %g\n", bfac_min, bfac_max);
    for (i = 1; (i < 12); i++)
    {
        fprintf(out,
                "%-6s%5d  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                "ATOM  ",
                natoms + 1 + i,
                "CA",
                "LEG",
                space,
                nres + 1,
                space,
                (xmin + (i * 0.12)) * 10,
                ymin * 10,
                zmin * 10,
                1.0,
                bfac_min + ((i - 1.0) * (bfac_max - bfac_min) / 10));
    }
}

static void visualize_images(const char* fn, PbcType pbcType, matrix box)
{
    t_atoms atoms;
    rvec*   img;
    char *  c, *ala;
    int     nat, i;

    nat = NTRICIMG + 1;
    init_t_atoms(&atoms, nat, FALSE);
    atoms.nr = nat;
    snew(img, nat);
    /* FIXME: Constness should not be cast away */
    c   = const_cast<char*>("C");
    ala = const_cast<char*>("ALA");
    for (i = 0; i < nat; i++)
    {
        atoms.atomname[i]        = &c;
        atoms.atom[i].resind     = i;
        atoms.resinfo[i].name    = &ala;
        atoms.resinfo[i].nr      = i + 1;
        atoms.resinfo[i].chainid = 'A' + i / NCUCVERT;
    }
    calc_triclinic_images(box, img + 1);

    write_sto_conf(fn, "Images", &atoms, img, nullptr, pbcType, box);

    done_atom(&atoms);
    sfree(img);
}

static void visualize_box(FILE* out, int a0, int r0, matrix box, const rvec gridsize)
{
    int*  edge;
    rvec *vert, shift;
    int   nx, ny, nz, nbox, nat;
    int   i, j, x, y, z;
    int   rectedge[24] = { 0, 1, 1, 3, 3, 2, 0, 2, 0, 4, 1, 5, 3, 7, 2, 6, 4, 5, 5, 7, 7, 6, 6, 4 };

    a0++;
    r0++;

    nx   = gmx::roundToInt(gridsize[XX]);
    ny   = gmx::roundToInt(gridsize[YY]);
    nz   = gmx::roundToInt(gridsize[ZZ]);
    nbox = nx * ny * nz;
    if (TRICLINIC(box))
    {
        nat = nbox * NCUCVERT;
        snew(vert, nat);
        calc_compact_unitcell_vertices(ecenterDEF, box, vert);
        j = 0;
        for (z = 0; z < nz; z++)
        {
            for (y = 0; y < ny; y++)
            {
                for (x = 0; x < nx; x++)
                {
                    for (i = 0; i < DIM; i++)
                    {
                        shift[i] = x * box[0][i] + y * box[1][i] + z * box[2][i];
                    }
                    for (i = 0; i < NCUCVERT; i++)
                    {
                        rvec_add(vert[i], shift, vert[j]);
                        j++;
                    }
                }
            }
        }

        for (i = 0; i < nat; i++)
        {
            gmx_fprintf_pdb_atomline(out,
                                     PdbRecordType::Atom,
                                     a0 + i,
                                     "C",
                                     ' ',
                                     "BOX",
                                     'K' + i / NCUCVERT,
                                     r0 + i,
                                     ' ',
                                     10 * vert[i][XX],
                                     10 * vert[i][YY],
                                     10 * vert[i][ZZ],
                                     1.0,
                                     0.0,
                                     "");
        }

        edge = compact_unitcell_edges();
        for (j = 0; j < nbox; j++)
        {
            for (i = 0; i < NCUCEDGE; i++)
            {
                fprintf(out,
                        "CONECT%5d%5d\n",
                        a0 + j * NCUCVERT + edge[2 * i],
                        a0 + j * NCUCVERT + edge[2 * i + 1]);
            }
        }

        sfree(vert);
    }
    else
    {
        i = 0;
        for (z = 0; z <= 1; z++)
        {
            for (y = 0; y <= 1; y++)
            {
                for (x = 0; x <= 1; x++)
                {
                    gmx_fprintf_pdb_atomline(out,
                                             PdbRecordType::Atom,
                                             a0 + i,
                                             "C",
                                             ' ',
                                             "BOX",
                                             'K' + i / 8,
                                             r0 + i,
                                             ' ',
                                             x * 10 * box[XX][XX],
                                             y * 10 * box[YY][YY],
                                             z * 10 * box[ZZ][ZZ],
                                             1.0,
                                             0.0,
                                             "");
                    i++;
                }
            }
        }
        for (i = 0; i < 24; i += 2)
        {
            fprintf(out, "CONECT%5d%5d\n", a0 + rectedge[i], a0 + rectedge[i + 1]);
        }
    }
}

static void calc_rotmatrix(rvec principal_axis, rvec targetvec, matrix rotmatrix)
{
    rvec rotvec;
    real ux, uy, uz, costheta, sintheta;

    costheta = cos_angle(principal_axis, targetvec);
    sintheta = std::sqrt(1.0 - costheta * costheta); /* sign is always positive since 0<theta<pi */

    /* Determine rotation from cross product with target vector */
    cprod(principal_axis, targetvec, rotvec);
    unitv(rotvec, rotvec);
    printf("Aligning %g %g %g to %g %g %g : xprod  %g %g %g\n",
           principal_axis[XX],
           principal_axis[YY],
           principal_axis[ZZ],
           targetvec[XX],
           targetvec[YY],
           targetvec[ZZ],
           rotvec[XX],
           rotvec[YY],
           rotvec[ZZ]);

    ux              = rotvec[XX];
    uy              = rotvec[YY];
    uz              = rotvec[ZZ];
    rotmatrix[0][0] = ux * ux + (1.0 - ux * ux) * costheta;
    rotmatrix[0][1] = ux * uy * (1 - costheta) - uz * sintheta;
    rotmatrix[0][2] = ux * uz * (1 - costheta) + uy * sintheta;
    rotmatrix[1][0] = ux * uy * (1 - costheta) + uz * sintheta;
    rotmatrix[1][1] = uy * uy + (1.0 - uy * uy) * costheta;
    rotmatrix[1][2] = uy * uz * (1 - costheta) - ux * sintheta;
    rotmatrix[2][0] = ux * uz * (1 - costheta) - uy * sintheta;
    rotmatrix[2][1] = uy * uz * (1 - costheta) + ux * sintheta;
    rotmatrix[2][2] = uz * uz + (1.0 - uz * uz) * costheta;

    printf("Rotation matrix: \n%g %g %g\n%g %g %g\n%g %g %g\n",
           rotmatrix[0][0],
           rotmatrix[0][1],
           rotmatrix[0][2],
           rotmatrix[1][0],
           rotmatrix[1][1],
           rotmatrix[1][2],
           rotmatrix[2][0],
           rotmatrix[2][1],
           rotmatrix[2][2]);
}

static void renum_resnr(t_atoms* atoms, int isize, const int* index, int resnr_start)
{
    int i, resind_prev, resind;

    resind_prev = -1;
    for (i = 0; i < isize; i++)
    {
        resind = atoms->atom[index == nullptr ? i : index[i]].resind;
        if (resind != resind_prev)
        {
            atoms->resinfo[resind].nr = resnr_start;
            resnr_start++;
        }
        resind_prev = resind;
    }
}

int gmx_editconf(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] converts generic structure format to [REF].gro[ref], [TT].g96[tt]",
        "or [REF].pdb[ref].",
        "[PAR]",
        "The box can be modified with options [TT]-box[tt], [TT]-d[tt] and",
        "[TT]-angles[tt]. Both [TT]-box[tt] and [TT]-d[tt]",
        "will center the system in the box, unless [TT]-noc[tt] is used.",
        "The [TT]-center[tt] option can be used to shift the geometric center",
        "of the system from the default of (x/2, y/2, z/2) implied by [TT]-c[tt]",
        "to some other value.",
        "[PAR]",
        "Option [TT]-bt[tt] determines the box type: [TT]triclinic[tt] is a",
        "triclinic box, [TT]cubic[tt] is a rectangular box with all sides equal",
        "[TT]dodecahedron[tt] represents a rhombic dodecahedron and",
        "[TT]octahedron[tt] is a truncated octahedron.",
        "The last two are special cases of a triclinic box.",
        "The length of the three box vectors of the truncated octahedron is the",
        "shortest distance between two opposite hexagons.",
        "Relative to a cubic box with some periodic image distance, the volume of a ",
        "dodecahedron with this same periodic distance is 0.71 times that of the cube, ",
        "and that of a truncated octahedron is 0.77 times.",
        "[PAR]",
        "Option [TT]-box[tt] requires only",
        "one value for a cubic, rhombic dodecahedral, or truncated octahedral box.",
        "[PAR]",
        "With [TT]-d[tt] and a [TT]triclinic[tt] box the size of the system in the [IT]x[it]-, ",
        "[IT]y[it]-,",
        "and [IT]z[it]-directions is used. With [TT]-d[tt] and [TT]cubic[tt],",
        "[TT]dodecahedron[tt] or [TT]octahedron[tt] boxes, the dimensions are set",
        "to the diameter of the system (largest distance between atoms) plus twice",
        "the specified distance.",
        "[PAR]",
        "Option [TT]-angles[tt] is only meaningful with option [TT]-box[tt] and",
        "a triclinic box and cannot be used with option [TT]-d[tt].",
        "[PAR]",
        "When [TT]-n[tt] or [TT]-ndef[tt] is set, a group",
        "can be selected for calculating the size and the geometric center,",
        "otherwise the whole system is used.",
        "[PAR]",
        "[TT]-rotate[tt] rotates the coordinates and velocities.",
        "[PAR]",
        "[TT]-princ[tt] aligns the principal axes of the system along the",
        "coordinate axes, with the longest axis aligned with the [IT]x[it]-axis. ",
        "This may allow you to decrease the box volume,",
        "but beware that molecules can rotate significantly in a nanosecond.",
        "[PAR]",
        "Scaling is applied before any of the other operations are",
        "performed. Boxes and coordinates can be scaled to give a certain density (option",
        "[TT]-density[tt]). Note that this may be inaccurate in case a [REF].gro[ref]",
        "file is given as input. A special feature of the scaling option is that when the",
        "factor -1 is given in one dimension, one obtains a mirror image,",
        "mirrored in one of the planes. When one uses -1 in three dimensions, ",
        "a point-mirror image is obtained.[PAR]",
        "Groups are selected after all operations have been applied.[PAR]",
        "Periodicity can be removed in a crude manner.",
        "It is important that the box vectors at the bottom of your input file",
        "are correct when the periodicity is to be removed.",
        "[PAR]",
        "When writing [REF].pdb[ref] files, B-factors can be",
        "added with the [TT]-bf[tt] option. B-factors are read",
        "from a file with with following format: first line states number of",
        "entries in the file, next lines state an index",
        "followed by a B-factor. The B-factors will be attached per residue",
        "unless the number of B-factors is larger than the number of the residues or unless the",
        "[TT]-atom[tt] option is set. Obviously, any type of numeric data can",
        "be added instead of B-factors. [TT]-legend[tt] will produce",
        "a row of CA atoms with B-factors ranging from the minimum to the",
        "maximum value found, effectively making a legend for viewing.",
        "[PAR]",
        "With the option [TT]-mead[tt] a special [REF].pdb[ref] ([REF].pqr[ref])",
        "file for the MEAD electrostatics",
        "program (Poisson-Boltzmann solver) can be made. A further prerequisite",
        "is that the input file is a run input file.",
        "The B-factor field is then filled with the Van der Waals radius",
        "of the atoms while the occupancy field will hold the charge.",
        "[PAR]",
        "The option [TT]-grasp[tt] is similar, but it puts the charges in the B-factor",
        "and the radius in the occupancy.",
        "[PAR]",
        "Option [TT]-align[tt] allows alignment",
        "of the principal axis of a specified group against the given vector, ",
        "with an optional center of rotation specified by [TT]-aligncenter[tt].",
        "[PAR]",
        "Finally, with option [TT]-label[tt], [TT]editconf[tt] can add a chain identifier",
        "to a [REF].pdb[ref] file, which can be useful for analysis with e.g. Rasmol.",
        "[PAR]",
        "To convert a truncated octrahedron file produced by a package which uses",
        "a cubic box with the corners cut off (such as GROMOS), use::",
        "",
        "  gmx editconf -f in -rotate 0 45 35.264 -bt o -box veclen -o out",
        "",
        "where [TT]veclen[tt] is the size of the cubic box times [SQRT]3[sqrt]/2."
    };
    const char* bugs[] = {
        "For complex molecules, the periodicity removal routine may break down, ",
        "in that case you can use [gmx-trjconv]."
    };
    static real     dist = 0.0;
    static gmx_bool bNDEF = FALSE, bRMPBC = FALSE, bCenter = FALSE, bReadVDW = FALSE, bCONECT = FALSE;
    static gmx_bool peratom = FALSE, bLegend = FALSE, bOrient = FALSE, bMead = FALSE,
                    bGrasp = FALSE, bSig56 = FALSE;
    static rvec scale = { 1, 1, 1 }, newbox = { 0, 0, 0 }, newang = { 90, 90, 90 };
    static real rho = 1000.0, rvdw = 0.12;
    static rvec center = { 0, 0, 0 }, translation = { 0, 0, 0 }, rotangles = { 0, 0, 0 },
                aligncenter = { 0, 0, 0 }, targetvec = { 0, 0, 0 };
    static const char *btype[] = { nullptr,        "triclinic",  "cubic",
                                   "dodecahedron", "octahedron", nullptr },
                      *label   = "A";
    static rvec visbox         = { 0, 0, 0 };
    static int  resnr_start    = -1;
    t_pargs     pa[]           = {
        { "-ndef", FALSE, etBOOL, { &bNDEF }, "Choose output from default index groups" },
        { "-visbox",
          FALSE,
          etRVEC,
          { visbox },
          "HIDDENVisualize a grid of boxes, -1 visualizes the 14 box images" },
        { "-bt", FALSE, etENUM, { btype }, "Box type for [TT]-box[tt] and [TT]-d[tt]" },
        { "-box", FALSE, etRVEC, { newbox }, "Box vector lengths (a,b,c)" },
        { "-angles", FALSE, etRVEC, { newang }, "Angles between the box vectors (bc,ac,ab)" },
        { "-d", FALSE, etREAL, { &dist }, "Distance between the solute and the box" },
        { "-c",
          FALSE,
          etBOOL,
          { &bCenter },
          "Center molecule in box (implied by [TT]-box[tt] and [TT]-d[tt])" },
        { "-center", FALSE, etRVEC, { center }, "Shift the geometrical center to (x,y,z)" },
        { "-aligncenter", FALSE, etRVEC, { aligncenter }, "Center of rotation for alignment" },
        { "-align", FALSE, etRVEC, { targetvec }, "Align to target vector" },
        { "-translate", FALSE, etRVEC, { translation }, "Translation" },
        { "-rotate",
          FALSE,
          etRVEC,
          { rotangles },
          "Rotation around the X, Y and Z axes in degrees" },
        { "-princ", FALSE, etBOOL, { &bOrient }, "Orient molecule(s) along their principal axes" },
        { "-scale", FALSE, etRVEC, { scale }, "Scaling factor" },
        { "-density",
          FALSE,
          etREAL,
          { &rho },
          "Density (g/L) of the output box achieved by scaling" },
        { "-pbc", FALSE, etBOOL, { &bRMPBC }, "Remove the periodicity (make molecule whole again)" },
        { "-resnr", FALSE, etINT, { &resnr_start }, " Renumber residues starting from resnr" },
        { "-grasp",
          FALSE,
          etBOOL,
          { &bGrasp },
          "Store the charge of the atom in the B-factor field and the radius of the atom in the "
          "occupancy field" },
        { "-rvdw",
          FALSE,
          etREAL,
          { &rvdw },
          "Default Van der Waals radius (in nm) if one can not be found in the database or if no "
          "parameters are present in the topology file" },
        { "-sig56",
          FALSE,
          etBOOL,
          { &bSig56 },
          "Use rmin/2 (minimum in the Van der Waals potential) rather than [GRK]sigma[grk]/2 " },
        { "-vdwread",
          FALSE,
          etBOOL,
          { &bReadVDW },
          "Read the Van der Waals radii from the file [TT]vdwradii.dat[tt] rather than computing "
          "the radii based on the force field" },
        { "-atom", FALSE, etBOOL, { &peratom }, "Force B-factor attachment per atom" },
        { "-legend", FALSE, etBOOL, { &bLegend }, "Make B-factor legend" },
        { "-label", FALSE, etSTR, { &label }, "Add chain label for all residues" },
        { "-conect",
          FALSE,
          etBOOL,
          { &bCONECT },
          "Add CONECT records to a [REF].pdb[ref] file when written. Can only be done when a "
          "topology (tpr file) is present" }
    };
#define NPA asize(pa)

    FILE*             out;
    const char *      infile, *outfile;
    int               outftp, inftp, natom, i, j, n_bfac, itype, ntype;
    double *          bfac    = nullptr, c6, c12;
    int*              bfac_nr = nullptr;
    t_topology*       top     = nullptr;
    char *            grpname, *sgrpname, *agrpname;
    int               isize, ssize, numAlignmentAtoms;
    int *             index, *sindex, *aindex;
    rvec *            x, *v, gc, rmin, rmax, size;
    PbcType           pbcType;
    matrix            box, rotmatrix, trans;
    rvec              princd, tmpvec;
    gmx_bool          bIndex, bSetSize, bSetAng, bDist, bSetCenter, bAlign;
    gmx_bool          bHaveV, bScale, bRho, bTranslate, bRotate, bCalcGeom, bCalcDiam;
    real              diam = 0, mass = 0, d, vdw;
    gmx_conect        conect;
    gmx_output_env_t* oenv;
    t_filenm          fnm[] = { { efSTX, "-f", nullptr, ffREAD },
                       { efNDX, "-n", nullptr, ffOPTRD },
                       { efSTO, nullptr, nullptr, ffOPTWR },
                       { efPQR, "-mead", "mead", ffOPTWR },
                       { efDAT, "-bf", "bfact", ffOPTRD } };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, NPA, pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    fprintf(stdout,
            "Note that major changes are planned in future for "
            "editconf, to improve usability and utility.\n");

    bIndex     = opt2bSet("-n", NFILE, fnm) || bNDEF;
    bMead      = opt2bSet("-mead", NFILE, fnm);
    bSetSize   = opt2parg_bSet("-box", NPA, pa);
    bSetAng    = opt2parg_bSet("-angles", NPA, pa);
    bSetCenter = opt2parg_bSet("-center", NPA, pa);
    bDist      = opt2parg_bSet("-d", NPA, pa);
    bAlign     = opt2parg_bSet("-align", NPA, pa);
    /* Only automatically turn on centering without -noc */
    if ((bDist || bSetSize || bSetCenter) && !opt2parg_bSet("-c", NPA, pa))
    {
        bCenter = TRUE;
    }
    bScale     = opt2parg_bSet("-scale", NPA, pa);
    bRho       = opt2parg_bSet("-density", NPA, pa);
    bTranslate = opt2parg_bSet("-translate", NPA, pa);
    bRotate    = opt2parg_bSet("-rotate", NPA, pa);
    if (bScale && bRho)
    {
        fprintf(stderr, "WARNING: setting -density overrides -scale\n");
    }
    bScale    = bScale || bRho;
    bCalcGeom = bCenter || bRotate || bOrient || bScale || bDist;

    GMX_RELEASE_ASSERT(btype[0] != nullptr, "Option setting inconsistency; btype[0] is NULL");

    bCalcDiam = (btype[0][0] == 'c' || btype[0][0] == 'd' || btype[0][0] == 'o');

    infile = ftp2fn(efSTX, NFILE, fnm);
    if (bMead)
    {
        outfile = ftp2fn(efPQR, NFILE, fnm);
    }
    else
    {
        outfile = ftp2fn(efSTO, NFILE, fnm);
    }
    outftp = fn2ftp(outfile);
    inftp  = fn2ftp(infile);

    AtomProperties aps;

    if (bMead && bGrasp)
    {
        printf("Incompatible options -mead and -grasp. Turning off -grasp\n");
        bGrasp = FALSE;
    }
    if ((bGrasp || bCONECT) && (outftp != efPDB))
    {
        gmx_fatal(FARGS,
                  "Output file should be a .pdb file"
                  " when using the -grasp or -conect options\n");
    }
    if ((bMead || bGrasp || bCONECT) && inftp != efTPR)
    {
        gmx_fatal(FARGS,
                  "Input file should be a .tpr file"
                  " when using the -mead, -grasp, or -conect options\n");
    }

    t_symtab symtab;
    char*    name;
    t_atoms  atoms;
    open_symtab(&symtab);
    readConfAndAtoms(infile, &symtab, &name, &atoms, &pbcType, &x, &v, box);
    natom = atoms.nr;
    if (atoms.pdbinfo == nullptr)
    {
        snew(atoms.pdbinfo, atoms.nr);
    }
    atoms.havePdbInfo = TRUE;

    if (fn2ftp(infile) == efPDB)
    {
        get_pdb_atomnumber(&atoms, &aps);
    }
    printf("Read %d atoms\n", atoms.nr);

    /* Get the element numbers if available in a pdb file */
    if (fn2ftp(infile) == efPDB)
    {
        get_pdb_atomnumber(&atoms, &aps);
    }

    if (pbcType != PbcType::No)
    {
        real vol = det(box);
        printf("Volume: %g nm^3, corresponds to roughly %d electrons\n",
               vol,
               100 * (static_cast<int>(vol * 4.5)));
    }

    if (bMead || bGrasp || bCONECT)
    {
        top = read_top(infile, nullptr);
    }

    if (bMead || bGrasp)
    {
        if (atoms.nr != top->atoms.nr)
        {
            gmx_fatal(FARGS, "Atom numbers don't match (%d vs. %d)", atoms.nr, top->atoms.nr);
        }
        snew(atoms.pdbinfo, top->atoms.nr);
        ntype = top->idef.atnr;
        for (i = 0; (i < atoms.nr); i++)
        {
            /* Determine the Van der Waals radius from the force field */
            if (bReadVDW)
            {
                if (!aps.setAtomProperty(epropVDW,
                                         *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                         *top->atoms.atomname[i],
                                         &vdw))
                {
                    vdw = rvdw;
                }
            }
            else
            {
                itype = top->atoms.atom[i].type;
                c12   = top->idef.iparams[itype * ntype + itype].lj.c12;
                c6    = top->idef.iparams[itype * ntype + itype].lj.c6;
                if ((c6 != 0) && (c12 != 0))
                {
                    real sig6;
                    if (bSig56)
                    {
                        sig6 = 2 * c12 / c6;
                    }
                    else
                    {
                        sig6 = c12 / c6;
                    }
                    vdw = 0.5 * gmx::sixthroot(sig6);
                }
                else
                {
                    vdw = rvdw;
                }
            }
            /* Factor of 10 for nm -> Angstroms */
            vdw *= 10;

            if (bMead)
            {
                atoms.pdbinfo[i].occup = top->atoms.atom[i].q;
                atoms.pdbinfo[i].bfac  = vdw;
            }
            else
            {
                atoms.pdbinfo[i].occup = vdw;
                atoms.pdbinfo[i].bfac  = top->atoms.atom[i].q;
            }
        }
    }
    bHaveV = FALSE;
    for (i = 0; (i < natom) && !bHaveV; i++)
    {
        for (j = 0; (j < DIM) && !bHaveV; j++)
        {
            bHaveV = bHaveV || (v[i][j] != 0);
        }
    }
    printf("%selocities found\n", bHaveV ? "V" : "No v");

    if (visbox[0] > 0)
    {
        if (bIndex)
        {
            gmx_fatal(FARGS, "Sorry, can not visualize box with index groups");
        }
        if (outftp != efPDB)
        {
            gmx_fatal(FARGS, "Sorry, can only visualize box with a pdb file");
        }
    }
    else if (visbox[0] == -1)
    {
        visualize_images("images.pdb", pbcType, box);
    }

    /* remove pbc */
    if (bRMPBC)
    {
        rm_gropbc(&atoms, x, box);
    }

    if (bCalcGeom)
    {
        if (bIndex)
        {
            fprintf(stderr, "\nSelect a group for determining the system size:\n");
            get_index(&atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &ssize, &sindex, &sgrpname);
        }
        else
        {
            ssize  = atoms.nr;
            sindex = nullptr;
        }
        diam = calc_geom(ssize, sindex, x, gc, rmin, rmax, bCalcDiam);
        rvec_sub(rmax, rmin, size);
        printf("    system size :%7.3f%7.3f%7.3f (nm)\n", size[XX], size[YY], size[ZZ]);
        if (bCalcDiam)
        {
            printf("    diameter    :%7.3f               (nm)\n", diam);
        }
        printf("    center      :%7.3f%7.3f%7.3f (nm)\n", gc[XX], gc[YY], gc[ZZ]);
        printf("    box vectors :%7.3f%7.3f%7.3f (nm)\n", norm(box[XX]), norm(box[YY]), norm(box[ZZ]));
        printf("    box angles  :%7.2f%7.2f%7.2f (degrees)\n",
               norm2(box[ZZ]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[YY], box[ZZ]),
               norm2(box[ZZ]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[XX], box[ZZ]),
               norm2(box[YY]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[XX], box[YY]));
        printf("    box volume  :%7.2f               (nm^3)\n", det(box));
    }

    if (bRho || bOrient || bAlign)
    {
        mass = calc_mass(&atoms, !fn2bTPX(infile), &aps);
    }

    if (bOrient)
    {
        int*  index;
        char* grpnames;

        /* Get a group for principal component analysis */
        fprintf(stderr, "\nSelect group for the determining the orientation\n");
        get_index(&atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, &grpnames);

        /* Orient the principal axes along the coordinate axes */
        orient_princ(&atoms, isize, index, natom, x, bHaveV ? v : nullptr, nullptr);
        sfree(index);
        sfree(grpnames);
    }

    if (bScale)
    {
        /* scale coordinates and box */
        if (bRho)
        {
            /* Compute scaling constant */
            real vol, dens;

            vol  = det(box);
            dens = (mass * gmx::c_amu) / (vol * gmx::c_nano * gmx::c_nano * gmx::c_nano);
            fprintf(stderr, "Volume  of input %g (nm^3)\n", vol);
            fprintf(stderr, "Mass    of input %g (a.m.u.)\n", mass);
            fprintf(stderr, "Density of input %g (g/l)\n", dens);
            if (vol == 0 || mass == 0)
            {
                gmx_fatal(FARGS,
                          "Cannot scale density with "
                          "zero mass (%g) or volume (%g)\n",
                          mass,
                          vol);
            }

            scale[XX] = scale[YY] = scale[ZZ] = std::cbrt(dens / rho);
            fprintf(stderr, "Scaling all box vectors by %g\n", scale[XX]);
        }
        scale_conf(atoms.nr, x, box, scale);
    }

    if (bAlign)
    {
        if (bIndex)
        {
            fprintf(stderr, "\nSelect a group that you want to align:\n");
            get_index(&atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &numAlignmentAtoms, &aindex, &agrpname);
        }
        else
        {
            numAlignmentAtoms = atoms.nr;
            snew(aindex, numAlignmentAtoms);
            for (i = 0; i < numAlignmentAtoms; i++)
            {
                aindex[i] = i;
            }
        }
        printf("Aligning %d atoms (out of %d) to %g %g %g, center of rotation %g %g %g\n",
               numAlignmentAtoms,
               natom,
               targetvec[XX],
               targetvec[YY],
               targetvec[ZZ],
               aligncenter[XX],
               aligncenter[YY],
               aligncenter[ZZ]);
        /*subtract out pivot point*/
        for (i = 0; i < numAlignmentAtoms; i++)
        {
            rvec_dec(x[aindex[i]], aligncenter);
        }
        /*now determine transform and rotate*/
        /*will this work?*/
        principal_comp(numAlignmentAtoms, aindex, atoms.atom, x, trans, princd);

        unitv(targetvec, targetvec);
        printf("Using %g %g %g as principal axis\n", trans[0][2], trans[1][2], trans[2][2]);
        tmpvec[XX] = trans[0][2];
        tmpvec[YY] = trans[1][2];
        tmpvec[ZZ] = trans[2][2];
        calc_rotmatrix(tmpvec, targetvec, rotmatrix);
        /* rotmatrix finished */

        for (i = 0; i < numAlignmentAtoms; ++i)
        {
            mvmul(rotmatrix, x[aindex[i]], tmpvec);
            copy_rvec(tmpvec, x[aindex[i]]);
        }

        /*add pivot point back*/
        for (i = 0; i < numAlignmentAtoms; i++)
        {
            rvec_inc(x[aindex[i]], aligncenter);
        }
        if (!bIndex)
        {
            sfree(aindex);
        }
    }

    if (bTranslate)
    {
        if (bIndex)
        {
            fprintf(stderr, "\nSelect a group that you want to translate:\n");
            get_index(&atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &ssize, &sindex, &sgrpname);
        }
        else
        {
            ssize  = atoms.nr;
            sindex = nullptr;
        }
        printf("Translating %d atoms (out of %d) by %g %g %g nm\n",
               ssize,
               natom,
               translation[XX],
               translation[YY],
               translation[ZZ]);
        if (sindex)
        {
            for (i = 0; i < ssize; i++)
            {
                rvec_inc(x[sindex[i]], translation);
            }
        }
        else
        {
            for (i = 0; i < natom; i++)
            {
                rvec_inc(x[i], translation);
            }
        }
    }
    if (bRotate)
    {
        /* Rotate */
        printf("Rotating %g, %g, %g degrees around the X, Y and Z axis respectively\n",
               rotangles[XX],
               rotangles[YY],
               rotangles[ZZ]);
        for (i = 0; i < DIM; i++)
        {
            rotangles[i] *= gmx::c_deg2Rad;
        }
        rotate_conf(natom, x, v, rotangles[XX], rotangles[YY], rotangles[ZZ]);
    }

    if (bCalcGeom)
    {
        /* recalc geometrical center and max and min coordinates and size */
        calc_geom(ssize, sindex, x, gc, rmin, rmax, FALSE);
        rvec_sub(rmax, rmin, size);
        if (bScale || bOrient || bRotate)
        {
            printf("new system size : %6.3f %6.3f %6.3f\n", size[XX], size[YY], size[ZZ]);
        }
    }

    if ((btype[0] != nullptr) && (bSetSize || bDist || (btype[0][0] == 't' && bSetAng)))
    {
        pbcType = PbcType::Xyz;
        if (!(bSetSize || bDist))
        {
            for (i = 0; i < DIM; i++)
            {
                newbox[i] = norm(box[i]);
            }
        }
        clear_mat(box);
        /* calculate new boxsize */
        switch (btype[0][0])
        {
            case 't':
                if (bDist)
                {
                    for (i = 0; i < DIM; i++)
                    {
                        newbox[i] = size[i] + 2 * dist;
                    }
                }
                if (!bSetAng)
                {
                    box[XX][XX] = newbox[XX];
                    box[YY][YY] = newbox[YY];
                    box[ZZ][ZZ] = newbox[ZZ];
                }
                else
                {
                    matrix_convert(box, newbox, newang);
                }
                break;
            case 'c':
            case 'd':
            case 'o':
                if (bSetSize)
                {
                    d = newbox[0];
                }
                else
                {
                    d = diam + 2 * dist;
                }
                if (btype[0][0] == 'c')
                {
                    for (i = 0; i < DIM; i++)
                    {
                        box[i][i] = d;
                    }
                }
                else if (btype[0][0] == 'd')
                {
                    box[XX][XX] = d;
                    box[YY][YY] = d;
                    box[ZZ][XX] = d / 2;
                    box[ZZ][YY] = d / 2;
                    box[ZZ][ZZ] = d * std::sqrt(2.0) / 2.0;
                }
                else
                {
                    box[XX][XX] = d;
                    box[YY][XX] = d / 3;
                    box[YY][YY] = d * std::sqrt(2.0) * 2.0 / 3.0;
                    box[ZZ][XX] = -d / 3;
                    box[ZZ][YY] = d * std::sqrt(2.0) / 3.0;
                    box[ZZ][ZZ] = d * std::sqrt(6.0) / 3.0;
                }
                break;
        }
    }

    /* calculate new coords for geometrical center */
    if (!bSetCenter)
    {
        calc_box_center(ecenterDEF, box, center);
    }

    /* center molecule on 'center' */
    if (bCenter)
    {
        center_conf(natom, x, center, gc);
    }

    /* print some */
    if (bCalcGeom)
    {
        calc_geom(ssize, sindex, x, gc, rmin, rmax, FALSE);
        printf("new center      :%7.3f%7.3f%7.3f (nm)\n", gc[XX], gc[YY], gc[ZZ]);
    }
    if (bOrient || bScale || bDist || bSetSize)
    {
        printf("new box vectors :%7.3f%7.3f%7.3f (nm)\n", norm(box[XX]), norm(box[YY]), norm(box[ZZ]));
        printf("new box angles  :%7.2f%7.2f%7.2f (degrees)\n",
               norm2(box[ZZ]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[YY], box[ZZ]),
               norm2(box[ZZ]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[XX], box[ZZ]),
               norm2(box[YY]) == 0 ? 0 : gmx::c_rad2Deg * gmx_angle(box[XX], box[YY]));
        printf("new box volume  :%7.2f               (nm^3)\n", det(box));
    }

    if (check_box(PbcType::Xyz, box))
    {
        printf("\nWARNING: %s\n"
               "See the GROMACS manual for a description of the requirements that\n"
               "must be satisfied by descriptions of simulation cells.\n",
               check_box(PbcType::Xyz, box));
    }

    if (bDist && btype[0][0] == 't')
    {
        if (TRICLINIC(box))
        {
            printf("\nWARNING: Your box is triclinic with non-orthogonal axes. In this case, the\n"
                   "distance from the solute to a box surface along the corresponding normal\n"
                   "vector might be somewhat smaller than your specified value %f.\n"
                   "You can check the actual value with g_mindist -pi\n",
                   dist);
        }
        else if (!opt2parg_bSet("-bt", NPA, pa))
        {
            printf("\nWARNING: No boxtype specified - distance condition applied in each "
                   "dimension.\n"
                   "If the molecule rotates the actual distance will be smaller. You might want\n"
                   "to use a cubic box instead, or why not try a dodecahedron today?\n");
        }
    }
    /* We have already checked that the output is a pdb file and the input a tpr file */
    if (bCONECT)
    {
        conect = gmx_conect_generate(top);
    }
    else
    {
        conect = nullptr;
    }

    if (bIndex)
    {
        fprintf(stderr, "\nSelect a group for output:\n");
        get_index(&atoms, opt2fn_null("-n", NFILE, fnm), 1, &isize, &index, &grpname);

        if (resnr_start >= 0)
        {
            renum_resnr(&atoms, isize, index, resnr_start);
        }

        if (opt2parg_bSet("-label", NPA, pa))
        {
            for (i = 0; (i < atoms.nr); i++)
            {
                atoms.resinfo[atoms.atom[i].resind].chainid = label[0];
            }
        }

        if (opt2bSet("-bf", NFILE, fnm) || bLegend)
        {
            gmx_fatal(FARGS, "Sorry, cannot do bfactors with an index group.");
        }

        if (outftp == efPDB)
        {
            out = gmx_ffopen(outfile, "w");
            write_pdbfile_indexed(out, name, &atoms, x, pbcType, box, ' ', 1, isize, index, conect, FALSE);
            gmx_ffclose(out);
        }
        else
        {
            write_sto_conf_indexed(
                    outfile, name, &atoms, x, bHaveV ? v : nullptr, pbcType, box, isize, index);
        }
        sfree(grpname);
        sfree(index);
    }
    else
    {
        if (resnr_start >= 0)
        {
            renum_resnr(&atoms, atoms.nr, nullptr, resnr_start);
        }

        if ((outftp == efPDB) || (outftp == efPQR))
        {
            out = gmx_ffopen(outfile, "w");
            if (bMead)
            {
                fprintf(out,
                        "REMARK    "
                        "The B-factors in this file hold atomic radii\n");
                fprintf(out,
                        "REMARK    "
                        "The occupancy in this file hold atomic charges\n");
            }
            else if (bGrasp)
            {
                fprintf(out, "GRASP PDB FILE\nFORMAT NUMBER=1\n");
                fprintf(out,
                        "REMARK    "
                        "The B-factors in this file hold atomic charges\n");
                fprintf(out,
                        "REMARK    "
                        "The occupancy in this file hold atomic radii\n");
            }
            else if (opt2bSet("-bf", NFILE, fnm))
            {
                read_bfac(opt2fn("-bf", NFILE, fnm), &n_bfac, &bfac, &bfac_nr);
                set_pdb_conf_bfac(atoms.nr, atoms.nres, &atoms, n_bfac, bfac, bfac_nr, peratom);
            }
            if (opt2parg_bSet("-label", NPA, pa))
            {
                for (i = 0; (i < atoms.nr); i++)
                {
                    atoms.resinfo[atoms.atom[i].resind].chainid = label[0];
                }
            }
            /* Need to bypass the regular write_pdbfile because I don't want to change
             * all instances to include the boolean flag for writing out PQR files.
             */
            int* index;
            snew(index, atoms.nr);
            for (int i = 0; i < atoms.nr; i++)
            {
                index[i] = i;
            }
            write_pdbfile_indexed(
                    out, name, &atoms, x, pbcType, box, ' ', -1, atoms.nr, index, conect, outftp == efPQR);
            sfree(index);
            if (bLegend)
            {
                pdb_legend(out, atoms.nr, atoms.nres, &atoms, x);
            }
            if (visbox[0] > 0)
            {
                visualize_box(out,
                              bLegend ? atoms.nr + 12 : atoms.nr,
                              bLegend ? atoms.nres = 12 : atoms.nres,
                              box,
                              visbox);
            }
            gmx_ffclose(out);
        }
        else
        {
            write_sto_conf(outfile, name, &atoms, x, bHaveV ? v : nullptr, pbcType, box);
        }
    }
    done_atom(&atoms);
    done_symtab(&symtab);
    sfree(name);
    if (x)
    {
        sfree(x);
    }
    if (v)
    {
        sfree(v);
    }
    do_view(oenv, outfile, nullptr);
    output_env_done(oenv);

    return 0;
}
