/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

enum {
    euSel, euRect, euTric, euCompact, euNR
};


static void calc_pbc_cluster(int ecenter, int nrefat, t_topology *top, int ePBC,
                             rvec x[], int index[], matrix box)
{
    int       m, i, j, j0, j1, jj, ai, aj;
    int       imin, jmin;
    real      fac, min_dist2;
    rvec      dx, xtest, box_center;
    int       nmol, imol_center;
    int      *molind;
    gmx_bool *bMol, *bTmp;
    rvec     *m_com, *m_shift;
    t_pbc     pbc;
    int      *cluster;
    int      *added;
    int       ncluster, nadded;
    real      tmp_r2;

    calc_box_center(ecenter, box, box_center);

    /* Initiate the pbc structure */
    std::memset(&pbc, 0, sizeof(pbc));
    set_pbc(&pbc, ePBC, box);

    /* Convert atom index to molecular */
    nmol   = top->mols.nr;
    molind = top->mols.index;
    snew(bMol, nmol);
    snew(m_com, nmol);
    snew(m_shift, nmol);
    snew(cluster, nmol);
    snew(added, nmol);
    snew(bTmp, top->atoms.nr);

    for (i = 0; (i < nrefat); i++)
    {
        /* Mark all molecules in the index */
        ai       = index[i];
        bTmp[ai] = TRUE;
        /* Binary search assuming the molecules are sorted */
        j0 = 0;
        j1 = nmol-1;
        while (j0 < j1)
        {
            if (ai < molind[j0+1])
            {
                j1 = j0;
            }
            else if (ai >= molind[j1])
            {
                j0 = j1;
            }
            else
            {
                jj = (j0+j1)/2;
                if (ai < molind[jj+1])
                {
                    j1 = jj;
                }
                else
                {
                    j0 = jj;
                }
            }
        }
        bMol[j0] = TRUE;
    }
    /* Double check whether all atoms in all molecules that are marked are part
     * of the cluster. Simultaneously compute the center of geometry.
     */
    min_dist2   = 10*gmx::square(trace(box));
    imol_center = -1;
    ncluster    = 0;
    for (i = 0; i < nmol; i++)
    {
        for (j = molind[i]; j < molind[i+1]; j++)
        {
            if (bMol[i] && !bTmp[j])
            {
                gmx_fatal(FARGS, "Molecule %d marked for clustering but not atom %d in it - check your index!", i+1, j+1);
            }
            else if (!bMol[i] && bTmp[j])
            {
                gmx_fatal(FARGS, "Atom %d marked for clustering but not molecule %d - this is an internal error...", j+1, i+1);
            }
            else if (bMol[i])
            {
                /* Make molecule whole, move 2nd and higher atom to same periodicity as 1st atom in molecule */
                if (j > molind[i])
                {
                    pbc_dx(&pbc, x[j], x[j-1], dx);
                    rvec_add(x[j-1], dx, x[j]);
                }
                /* Compute center of geometry of molecule - m_com[i] was zeroed when we did snew() on it! */
                rvec_inc(m_com[i], x[j]);
            }
        }
        if (bMol[i])
        {
            /* Normalize center of geometry */
            fac = 1.0/(molind[i+1]-molind[i]);
            for (m = 0; (m < DIM); m++)
            {
                m_com[i][m] *= fac;
            }
            /* Determine which molecule is closest to the center of the box */
            pbc_dx(&pbc, box_center, m_com[i], dx);
            tmp_r2 = iprod(dx, dx);

            if (tmp_r2 < min_dist2)
            {
                min_dist2   = tmp_r2;
                imol_center = i;
            }
            cluster[ncluster++] = i;
        }
    }
    sfree(bTmp);

    if (ncluster <= 0)
    {
        fprintf(stderr, "No molecules selected in the cluster\n");
        return;
    }
    else if (imol_center == -1)
    {
        fprintf(stderr, "No central molecules could be found\n");
        return;
    }

    nadded            = 0;
    added[nadded++]   = imol_center;
    bMol[imol_center] = FALSE;

    while (nadded < ncluster)
    {
        /* Find min distance between cluster molecules and those remaining to be added */
        min_dist2   = 10*gmx::square(trace(box));
        imin        = -1;
        jmin        = -1;
        /* Loop over added mols */
        for (i = 0; i < nadded; i++)
        {
            ai = added[i];
            /* Loop over all mols */
            for (j = 0; j < ncluster; j++)
            {
                aj = cluster[j];
                /* check those remaining to be added */
                if (bMol[aj])
                {
                    pbc_dx(&pbc, m_com[aj], m_com[ai], dx);
                    tmp_r2 = iprod(dx, dx);
                    if (tmp_r2 < min_dist2)
                    {
                        min_dist2   = tmp_r2;
                        imin        = ai;
                        jmin        = aj;
                    }
                }
            }
        }

        /* Add the best molecule */
        added[nadded++]   = jmin;
        bMol[jmin]        = FALSE;
        /* Calculate the shift from the ai molecule */
        pbc_dx(&pbc, m_com[jmin], m_com[imin], dx);
        rvec_add(m_com[imin], dx, xtest);
        rvec_sub(xtest, m_com[jmin], m_shift[jmin]);
        rvec_inc(m_com[jmin], m_shift[jmin]);

        for (j = molind[jmin]; j < molind[jmin+1]; j++)
        {
            rvec_inc(x[j], m_shift[jmin]);
        }
        fprintf(stdout, "\rClustering iteration %d of %d...", nadded, ncluster);
        fflush(stdout);
    }

    sfree(added);
    sfree(cluster);
    sfree(bMol);
    sfree(m_com);
    sfree(m_shift);

    fprintf(stdout, "\n");
}

static void put_molecule_com_in_box(int unitcell_enum, int ecenter,
                                    t_block *mols,
                                    int natoms, t_atom atom[],
                                    int ePBC, matrix box, rvec x[])
{
    int     i, j;
    int     d;
    rvec    com, new_com, shift, box_center;
    real    m;
    double  mtot;
    t_pbc   pbc;

    calc_box_center(ecenter, box, box_center);
    set_pbc(&pbc, ePBC, box);
    if (mols->nr <= 0)
    {
        gmx_fatal(FARGS, "There are no molecule descriptions. I need a .tpr file for this pbc option.");
    }
    for (i = 0; (i < mols->nr); i++)
    {
        /* calc COM */
        clear_rvec(com);
        mtot = 0;
        for (j = mols->index[i]; (j < mols->index[i+1] && j < natoms); j++)
        {
            m = atom[j].m;
            for (d = 0; d < DIM; d++)
            {
                com[d] += m*x[j][d];
            }
            mtot += m;
        }
        /* calculate final COM */
        svmul(1.0/mtot, com, com);

        /* check if COM is outside box */
        copy_rvec(com, new_com);
        switch (unitcell_enum)
        {
            case euRect:
                put_atoms_in_box(ePBC, box, 1, &new_com);
                break;
            case euTric:
                put_atoms_in_triclinic_unitcell(ecenter, box, 1, &new_com);
                break;
            case euCompact:
                put_atoms_in_compact_unitcell(ePBC, ecenter, box, 1, &new_com);
                break;
        }
        rvec_sub(new_com, com, shift);
        if (norm2(shift) > 0)
        {
            if (debug)
            {
                fprintf(debug, "\nShifting position of molecule %d "
                        "by %8.3f  %8.3f  %8.3f\n", i+1,
                        shift[XX], shift[YY], shift[ZZ]);
            }
            for (j = mols->index[i]; (j < mols->index[i+1] && j < natoms); j++)
            {
                rvec_inc(x[j], shift);
            }
        }
    }
}

static void put_residue_com_in_box(int unitcell_enum, int ecenter,
                                   int natoms, t_atom atom[],
                                   int ePBC, matrix box, rvec x[])
{
    int              i, j, res_start, res_end;
    int              d, presnr;
    real             m;
    double           mtot;
    rvec             box_center, com, new_com, shift;
    static const int NOTSET = -12347;
    calc_box_center(ecenter, box, box_center);

    presnr    = NOTSET;
    res_start = 0;
    clear_rvec(com);
    mtot = 0;
    for (i = 0; i < natoms+1; i++)
    {
        if (i == natoms || (presnr != atom[i].resind && presnr != NOTSET))
        {
            /* calculate final COM */
            res_end = i;
            svmul(1.0/mtot, com, com);

            /* check if COM is outside box */
            copy_rvec(com, new_com);
            switch (unitcell_enum)
            {
                case euRect:
                    put_atoms_in_box(ePBC, box, 1, &new_com);
                    break;
                case euTric:
                    put_atoms_in_triclinic_unitcell(ecenter, box, 1, &new_com);
                    break;
                case euCompact:
                    put_atoms_in_compact_unitcell(ePBC, ecenter, box, 1, &new_com);
                    break;
            }
            rvec_sub(new_com, com, shift);
            if (norm2(shift))
            {
                if (debug)
                {
                    fprintf(debug, "\nShifting position of residue %d (atoms %d-%d) "
                            "by %g,%g,%g\n", atom[res_start].resind+1,
                            res_start+1, res_end+1, shift[XX], shift[YY], shift[ZZ]);
                }
                for (j = res_start; j < res_end; j++)
                {
                    rvec_inc(x[j], shift);
                }
            }
            clear_rvec(com);
            mtot = 0;

            /* remember start of new residue */
            res_start = i;
        }
        if (i < natoms)
        {
            /* calc COM */
            m = atom[i].m;
            for (d = 0; d < DIM; d++)
            {
                com[d] += m*x[i][d];
            }
            mtot += m;

            presnr = atom[i].resind;
        }
    }
}

static void center_x(int ecenter, rvec x[], matrix box, int n, int nc, int ci[])
{
    int  i, m, ai;
    rvec cmin, cmax, box_center, dx;

    if (nc > 0)
    {
        copy_rvec(x[ci[0]], cmin);
        copy_rvec(x[ci[0]], cmax);
        for (i = 0; i < nc; i++)
        {
            ai = ci[i];
            for (m = 0; m < DIM; m++)
            {
                if (x[ai][m] < cmin[m])
                {
                    cmin[m] = x[ai][m];
                }
                else if (x[ai][m] > cmax[m])
                {
                    cmax[m] = x[ai][m];
                }
            }
        }
        calc_box_center(ecenter, box, box_center);
        for (m = 0; m < DIM; m++)
        {
            dx[m] = box_center[m]-(cmin[m]+cmax[m])*0.5;
        }

        for (i = 0; i < n; i++)
        {
            rvec_inc(x[i], dx);
        }
    }
}

static void mk_filenm(char *base, const char *ext, int ndigit, int file_nr,
                      char out_file[])
{
    char nbuf[128];
    int  nd = 0, fnr;

    std::strcpy(out_file, base);
    fnr = file_nr;
    do
    {
        fnr /= 10;
        nd++;
    }
    while (fnr > 0);

    if (nd < ndigit)
    {
        std::strncat(out_file, "00000000000", ndigit-nd);
    }
    sprintf(nbuf, "%d.", file_nr);
    std::strcat(out_file, nbuf);
    std::strcat(out_file, ext);
}

void check_trr(const char *fn)
{
    if (fn2ftp(fn) != efTRR)
    {
        gmx_fatal(FARGS, "%s is not a trajectory file, exiting\n", fn);
    }
}

void do_trunc(const char *fn, real t0)
{
    t_fileio        *in;
    FILE            *fp;
    gmx_bool         bStop, bOK;
    gmx_trr_header_t sh;
    gmx_off_t        fpos;
    char             yesno[256];
    int              j;
    real             t = 0;

    if (t0 == -1)
    {
        gmx_fatal(FARGS, "You forgot to set the truncation time");
    }

    /* Check whether this is a .trr file */
    check_trr(fn);

    in   = gmx_trr_open(fn, "r");
    fp   = gmx_fio_getfp(in);
    if (fp == nullptr)
    {
        fprintf(stderr, "Sorry, can not trunc %s, truncation of this filetype is not supported\n", fn);
        gmx_trr_close(in);
    }
    else
    {
        j     = 0;
        fpos  = gmx_fio_ftell(in);
        bStop = FALSE;
        while (!bStop && gmx_trr_read_frame_header(in, &sh, &bOK))
        {
            gmx_trr_read_frame_data(in, &sh, nullptr, nullptr, nullptr, nullptr);
            fpos = gmx_ftell(fp);
            t    = sh.t;
            if (t >= t0)
            {
                gmx_fseek(fp, fpos, SEEK_SET);
                bStop = TRUE;
            }
        }
        if (bStop)
        {
            fprintf(stderr, "Do you REALLY want to truncate this trajectory (%s) at:\n"
                    "frame %d, time %g, bytes %ld ??? (type YES if so)\n",
                    fn, j, t, (long int)fpos);
            if (1 != scanf("%s", yesno))
            {
                gmx_fatal(FARGS, "Error reading user input");
            }
            if (std::strcmp(yesno, "YES") == 0)
            {
                fprintf(stderr, "Once again, I'm gonna DO this...\n");
                gmx_trr_close(in);
                if (0 != gmx_truncate(fn, fpos))
                {
                    gmx_fatal(FARGS, "Error truncating file %s", fn);
                }
            }
            else
            {
                fprintf(stderr, "Ok, I'll forget about it\n");
            }
        }
        else
        {
            fprintf(stderr, "Already at end of file (t=%g)...\n", t);
            gmx_trr_close(in);
        }
    }
}

/*! \brief Read a full molecular topology if useful and available.
 *
 * If the input trajectory file is not in TNG format, and the output
 * file is in TNG format, then we want to try to read a full topology
 * (if available), so that we can write molecule information to the
 * output file. The full topology provides better molecule information
 * than is available from the normal t_topology data used by GROMACS
 * tools.
 *
 * Also, the t_topology is only read under (different) particular
 * conditions. If both apply, then a .tpr file might be read
 * twice. Trying to fix this redundancy while trjconv is still an
 * all-purpose tool does not seem worthwhile.
 *
 * Because of the way gmx_prepare_tng_writing is implemented, the case
 * where the input TNG file has no molecule information will never
 * lead to an output TNG file having molecule information. Since
 * molecule information will generally be present if the input TNG
 * file was written by a GROMACS tool, this seems like reasonable
 * behaviour. */
static gmx_mtop_t *read_mtop_for_tng(const char *tps_file,
                                     const char *input_file,
                                     const char *output_file)
{
    gmx_mtop_t *mtop = nullptr;

    if (fn2bTPX(tps_file) &&
        efTNG != fn2ftp(input_file) &&
        efTNG == fn2ftp(output_file))
    {
        int temp_natoms = -1;
        snew(mtop, 1);
        read_tpx(tps_file, nullptr, nullptr, &temp_natoms,
                 nullptr, nullptr, mtop);
    }

    return mtop;
}

int gmx_trjconv(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] can convert trajectory files in many ways:",
        "",
        "* from one format to another",
        "* select a subset of atoms",
        "* change the periodicity representation",
        "* keep multimeric molecules together",
        "* center atoms in the box",
        "* fit atoms to reference structure",
        "* reduce the number of frames",
        "* change the timestamps of the frames ([TT]-t0[tt] and [TT]-timestep[tt])",
        "* cut the trajectory in small subtrajectories according",
        "  to information in an index file. This allows subsequent analysis of",
        "  the subtrajectories that could, for example, be the result of a",
        "  cluster analysis. Use option [TT]-sub[tt].",
        "  This assumes that the entries in the index file are frame numbers and",
        "  dumps each group in the index file to a separate trajectory file.",
        "* select frames within a certain range of a quantity given",
        "  in an [REF].xvg[ref] file.",
        "",
        "[gmx-trjcat] is better suited for concatenating multiple trajectory files.",
        "[PAR]",

        "The following formats are supported for input and output:",
        "[REF].xtc[ref], [REF].trr[ref], [REF].gro[ref], [TT].g96[tt]",
        "and [REF].pdb[ref].",
        "The file formats are detected from the file extension.",
        "The precision of [REF].xtc[ref] and [REF].gro[ref] output is taken from the",
        "input file for [REF].xtc[ref], [REF].gro[ref] and [REF].pdb[ref],",
        "and from the [TT]-ndec[tt] option for other input formats. The precision",
        "is always taken from [TT]-ndec[tt], when this option is set.",
        "All other formats have fixed precision. [REF].trr[ref]",
        "output can be single or double precision, depending on the precision",
        "of the [THISMODULE] binary.",
        "Note that velocities are only supported in",
        "[REF].trr[ref], [REF].gro[ref] and [TT].g96[tt] files.[PAR]",

        "Option [TT]-sep[tt] can be used to write every frame to a separate",
        "[TT].gro, .g96[tt] or [REF].pdb[ref] file. By default, all frames all written to one file.",
        "[REF].pdb[ref] files with all frames concatenated can be viewed with",
        "[TT]rasmol -nmrpdb[tt].[PAR]",

        "It is possible to select part of your trajectory and write it out",
        "to a new trajectory file in order to save disk space, e.g. for leaving",
        "out the water from a trajectory of a protein in water.",
        "[BB]ALWAYS[bb] put the original trajectory on tape!",
        "We recommend to use the portable [REF].xtc[ref] format for your analysis",
        "to save disk space and to have portable files.[PAR]",

        "There are two options for fitting the trajectory to a reference",
        "either for essential dynamics analysis, etc.",
        "The first option is just plain fitting to a reference structure",
        "in the structure file. The second option is a progressive fit",
        "in which the first timeframe is fitted to the reference structure ",
        "in the structure file to obtain and each subsequent timeframe is ",
        "fitted to the previously fitted structure. This way a continuous",
        "trajectory is generated, which might not be the case when using the",
        "regular fit method, e.g. when your protein undergoes large",
        "conformational transitions.[PAR]",

        "Option [TT]-pbc[tt] sets the type of periodic boundary condition",
        "treatment:",
        "",
        " * [TT]mol[tt] puts the center of mass of molecules in the box,",
        "   and requires a run input file to be supplied with [TT]-s[tt].",
        " * [TT]res[tt] puts the center of mass of residues in the box.",
        " * [TT]atom[tt] puts all the atoms in the box.",
        " * [TT]nojump[tt] checks if atoms jump across the box and then puts",
        "   them back. This has the effect that all molecules",
        "   will remain whole (provided they were whole in the initial",
        "   conformation). [BB]Note[bb] that this ensures a continuous trajectory but",
        "   molecules may diffuse out of the box. The starting configuration",
        "   for this procedure is taken from the structure file, if one is",
        "   supplied, otherwise it is the first frame.",
        " * [TT]cluster[tt] clusters all the atoms in the selected index",
        "   such that they are all closest to the center of mass of the cluster,",
        "   which is iteratively updated. [BB]Note[bb] that this will only give meaningful",
        "   results if you in fact have a cluster. Luckily that can be checked",
        "   afterwards using a trajectory viewer. Note also that if your molecules",
        "   are broken this will not work either.",
        " * [TT]whole[tt] only makes broken molecules whole.",
        "",

        "Option [TT]-ur[tt] sets the unit cell representation for options",
        "[TT]mol[tt], [TT]res[tt] and [TT]atom[tt] of [TT]-pbc[tt].",
        "All three options give different results for triclinic boxes and",
        "identical results for rectangular boxes.",
        "[TT]rect[tt] is the ordinary brick shape.",
        "[TT]tric[tt] is the triclinic unit cell.",
        "[TT]compact[tt] puts all atoms at the closest distance from the center",
        "of the box. This can be useful for visualizing e.g. truncated octahedra",
        "or rhombic dodecahedra. The center for options [TT]tric[tt] and [TT]compact[tt]",
        "is [TT]tric[tt] (see below), unless the option [TT]-boxcenter[tt]",
        "is set differently.[PAR]",

        "Option [TT]-center[tt] centers the system in the box. The user can",
        "select the group which is used to determine the geometrical center.",
        "Option [TT]-boxcenter[tt] sets the location of the center of the box",
        "for options [TT]-pbc[tt] and [TT]-center[tt]. The center options are:",
        "[TT]tric[tt]: half of the sum of the box vectors,",
        "[TT]rect[tt]: half of the box diagonal,",
        "[TT]zero[tt]: zero.",
        "Use option [TT]-pbc mol[tt] in addition to [TT]-center[tt] when you",
        "want all molecules in the box after the centering.[PAR]",

        "Option [TT]-box[tt] sets the size of the new box. This option only works",
        "for leading dimensions and is thus generally only useful for rectangular boxes.",
        "If you want to modify only some of the dimensions, e.g. when reading from",
        "a trajectory, you can use -1 for those dimensions that should stay the same",

        "It is not always possible to use combinations of [TT]-pbc[tt],",
        "[TT]-fit[tt], [TT]-ur[tt] and [TT]-center[tt] to do exactly what",
        "you want in one call to [THISMODULE]. Consider using multiple",
        "calls, and check out the GROMACS website for suggestions.[PAR]",

        "With [TT]-dt[tt], it is possible to reduce the number of ",
        "frames in the output. This option relies on the accuracy of the times",
        "in your input trajectory, so if these are inaccurate use the",
        "[TT]-timestep[tt] option to modify the time (this can be done",
        "simultaneously). For making smooth movies, the program [gmx-filter]",
        "can reduce the number of frames while using low-pass frequency",
        "filtering, this reduces aliasing of high frequency motions.[PAR]",

        "Using [TT]-trunc[tt] [THISMODULE] can truncate [REF].trr[ref] in place, i.e.",
        "without copying the file. This is useful when a run has crashed",
        "during disk I/O (i.e. full disk), or when two contiguous",
        "trajectories must be concatenated without having double frames.[PAR]",

        "Option [TT]-dump[tt] can be used to extract a frame at or near",
        "one specific time from your trajectory, but only works reliably",
        "if the time interval between frames is uniform.[PAR]",

        "Option [TT]-drop[tt] reads an [REF].xvg[ref] file with times and values.",
        "When options [TT]-dropunder[tt] and/or [TT]-dropover[tt] are set,",
        "frames with a value below and above the value of the respective options",
        "will not be written."
    };

    int         pbc_enum;
    enum
    {
        epSel,
        epNone,
        epComMol,
        epComRes,
        epComAtom,
        epNojump,
        epCluster,
        epWhole,
        epNR
    };
    const char *pbc_opt[epNR + 1] =
    {
        nullptr, "none", "mol", "res", "atom", "nojump", "cluster", "whole",
        nullptr
    };

    int         unitcell_enum;
    const char *unitcell_opt[euNR+1] =
    { nullptr, "rect", "tric", "compact", nullptr };

    enum
    {
        ecSel, ecTric, ecRect, ecZero, ecNR
    };
    const char *center_opt[ecNR+1] =
    { nullptr, "tric", "rect", "zero", nullptr };
    int         ecenter;

    int         fit_enum;
    enum
    {
        efSel, efNone, efFit, efFitXY, efReset, efResetXY, efPFit, efNR
    };
    const char *fit[efNR + 1] =
    {
        nullptr, "none", "rot+trans", "rotxy+transxy", "translation", "transxy",
        "progressive", nullptr
    };

    static gmx_bool  bSeparate     = FALSE, bVels = TRUE, bForce = FALSE, bCONECT = FALSE;
    static gmx_bool  bCenter       = FALSE;
    static int       skip_nr       = 1, ndec = 3, nzero = 0;
    static real      tzero         = 0, delta_t = 0, timestep = 0, ttrunc = -1, tdump = -1, split_t = 0;
    static rvec      newbox        = {0, 0, 0}, shift = {0, 0, 0}, trans = {0, 0, 0};
    static char     *exec_command  = nullptr;
    static real      dropunder     = 0, dropover = 0;
    static gmx_bool  bRound        = FALSE;

    t_pargs
        pa[] =
    {
        { "-skip", FALSE, etINT,
          { &skip_nr }, "Only write every nr-th frame" },
        { "-dt", FALSE, etTIME,
          { &delta_t },
          "Only write frame when t MOD dt = first time (%t)" },
        { "-round", FALSE, etBOOL,
          { &bRound }, "Round measurements to nearest picosecond"},
        { "-dump", FALSE, etTIME,
          { &tdump }, "Dump frame nearest specified time (%t)" },
        { "-t0", FALSE, etTIME,
          { &tzero },
          "Starting time (%t) (default: don't change)" },
        { "-timestep", FALSE, etTIME,
          { &timestep },
          "Change time step between input frames (%t)" },
        { "-pbc", FALSE, etENUM,
          { pbc_opt },
          "PBC treatment (see help text for full description)" },
        { "-ur", FALSE, etENUM,
          { unitcell_opt }, "Unit-cell representation" },
        { "-center", FALSE, etBOOL,
          { &bCenter }, "Center atoms in box" },
        { "-boxcenter", FALSE, etENUM,
          { center_opt }, "Center for -pbc and -center" },
        { "-box", FALSE, etRVEC,
          { newbox },
          "Size for new cubic box (default: read from input)" },
        { "-trans", FALSE, etRVEC,
          { trans },
          "All coordinates will be translated by trans. This "
          "can advantageously be combined with -pbc mol -ur "
          "compact." },
        { "-shift", FALSE, etRVEC,
          { shift },
          "All coordinates will be shifted by framenr*shift" },
        { "-fit", FALSE, etENUM,
          { fit },
          "Fit molecule to ref structure in the structure file" },
        { "-ndec", FALSE, etINT,
          { &ndec },
          "Precision for .xtc and .gro writing in number of "
          "decimal places" },
        { "-vel", FALSE, etBOOL,
          { &bVels }, "Read and write velocities if possible" },
        { "-force", FALSE, etBOOL,
          { &bForce }, "Read and write forces if possible" },
        { "-trunc", FALSE, etTIME,
          { &ttrunc },
          "Truncate input trajectory file after this time (%t)" },
        { "-exec", FALSE, etSTR,
          { &exec_command },
          "Execute command for every output frame with the "
          "frame number as argument" },
        { "-split", FALSE, etTIME,
          { &split_t },
          "Start writing new file when t MOD split = first "
          "time (%t)" },
        { "-sep", FALSE, etBOOL,
          { &bSeparate },
          "Write each frame to a separate .gro, .g96 or .pdb "
          "file" },
        { "-nzero", FALSE, etINT,
          { &nzero },
          "If the -sep flag is set, use these many digits "
          "for the file numbers and prepend zeros as needed" },
        { "-dropunder", FALSE, etREAL,
          { &dropunder }, "Drop all frames below this value" },
        { "-dropover", FALSE, etREAL,
          { &dropover }, "Drop all frames above this value" },
        { "-conect", FALSE, etBOOL,
          { &bCONECT },
          "Add conect records when writing [REF].pdb[ref] files. Useful "
          "for visualization of non-standard molecules, e.g. "
          "coarse grained ones" }
    };
#define NPA asize(pa)

    FILE             *out    = nullptr;
    t_trxstatus      *trxout = nullptr;
    t_trxstatus      *trxin;
    int               file_nr;
    t_trxframe        fr, frout;
    int               flags;
    rvec             *xmem  = nullptr, *vmem = nullptr, *fmem = nullptr;
    rvec             *xp    = nullptr, x_shift, hbox;
    real             *w_rls = nullptr;
    int               m, i, d, frame, outframe, natoms, nout, ncent, newstep = 0, model_nr;
#define SKIP 10
    t_topology        top;
    gmx_mtop_t       *mtop  = nullptr;
    gmx_conect        gc    = nullptr;
    int               ePBC  = -1;
    t_atoms          *atoms = nullptr, useatoms;
    matrix            top_box;
    int              *index, *cindex;
    char             *grpnm;
    int              *frindex, nrfri;
    char             *frname;
    int               ifit, my_clust = -1;
    int              *ind_fit;
    char             *gn_fit;
    t_cluster_ndx    *clust           = nullptr;
    t_trxstatus     **clust_status    = nullptr;
    int              *clust_status_id = nullptr;
    int               ntrxopen        = 0;
    int              *nfwritten       = nullptr;
    int               ndrop           = 0, ncol, drop0 = 0, drop1 = 0, dropuse = 0;
    double          **dropval;
    real              tshift = 0, t0 = -1, dt = 0.001, prec;
    gmx_bool          bFit, bPFit, bReset;
    int               nfitdim;
    gmx_rmpbc_t       gpbc = nullptr;
    gmx_bool          bRmPBC, bPBCWhole, bPBCcomRes, bPBCcomMol, bPBCcomAtom, bPBC, bNoJump, bCluster;
    gmx_bool          bCopy, bDoIt, bIndex, bTDump, bSetTime, bTPS = FALSE, bDTset = FALSE;
    gmx_bool          bExec, bTimeStep = FALSE, bDumpFrame = FALSE, bSetPrec, bNeedPrec;
    gmx_bool          bHaveFirstFrame, bHaveNextFrame, bSetBox, bSetUR, bSplit = FALSE;
    gmx_bool          bSubTraj = FALSE, bDropUnder = FALSE, bDropOver = FALSE, bTrans = FALSE;
    gmx_bool          bWriteFrame, bSplitHere;
    const char       *top_file, *in_file, *out_file = nullptr;
    char              out_file2[256], *charpt;
    char             *outf_base = nullptr;
    const char       *outf_ext  = nullptr;
    char              top_title[256], title[256], filemode[5];
    gmx_output_env_t *oenv;

    t_filenm          fnm[] = {
        { efTRX, "-f",   nullptr,      ffREAD  },
        { efTRO, "-o",   nullptr,      ffWRITE },
        { efTPS, nullptr,   nullptr,      ffOPTRD },
        { efNDX, nullptr,   nullptr,      ffOPTRD },
        { efNDX, "-fr",  "frames",  ffOPTRD },
        { efNDX, "-sub", "cluster", ffOPTRD },
        { efXVG, "-drop", "drop",    ffOPTRD }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW |
                           PCA_TIME_UNIT,
                           NFILE, fnm, NPA, pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    top_file = ftp2fn(efTPS, NFILE, fnm);
    init_top(&top);

    /* Check command line */
    in_file = opt2fn("-f", NFILE, fnm);

    if (ttrunc != -1)
    {
        do_trunc(in_file, ttrunc);
    }
    else
    {
        /* mark active cmdline options */
        bSetBox    = opt2parg_bSet("-box", NPA, pa);
        bSetTime   = opt2parg_bSet("-t0", NPA, pa);
        bSetPrec   = opt2parg_bSet("-ndec", NPA, pa);
        bSetUR     = opt2parg_bSet("-ur", NPA, pa);
        bExec      = opt2parg_bSet("-exec", NPA, pa);
        bTimeStep  = opt2parg_bSet("-timestep", NPA, pa);
        bTDump     = opt2parg_bSet("-dump", NPA, pa);
        bDropUnder = opt2parg_bSet("-dropunder", NPA, pa);
        bDropOver  = opt2parg_bSet("-dropover", NPA, pa);
        bTrans     = opt2parg_bSet("-trans", NPA, pa);
        bSplit     = (split_t != 0);

        /* parse enum options */
        fit_enum      = nenum(fit);
        bFit          = (fit_enum == efFit || fit_enum == efFitXY);
        bReset        = (fit_enum == efReset || fit_enum == efResetXY);
        bPFit         = fit_enum == efPFit;
        pbc_enum      = nenum(pbc_opt);
        bPBCWhole     = pbc_enum == epWhole;
        bPBCcomRes    = pbc_enum == epComRes;
        bPBCcomMol    = pbc_enum == epComMol;
        bPBCcomAtom   = pbc_enum == epComAtom;
        bNoJump       = pbc_enum == epNojump;
        bCluster      = pbc_enum == epCluster;
        bPBC          = pbc_enum != epNone;
        unitcell_enum = nenum(unitcell_opt);
        ecenter       = nenum(center_opt) - ecTric;

        /* set and check option dependencies */
        if (bPFit)
        {
            bFit = TRUE;        /* for pfit, fit *must* be set */
        }
        if (bFit)
        {
            bReset = TRUE;       /* for fit, reset *must* be set */
        }
        nfitdim = 0;
        if (bFit || bReset)
        {
            nfitdim = (fit_enum == efFitXY || fit_enum == efResetXY) ? 2 : 3;
        }
        bRmPBC = bFit || bPBCWhole || bPBCcomRes || bPBCcomMol;

        if (bSetUR)
        {
            if (!(bPBCcomRes || bPBCcomMol ||  bPBCcomAtom))
            {
                fprintf(stderr,
                        "WARNING: Option for unitcell representation (-ur %s)\n"
                        "         only has effect in combination with -pbc %s, %s or %s.\n"
                        "         Ingoring unitcell representation.\n\n",
                        unitcell_opt[0], pbc_opt[2], pbc_opt[3], pbc_opt[4]);
            }
        }
        if (bFit && bPBC)
        {
            gmx_fatal(FARGS, "PBC condition treatment does not work together with rotational fit.\n"
                      "Please do the PBC condition treatment first and then run trjconv in a second step\n"
                      "for the rotational fit.\n"
                      "First doing the rotational fit and then doing the PBC treatment gives incorrect\n"
                      "results!");
        }

        /* ndec is in nr of decimal places, prec is a multiplication factor: */
        prec = 1;
        for (i = 0; i < ndec; i++)
        {
            prec *= 10;
        }

        bIndex = ftp2bSet(efNDX, NFILE, fnm);


        /* Determine output type */
        out_file = opt2fn("-o", NFILE, fnm);
        int ftp  = fn2ftp(out_file);
        fprintf(stderr, "Will write %s: %s\n", ftp2ext(ftp), ftp2desc(ftp));
        bNeedPrec = (ftp == efXTC || ftp == efGRO);
        int ftpin = fn2ftp(in_file);
        if (bVels)
        {
            /* check if velocities are possible in input and output files */
            bVels = (ftp == efTRR || ftp == efGRO ||
                     ftp == efG96 || ftp == efTNG)
                && (ftpin == efTRR || ftpin == efGRO ||
                    ftpin == efG96 || ftpin == efTNG || ftpin == efCPT);
        }
        if (bSeparate || bSplit)
        {
            outf_ext = std::strrchr(out_file, '.');
            if (outf_ext == nullptr)
            {
                gmx_fatal(FARGS, "Output file name '%s' does not contain a '.'", out_file);
            }
            outf_base = gmx_strdup(out_file);
            outf_base[outf_ext - out_file] = '\0';
        }

        bSubTraj = opt2bSet("-sub", NFILE, fnm);
        if (bSubTraj)
        {
            if ((ftp != efXTC) && (ftp != efTRR))
            {
                /* It seems likely that other trajectory file types
                 * could work here. */
                gmx_fatal(FARGS, "Can only use the sub option with output file types "
                          "xtc and trr");
            }
            clust = cluster_index(nullptr, opt2fn("-sub", NFILE, fnm));

            /* Check for number of files disabled, as FOPEN_MAX is not the correct
             * number to check for. In my linux box it is only 16.
             */
            if (0 && (clust->clust->nr > FOPEN_MAX-4))
            {
                gmx_fatal(FARGS, "Can not open enough (%d) files to write all the"
                          " trajectories.\ntry splitting the index file in %d parts.\n"
                          "FOPEN_MAX = %d",
                          clust->clust->nr, 1+clust->clust->nr/FOPEN_MAX, FOPEN_MAX);
            }
            gmx_warning("The -sub option could require as many open output files as there are\n"
                        "index groups in the file (%d). If you get I/O errors opening new files,\n"
                        "try reducing the number of index groups in the file, and perhaps\n"
                        "using trjconv -sub several times on different chunks of your index file.\n",
                        clust->clust->nr);

            snew(clust_status, clust->clust->nr);
            snew(clust_status_id, clust->clust->nr);
            snew(nfwritten, clust->clust->nr);
            for (i = 0; (i < clust->clust->nr); i++)
            {
                clust_status[i]    = nullptr;
                clust_status_id[i] = -1;
            }
            bSeparate = bSplit = FALSE;
        }
        /* skipping */
        if (skip_nr <= 0)
        {
        }

        mtop = read_mtop_for_tng(top_file, in_file, out_file);

        /* Determine whether to read a topology */
        bTPS = (ftp2bSet(efTPS, NFILE, fnm) ||
                bRmPBC || bReset || bPBCcomMol || bCluster ||
                (ftp == efGRO) || (ftp == efPDB) || bCONECT);

        /* Determine if when can read index groups */
        bIndex = (bIndex || bTPS);

        if (bTPS)
        {
            read_tps_conf(top_file, &top, &ePBC, &xp, nullptr, top_box,
                          bReset || bPBCcomRes);
            std::strncpy(top_title, *top.name, 255);
            top_title[255] = '\0';
            atoms          = &top.atoms;

            if (0 == top.mols.nr && (bCluster || bPBCcomMol))
            {
                gmx_fatal(FARGS, "Option -pbc %s requires a .tpr file for the -s option", pbc_opt[pbc_enum]);
            }

            /* top_title is only used for gro and pdb,
             * the header in such a file is top_title t= ...
             * to prevent a double t=, remove it from top_title
             */
            if ((charpt = std::strstr(top_title, " t= ")))
            {
                charpt[0] = '\0';
            }

            if (bCONECT)
            {
                gc = gmx_conect_generate(&top);
            }
            if (bRmPBC)
            {
                gpbc = gmx_rmpbc_init(&top.idef, ePBC, top.atoms.nr);
            }
        }

        /* get frame number index */
        frindex = nullptr;
        if (opt2bSet("-fr", NFILE, fnm))
        {
            printf("Select groups of frame number indices:\n");
            rd_index(opt2fn("-fr", NFILE, fnm), 1, &nrfri, (int **)&frindex, &frname);
            if (debug)
            {
                for (i = 0; i < nrfri; i++)
                {
                    fprintf(debug, "frindex[%4d]=%4d\n", i, frindex[i]);
                }
            }
        }

        /* get index groups etc. */
        if (bReset)
        {
            printf("Select group for %s fit\n",
                   bFit ? "least squares" : "translational");
            get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm),
                      1, &ifit, &ind_fit, &gn_fit);

            if (bFit)
            {
                if (ifit < 2)
                {
                    gmx_fatal(FARGS, "Need at least 2 atoms to fit!\n");
                }
                else if (ifit == 3)
                {
                    fprintf(stderr, "WARNING: fitting with only 2 atoms is not unique\n");
                }
            }
        }
        else if (bCluster)
        {
            printf("Select group for clustering\n");
            get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm),
                      1, &ifit, &ind_fit, &gn_fit);
        }

        if (bIndex)
        {
            if (bCenter)
            {
                printf("Select group for centering\n");
                get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm),
                          1, &ncent, &cindex, &grpnm);
            }
            printf("Select group for output\n");
            get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm),
                      1, &nout, &index, &grpnm);
        }
        else
        {
            /* no index file, so read natoms from TRX */
            if (!read_first_frame(oenv, &trxin, in_file, &fr, TRX_DONT_SKIP))
            {
                gmx_fatal(FARGS, "Could not read a frame from %s", in_file);
            }
            natoms = fr.natoms;
            close_trj(trxin);
            sfree(fr.x);
            snew(index, natoms);
            for (i = 0; i < natoms; i++)
            {
                index[i] = i;
            }
            nout = natoms;
            if (bCenter)
            {
                ncent  = nout;
                cindex = index;
            }
        }

        if (bReset)
        {
            snew(w_rls, atoms->nr);
            for (i = 0; (i < ifit); i++)
            {
                w_rls[ind_fit[i]] = atoms->atom[ind_fit[i]].m;
            }

            /* Restore reference structure and set to origin,
               store original location (to put structure back) */
            if (bRmPBC)
            {
                gmx_rmpbc(gpbc, top.atoms.nr, top_box, xp);
            }
            copy_rvec(xp[index[0]], x_shift);
            reset_x_ndim(nfitdim, ifit, ind_fit, atoms->nr, nullptr, xp, w_rls);
            rvec_dec(x_shift, xp[index[0]]);
        }
        else
        {
            clear_rvec(x_shift);
        }

        if (bDropUnder || bDropOver)
        {
            /* Read the .xvg file with the drop values */
            fprintf(stderr, "\nReading drop file ...");
            ndrop = read_xvg(opt2fn("-drop", NFILE, fnm), &dropval, &ncol);
            fprintf(stderr, " %d time points\n", ndrop);
            if (ndrop == 0 || ncol < 2)
            {
                gmx_fatal(FARGS, "Found no data points in %s",
                          opt2fn("-drop", NFILE, fnm));
            }
            drop0 = 0;
            drop1 = 0;
        }

        /* Make atoms struct for output in GRO or PDB files */
        if ((ftp == efGRO) || ((ftp == efG96) && bTPS) || (ftp == efPDB))
        {
            /* get memory for stuff to go in .pdb file, and initialize
             * the pdbinfo structure part if the input has it.
             */
            init_t_atoms(&useatoms, atoms->nr, atoms->havePdbInfo);
            sfree(useatoms.resinfo);
            useatoms.resinfo = atoms->resinfo;
            for (i = 0; (i < nout); i++)
            {
                useatoms.atomname[i] = atoms->atomname[index[i]];
                useatoms.atom[i]     = atoms->atom[index[i]];
                if (atoms->havePdbInfo)
                {
                    useatoms.pdbinfo[i]  = atoms->pdbinfo[index[i]];
                }
                useatoms.nres        = std::max(useatoms.nres, useatoms.atom[i].resind+1);
            }
            useatoms.nr = nout;
        }
        /* select what to read */
        if (ftp == efTRR)
        {
            flags = TRX_READ_X;
        }
        else
        {
            flags = TRX_NEED_X;
        }
        if (bVels)
        {
            flags = flags | TRX_READ_V;
        }
        if (bForce)
        {
            flags = flags | TRX_READ_F;
        }

        /* open trx file for reading */
        bHaveFirstFrame = read_first_frame(oenv, &trxin, in_file, &fr, flags);
        if (fr.bPrec)
        {
            fprintf(stderr, "\nPrecision of %s is %g (nm)\n", in_file, 1/fr.prec);
        }
        if (bNeedPrec)
        {
            if (bSetPrec || !fr.bPrec)
            {
                fprintf(stderr, "\nSetting output precision to %g (nm)\n", 1/prec);
            }
            else
            {
                fprintf(stderr, "Using output precision of %g (nm)\n", 1/prec);
            }
        }

        if (bHaveFirstFrame)
        {
            set_trxframe_ePBC(&fr, ePBC);

            natoms = fr.natoms;

            if (bSetTime)
            {
                tshift = tzero-fr.time;
            }
            else
            {
                tzero = fr.time;
            }

            bCopy = FALSE;
            if (bIndex)
            {
                /* check if index is meaningful */
                for (i = 0; i < nout; i++)
                {
                    if (index[i] >= natoms)
                    {
                        gmx_fatal(FARGS,
                                  "Index[%d] %d is larger than the number of atoms in the\n"
                                  "trajectory file (%d). There is a mismatch in the contents\n"
                                  "of your -f, -s and/or -n files.", i, index[i]+1, natoms);
                    }
                    bCopy = bCopy || (i != index[i]);
                }
            }

            /* open output for writing */
            std::strcpy(filemode, "w");
            switch (ftp)
            {
                case efTNG:
                    trjtools_gmx_prepare_tng_writing(out_file,
                                                     filemode[0],
                                                     trxin,
                                                     &trxout,
                                                     nullptr,
                                                     nout,
                                                     mtop,
                                                     index,
                                                     grpnm);
                    break;
                case efXTC:
                case efTRR:
                    out = nullptr;
                    if (!bSplit && !bSubTraj)
                    {
                        trxout = open_trx(out_file, filemode);
                    }
                    break;
                case efGRO:
                case efG96:
                case efPDB:
                    if (( !bSeparate && !bSplit ) && !bSubTraj)
                    {
                        out = gmx_ffopen(out_file, filemode);
                    }
                    break;
                default:
                    gmx_incons("Illegal output file format");
            }

            if (bCopy)
            {
                snew(xmem, nout);
                if (bVels)
                {
                    snew(vmem, nout);
                }
                if (bForce)
                {
                    snew(fmem, nout);
                }
            }

            /* Start the big loop over frames */
            file_nr  =  0;
            frame    =  0;
            outframe =  0;
            model_nr =  0;
            bDTset   = FALSE;

            /* Main loop over frames */
            do
            {
                if (!fr.bStep)
                {
                    /* set the step */
                    fr.step = newstep;
                    newstep++;
                }
                if (bSubTraj)
                {
                    /*if (frame >= clust->clust->nra)
                       gmx_fatal(FARGS,"There are more frames in the trajectory than in the cluster index file\n");*/
                    if (frame > clust->maxframe)
                    {
                        my_clust = -1;
                    }
                    else
                    {
                        my_clust = clust->inv_clust[frame];
                    }
                    if ((my_clust < 0) || (my_clust >= clust->clust->nr) ||
                        (my_clust == -1))
                    {
                        my_clust = -1;
                    }
                }

                if (bSetBox)
                {
                    /* generate new box */
                    if (fr.bBox == FALSE)
                    {
                        clear_mat(fr.box);
                    }
                    for (m = 0; m < DIM; m++)
                    {
                        if (newbox[m] >= 0)
                        {
                            fr.box[m][m] = newbox[m];
                        }
                        else
                        {
                            if (fr.bBox == FALSE)
                            {
                                gmx_fatal(FARGS, "Cannot preserve a box that does not exist.\n");
                            }
                        }
                    }
                }

                if (bTrans)
                {
                    for (i = 0; i < natoms; i++)
                    {
                        rvec_inc(fr.x[i], trans);
                    }
                }

                if (bTDump)
                {
                    /* determine timestep */
                    if (t0 == -1)
                    {
                        t0 = fr.time;
                    }
                    else
                    {
                        if (!bDTset)
                        {
                            dt     = fr.time-t0;
                            bDTset = TRUE;
                        }
                    }
                    /* This is not very elegant, as one can not dump a frame after
                     * a timestep with is more than twice as small as the first one. */
                    bDumpFrame = (fr.time > tdump-0.5*dt) && (fr.time <= tdump+0.5*dt);
                }
                else
                {
                    bDumpFrame = FALSE;
                }

                /* determine if an atom jumped across the box and reset it if so */
                if (bNoJump && (bTPS || frame != 0))
                {
                    for (d = 0; d < DIM; d++)
                    {
                        hbox[d] = 0.5*fr.box[d][d];
                    }
                    for (i = 0; i < natoms; i++)
                    {
                        if (bReset)
                        {
                            rvec_dec(fr.x[i], x_shift);
                        }
                        for (m = DIM-1; m >= 0; m--)
                        {
                            if (hbox[m] > 0)
                            {
                                while (fr.x[i][m]-xp[i][m] <= -hbox[m])
                                {
                                    for (d = 0; d <= m; d++)
                                    {
                                        fr.x[i][d] += fr.box[m][d];
                                    }
                                }
                                while (fr.x[i][m]-xp[i][m] > hbox[m])
                                {
                                    for (d = 0; d <= m; d++)
                                    {
                                        fr.x[i][d] -= fr.box[m][d];
                                    }
                                }
                            }
                        }
                    }
                }
                else if (bCluster)
                {
                    calc_pbc_cluster(ecenter, ifit, &top, ePBC, fr.x, ind_fit, fr.box);
                }

                if (bPFit)
                {
                    /* Now modify the coords according to the flags,
                       for normal fit, this is only done for output frames */
                    if (bRmPBC)
                    {
                        gmx_rmpbc_trxfr(gpbc, &fr);
                    }

                    reset_x_ndim(nfitdim, ifit, ind_fit, natoms, nullptr, fr.x, w_rls);
                    do_fit(natoms, w_rls, xp, fr.x);
                }

                /* store this set of coordinates for future use */
                if (bPFit || bNoJump)
                {
                    if (xp == nullptr)
                    {
                        snew(xp, natoms);
                    }
                    for (i = 0; (i < natoms); i++)
                    {
                        copy_rvec(fr.x[i], xp[i]);
                        rvec_inc(fr.x[i], x_shift);
                    }
                }

                if (frindex)
                {
                    /* see if we have a frame from the frame index group */
                    for (i = 0; i < nrfri && !bDumpFrame; i++)
                    {
                        bDumpFrame = frame == frindex[i];
                    }
                }
                if (debug && bDumpFrame)
                {
                    fprintf(debug, "dumping %d\n", frame);
                }

                bWriteFrame =
                    ( ( !bTDump && !frindex && frame % skip_nr == 0 ) || bDumpFrame );

                if (bWriteFrame && (bDropUnder || bDropOver))
                {
                    while (dropval[0][drop1] < fr.time && drop1+1 < ndrop)
                    {
                        drop0 = drop1;
                        drop1++;
                    }
                    if (std::abs(dropval[0][drop0] - fr.time)
                        < std::abs(dropval[0][drop1] - fr.time))
                    {
                        dropuse = drop0;
                    }
                    else
                    {
                        dropuse = drop1;
                    }
                    if ((bDropUnder && dropval[1][dropuse] < dropunder) ||
                        (bDropOver && dropval[1][dropuse] > dropover))
                    {
                        bWriteFrame = FALSE;
                    }
                }

                if (bWriteFrame)
                {
                    /* We should avoid modifying the input frame,
                     * but since here we don't have the output frame yet,
                     * we introduce a temporary output frame time variable.
                     */
                    real frout_time;

                    frout_time = fr.time;

                    /* calc new time */
                    if (bTimeStep)
                    {
                        frout_time = tzero + frame*timestep;
                    }
                    else
                    if (bSetTime)
                    {
                        frout_time += tshift;
                    }

                    if (bTDump)
                    {
                        fprintf(stderr, "\nDumping frame at t= %g %s\n",
                                output_env_conv_time(oenv, frout_time), output_env_get_time_unit(oenv).c_str());
                    }

                    /* check for writing at each delta_t */
                    bDoIt = (delta_t == 0);
                    if (!bDoIt)
                    {
                        if (!bRound)
                        {
                            bDoIt = bRmod(frout_time, tzero, delta_t);
                        }
                        else
                        {
                            /* round() is not C89 compatible, so we do this:  */
                            bDoIt = bRmod(std::floor(frout_time+0.5), std::floor(tzero+0.5),
                                          std::floor(delta_t+0.5));
                        }
                    }

                    if (bDoIt || bTDump)
                    {
                        /* print sometimes */
                        if ( ((outframe % SKIP) == 0) || (outframe < SKIP) )
                        {
                            fprintf(stderr, " ->  frame %6d time %8.3f      \r",
                                    outframe, output_env_conv_time(oenv, frout_time));
                            fflush(stderr);
                        }

                        if (!bPFit)
                        {
                            /* Now modify the coords according to the flags,
                               for PFit we did this already! */

                            if (bRmPBC)
                            {
                                gmx_rmpbc_trxfr(gpbc, &fr);
                            }

                            if (bReset)
                            {
                                reset_x_ndim(nfitdim, ifit, ind_fit, natoms, nullptr, fr.x, w_rls);
                                if (bFit)
                                {
                                    do_fit_ndim(nfitdim, natoms, w_rls, xp, fr.x);
                                }
                                if (!bCenter)
                                {
                                    for (i = 0; i < natoms; i++)
                                    {
                                        rvec_inc(fr.x[i], x_shift);
                                    }
                                }
                            }

                            if (bCenter)
                            {
                                center_x(ecenter, fr.x, fr.box, natoms, ncent, cindex);
                            }
                        }

                        if (bPBCcomAtom)
                        {
                            switch (unitcell_enum)
                            {
                                case euRect:
                                    put_atoms_in_box(ePBC, fr.box, natoms, fr.x);
                                    break;
                                case euTric:
                                    put_atoms_in_triclinic_unitcell(ecenter, fr.box, natoms, fr.x);
                                    break;
                                case euCompact:
                                    put_atoms_in_compact_unitcell(ePBC, ecenter, fr.box,
                                                                  natoms, fr.x);
                                    break;
                            }
                        }
                        if (bPBCcomRes)
                        {
                            put_residue_com_in_box(unitcell_enum, ecenter,
                                                   natoms, atoms->atom, ePBC, fr.box, fr.x);
                        }
                        if (bPBCcomMol)
                        {
                            put_molecule_com_in_box(unitcell_enum, ecenter,
                                                    &top.mols,
                                                    natoms, atoms->atom, ePBC, fr.box, fr.x);
                        }
                        /* Copy the input trxframe struct to the output trxframe struct */
                        frout        = fr;
                        frout.time   = frout_time;
                        frout.bV     = (frout.bV && bVels);
                        frout.bF     = (frout.bF && bForce);
                        frout.natoms = nout;
                        if (bNeedPrec && (bSetPrec || !fr.bPrec))
                        {
                            frout.bPrec = TRUE;
                            frout.prec  = prec;
                        }
                        if (bCopy)
                        {
                            frout.x = xmem;
                            if (frout.bV)
                            {
                                frout.v = vmem;
                            }
                            if (frout.bF)
                            {
                                frout.f = fmem;
                            }
                            for (i = 0; i < nout; i++)
                            {
                                copy_rvec(fr.x[index[i]], frout.x[i]);
                                if (frout.bV)
                                {
                                    copy_rvec(fr.v[index[i]], frout.v[i]);
                                }
                                if (frout.bF)
                                {
                                    copy_rvec(fr.f[index[i]], frout.f[i]);
                                }
                            }
                        }

                        if (opt2parg_bSet("-shift", NPA, pa))
                        {
                            for (i = 0; i < nout; i++)
                            {
                                for (d = 0; d < DIM; d++)
                                {
                                    frout.x[i][d] += outframe*shift[d];
                                }
                            }
                        }

                        if (!bRound)
                        {
                            bSplitHere = bSplit && bRmod(frout.time, tzero, split_t);
                        }
                        else
                        {
                            /* round() is not C89 compatible, so we do this: */
                            bSplitHere = bSplit && bRmod(std::floor(frout.time+0.5),
                                                         std::floor(tzero+0.5),
                                                         std::floor(split_t+0.5));
                        }
                        if (bSeparate || bSplitHere)
                        {
                            mk_filenm(outf_base, ftp2ext(ftp), nzero, file_nr, out_file2);
                        }

                        switch (ftp)
                        {
                            case efTNG:
                                write_tng_frame(trxout, &frout);
                                // TODO when trjconv behaves better: work how to read and write lambda
                                break;
                            case efTRR:
                            case efXTC:
                                if (bSplitHere)
                                {
                                    if (trxout)
                                    {
                                        close_trx(trxout);
                                    }
                                    trxout = open_trx(out_file2, filemode);
                                }
                                if (bSubTraj)
                                {
                                    if (my_clust != -1)
                                    {
                                        char buf[STRLEN];
                                        if (clust_status_id[my_clust] == -1)
                                        {
                                            sprintf(buf, "%s.%s", clust->grpname[my_clust], ftp2ext(ftp));
                                            clust_status[my_clust]    = open_trx(buf, "w");
                                            clust_status_id[my_clust] = 1;
                                            ntrxopen++;
                                        }
                                        else if (clust_status_id[my_clust] == -2)
                                        {
                                            gmx_fatal(FARGS, "File %s.xtc should still be open (%d open .xtc files)\n" "in order to write frame %d. my_clust = %d",
                                                      clust->grpname[my_clust], ntrxopen, frame,
                                                      my_clust);
                                        }
                                        write_trxframe(clust_status[my_clust], &frout, gc);
                                        nfwritten[my_clust]++;
                                        if (nfwritten[my_clust] ==
                                            (clust->clust->index[my_clust+1]-
                                             clust->clust->index[my_clust]))
                                        {
                                            close_trx(clust_status[my_clust]);
                                            clust_status[my_clust]    = nullptr;
                                            clust_status_id[my_clust] = -2;
                                            ntrxopen--;
                                            if (ntrxopen < 0)
                                            {
                                                gmx_fatal(FARGS, "Less than zero open .xtc files!");
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    write_trxframe(trxout, &frout, gc);
                                }
                                break;
                            case efGRO:
                            case efG96:
                            case efPDB:
                                sprintf(title, "Generated by trjconv : %s t= %9.5f",
                                        top_title, frout.time);
                                if (bSeparate || bSplitHere)
                                {
                                    out = gmx_ffopen(out_file2, "w");
                                }
                                switch (ftp)
                                {
                                    case efGRO:
                                        write_hconf_p(out, title, &useatoms,
                                                      frout.x, frout.bV ? frout.v : nullptr, frout.box);
                                        break;
                                    case efPDB:
                                        fprintf(out, "REMARK    GENERATED BY TRJCONV\n");
                                        sprintf(title, "%s t= %9.5f", top_title, frout.time);
                                        /* if reading from pdb, we want to keep the original
                                           model numbering else we write the output frame
                                           number plus one, because model 0 is not allowed in pdb */
                                        if (ftpin == efPDB && fr.bStep && fr.step > model_nr)
                                        {
                                            model_nr = fr.step;
                                        }
                                        else
                                        {
                                            model_nr++;
                                        }
                                        write_pdbfile(out, title, &useatoms, frout.x,
                                                      frout.ePBC, frout.box, ' ', model_nr, gc, TRUE);
                                        break;
                                    case efG96:
                                        frout.title = title;
                                        if (bSeparate || bTDump)
                                        {
                                            frout.bTitle = TRUE;
                                            if (bTPS)
                                            {
                                                frout.bAtoms = TRUE;
                                            }
                                            frout.atoms  = &useatoms;
                                            frout.bStep  = FALSE;
                                            frout.bTime  = FALSE;
                                        }
                                        else
                                        {
                                            frout.bTitle = (outframe == 0);
                                            frout.bAtoms = FALSE;
                                            frout.bStep  = TRUE;
                                            frout.bTime  = TRUE;
                                        }
                                        write_g96_conf(out, &frout, -1, nullptr);
                                }
                                if (bSeparate || bSplitHere)
                                {
                                    gmx_ffclose(out);
                                    out = nullptr;
                                }
                                break;
                            default:
                                gmx_fatal(FARGS, "DHE, ftp=%d\n", ftp);
                        }
                        if (bSeparate || bSplitHere)
                        {
                            file_nr++;
                        }

                        /* execute command */
                        if (bExec)
                        {
                            char c[255];
                            sprintf(c, "%s  %d", exec_command, file_nr-1);
                            /*fprintf(stderr,"Executing '%s'\n",c);*/
                            if (0 != system(c))
                            {
                                gmx_fatal(FARGS, "Error executing command: %s", c);
                            }
                        }
                        outframe++;
                    }
                }
                frame++;
                bHaveNextFrame = read_next_frame(oenv, trxin, &fr);
            }
            while (!(bTDump && bDumpFrame) && bHaveNextFrame);
        }

        if (!bHaveFirstFrame || (bTDump && !bDumpFrame))
        {
            fprintf(stderr, "\nWARNING no output, "
                    "last frame read at t=%g\n", fr.time);
        }
        fprintf(stderr, "\n");

        close_trj(trxin);
        sfree(outf_base);

        if (bRmPBC)
        {
            gmx_rmpbc_done(gpbc);
        }

        if (trxout)
        {
            close_trx(trxout);
        }
        else if (out != nullptr)
        {
            gmx_ffclose(out);
        }
        if (bSubTraj)
        {
            for (i = 0; (i < clust->clust->nr); i++)
            {
                if (clust_status_id[i] >= 0)
                {
                    close_trx(clust_status[i]);
                }
            }
        }
    }

    sfree(mtop);

    do_view(oenv, out_file, nullptr);

    return 0;
}
