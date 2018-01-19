/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/utility/arrayref.h"
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
    rvec    com, shift, box_center;
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
        gmx::RVec newCom;
        copy_rvec(com, newCom);
        auto      newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
        switch (unitcell_enum)
        {
            case euRect:
                put_atoms_in_box(ePBC, box, newComArrayRef);
                break;
            case euTric:
                put_atoms_in_triclinic_unitcell(ecenter, box, newComArrayRef);
                break;
            case euCompact:
                put_atoms_in_compact_unitcell(ePBC, ecenter, box, newComArrayRef);
                break;
        }
        rvec_sub(newCom, com, shift);
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
    rvec             box_center, com, shift;
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
            gmx::RVec newCom;
            copy_rvec(com, newCom);
            auto      newComArrayRef = gmx::arrayRefFromArray(&newCom, 1);
            switch (unitcell_enum)
            {
                case euRect:
                    put_atoms_in_box(ePBC, box, newComArrayRef);
                    break;
                case euTric:
                    put_atoms_in_triclinic_unitcell(ecenter, box, newComArrayRef);
                    break;
                case euCompact:
                    put_atoms_in_compact_unitcell(ePBC, ecenter, box, newComArrayRef);
                    break;
            }
            rvec_sub(newCom, com, shift);
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

int gmx_trjpbc(int argc, char *argv[])
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

    static gmx_bool  bVels = TRUE, bForce = FALSE, bCONECT = FALSE;
    static gmx_bool  bCenter       = FALSE;
    static int       skip_nr       = 1, ndec = 3, nzero = 0;
    static real      delta_t = 0;
    static rvec      newbox        = {0, 0, 0}, shift = {0, 0, 0}, trans = {0, 0, 0};
    static char     *exec_command  = nullptr;
    bool    bdump;

    t_pargs
        pa[] =
    {
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
          "decimal places. Ignored for all other formats" },
        { "-vel", FALSE, etBOOL,
          { &bVels }, "Read and write velocities if possible" },
        { "-force", FALSE, etBOOL,
          { &bForce }, "Read and write forces if possible" },
        { "-exec", FALSE, etSTR,
          { &exec_command },
          "Execute command for every output frame with the "
          "frame number as argument" },
        { "-nzero", FALSE, etINT,
          { &nzero },
          "If the -sep flag is set, use these many digits "
          "for the file numbers and prepend zeros as needed" },
        { "-conect", FALSE, etBOOL,
          { &bCONECT },
          "Add conect records when writing [REF].pdb[ref] files. Useful "
          "for visualization of non-standard molecules, e.g. "
          "coarse grained ones. Ignored for other output types" },
        { "-dump", FALSE, etBOOL,
          { &bdump }, "Not doing anything, just used for backwards compatability" },
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
    int               m, i, d, frame, outframe, natoms, nout, ncent, newstep = 0;
#define SKIP 10
    t_topology        top;
    gmx_mtop_t       *mtop  = nullptr;
    gmx_conect        gc    = nullptr;
    int               ePBC  = -1;
    t_atoms          *atoms = nullptr, useatoms;
    matrix            top_box;
    int              *index = nullptr, *cindex = nullptr;
    char             *grpnm = nullptr;
    int               ifit;
    int              *ind_fit;
    char             *gn_fit;
    real              prec;
    gmx_bool          bFit, bPFit, bReset;
    int               nfitdim;
    gmx_rmpbc_t       gpbc = nullptr;
    gmx_bool          bRmPBC, bPBCWhole, bPBCcomRes, bPBCcomMol, bPBCcomAtom, bPBC, bNoJump, bCluster;
    gmx_bool          bCopy, bDoIt, bIndex, bTDump, bTPS = FALSE;
    gmx_bool          bExec, bSetPrec, bNeedPrec;
    gmx_bool          bHaveFirstFrame, bHaveNextFrame, bSetBox, bSetUR;
    gmx_bool          bTrans = FALSE;
    gmx_bool          bWriteFrame;
    const char       *top_file, *in_file, *out_file = nullptr;
    char             *charpt;
    char             *outf_base = nullptr;
    char              top_title[256], title[256], timestr[32], stepstr[32], filemode[5];
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

    /* Check command line */
    in_file = opt2fn("-f", NFILE, fnm);


        /* mark active cmdline options */
        bSetBox    = opt2parg_bSet("-box", NPA, pa);
        bSetPrec   = opt2parg_bSet("-ndec", NPA, pa);
        bSetUR     = opt2parg_bSet("-ur", NPA, pa);
        bExec      = opt2parg_bSet("-exec", NPA, pa);
        bTDump     = opt2parg_bSet("-dump", NPA, pa);
        bTrans     = opt2parg_bSet("-trans", NPA, pa);

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


        /* Output type is equal to input type 
	 * ignores fileypes given by user, but gives warning because this is new
	 */
        out_file = opt2fn("-o", NFILE, fnm);
        int ftp  = fn2ftp(out_file);
        int ftpin = fn2ftp(in_file);
	if (ftpin != ftp)
{
	ftp = ftpin;

        // hack the fnm structure to update the output file name
        sprintf(*fnm[1].fns, "%s.%s", *fnm[1].fns,ftp2ext(ftpin));

        out_file = opt2fn("-o", NFILE, fnm);
	fprintf(stderr, "Note! You defined the output file type, but trjpbc always\n"
			"keeps the input file format. Please use trjconv to change file formats!\n");
}
        fprintf(stderr, "Will write %s: %s\n", ftp2ext(ftp), ftp2desc(ftp));
        bNeedPrec = (ftp == efXTC || ftp == efGRO);

        if (bVels)
        {
            /* check if velocities are possible in input and output files */
            bVels = (ftp == efTRR || ftp == efGRO ||
                     ftp == efG96 || ftp == efTNG)
                && (ftpin == efTRR || ftpin == efGRO ||
                    ftpin == efG96 || ftpin == efTNG || ftpin == efCPT);
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
             * the header in such a file is top_title, followed by
             * t= ... and/or step= ...
             * to prevent double t= or step=, remove it from top_title.
             * From GROMACS-2018 we only write t/step when the frame actually
             * has a valid time/step, so we need to check for both separately.
             */
            if ((charpt = std::strstr(top_title, " t= ")))
            {
                charpt[0] = '\0';
            }
            if ((charpt = std::strstr(top_title, " step= ")))
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
            close_trx(trxin);
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
                    trxout = open_trx(out_file, filemode);
                    break;
                case efGRO:
                case efG96:
                case efPDB:
                    out = gmx_ffopen(out_file, filemode);
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

            /* Main loop over frames */
            do
            {
                if (!fr.bStep)
                {
                    /* set the step */
                    fr.step = newstep;
                    newstep++;
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

                bWriteFrame =
                    ( ( frame % skip_nr == 0 ) );

                if (bWriteFrame)
                {
                    /* We should avoid modifying the input frame,
                     * but since here we don't have the output frame yet,
                     * we introduce a temporary output frame time variable.
                     */
                    real frout_time;

                    frout_time = fr.time;

                    /* check for writing at each delta_t */
                    bDoIt = (delta_t == 0);

                    if (bDoIt || bTDump)
                    {

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

                        auto positionsArrayRef = gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(fr.x), natoms);
                        if (bPBCcomAtom)
                        {
                            switch (unitcell_enum)
                            {
                                case euRect:
                                    put_atoms_in_box(ePBC, fr.box, positionsArrayRef);
                                    break;
                                case euTric:
                                    put_atoms_in_triclinic_unitcell(ecenter, fr.box, positionsArrayRef);
                                    break;
                                case euCompact:
                                    put_atoms_in_compact_unitcell(ePBC, ecenter, fr.box,
                                                                  positionsArrayRef);
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

                        switch (ftp)
                        {
                            case efTNG:
                                write_tng_frame(trxout, &frout);
                                // TODO when trjconv behaves better: work how to read and write lambda
                                break;
                            case efTRR:
                            case efXTC:

                                {
                                    write_trxframe(trxout, &frout, gc);
                                }
                                break;
                            case efGRO:
                            case efG96:
                            case efPDB:
                                // Only add a generator statement if title is empty,
                                // to avoid multiple generated-by statements from various programs
                                if (std::strlen(top_title) == 0)
                                {
                                    sprintf(top_title, "Generated by trjconv");
                                }
                                if (frout.bTime)
                                {
                                    sprintf(timestr, " t= %9.5f", frout.time);
                                }
                                else
                                {
                                    std::strcpy(timestr, "");
                                }
                                if (frout.bStep)
                                {
                                    sprintf(stepstr, " step= %" GMX_PRId64, frout.step);
                                }
                                else
                                {
                                    std::strcpy(stepstr, "");
                                }
                                snprintf(title, 256, "%s%s%s", top_title, timestr, stepstr);
                                break;
                            default:
                                gmx_fatal(FARGS, "DHE, ftp=%d\n", ftp);
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
            while (bHaveNextFrame);
        }

        if (!bHaveFirstFrame)
        {
            fprintf(stderr, "\nWARNING no output, "
                    "last frame read at t=%g\n", fr.time);
        }
        fprintf(stderr, "\n");

        close_trx(trxin);
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
    }

    sfree(mtop);
    done_top(&top);
    sfree(xp);
    sfree(xmem);
    sfree(vmem);
    sfree(fmem);
    sfree(grpnm);
    sfree(index);
    sfree(cindex);
    done_filenms(NFILE, fnm);
    done_frame(&fr);

    do_view(oenv, out_file, nullptr);

    output_env_done(oenv);
    return 0;
}
