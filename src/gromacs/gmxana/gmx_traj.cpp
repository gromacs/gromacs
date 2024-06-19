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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/nrjac.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;

static void
low_print_data(FILE* fp, real time, rvec x[], int n, const int* index, const gmx_bool bDim[], const char* sffmt)
{
    int i, ii, d;

    fprintf(fp, " %g", time);
    for (i = 0; i < n; i++)
    {
        if (index != nullptr)
        {
            ii = index[i];
        }
        else
        {
            ii = i;
        }
        for (d = 0; d < DIM; d++)
        {
            if (bDim[d])
            {
                fprintf(fp, sffmt, x[ii][d]);
            }
        }
        if (bDim[DIM])
        {
            fprintf(fp, sffmt, norm(x[ii]));
        }
    }
    fprintf(fp, "\n");
}

static void average_data(rvec x[], rvec xav[], const real* mass, int ngrps, const int isize[], int** index)
{
    int    g, i, ind, d;
    real   m;
    rvec   tmp;
    double sum[DIM], mtot;

    for (g = 0; g < ngrps; g++)
    {
        clear_dvec(sum);
        clear_rvec(xav[g]);
        mtot = 0;
        for (i = 0; i < isize[g]; i++)
        {
            ind = index[g][i];
            if (mass != nullptr)
            {
                m = mass[ind];
                svmul(m, x[ind], tmp);
                for (d = 0; d < DIM; d++)
                {
                    sum[d] += tmp[d];
                }
                mtot += m;
            }
            else
            {
                for (d = 0; d < DIM; d++)
                {
                    sum[d] += x[ind][d];
                }
            }
        }
        if (mass != nullptr)
        {
            for (d = 0; d < DIM; d++)
            {
                xav[g][d] = sum[d] / mtot;
            }
        }
        else
        {
            /* mass=NULL, so these are forces: sum only (do not average) */
            for (d = 0; d < DIM; d++)
            {
                xav[g][d] = sum[d];
            }
        }
    }
}

static void print_data(FILE*       fp,
                       real        time,
                       rvec        x[],
                       real*       mass,
                       gmx_bool    bCom,
                       int         ngrps,
                       int         isize[],
                       int**       index,
                       gmx_bool    bDim[],
                       const char* sffmt)
{
    static std::vector<gmx::RVec> xav;

    if (bCom)
    {
        if (xav.empty())
        {
            xav.resize(ngrps);
        }
        average_data(x, as_rvec_array(xav.data()), mass, ngrps, isize, index);
        low_print_data(fp, time, as_rvec_array(xav.data()), ngrps, nullptr, bDim, sffmt);
    }
    else
    {
        low_print_data(fp, time, x, isize[0], index[0], bDim, sffmt);
    }
}

static void write_trx_x(t_trxstatus*      status,
                        const t_trxframe* fr,
                        real*             mass,
                        gmx_bool          bCom,
                        int               ngrps,
                        int               isize[],
                        int**             index)
{
    static std::vector<gmx::RVec> xav;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static t_atoms* atoms = nullptr;
    t_trxframe      fr_av;
    int             i;

    if (bCom)
    {
        if (xav.empty())
        {
            xav.resize(ngrps);
            snew(atoms, 1);
            *atoms = *fr->atoms;
            snew(atoms->atom, ngrps);
            atoms->nr = ngrps;
            /* Note that atom and residue names will be the ones
             * of the first atom in each group.
             */
            for (i = 0; i < ngrps; i++)
            {
                atoms->atom[i]     = fr->atoms->atom[index[i][0]];
                atoms->atomname[i] = fr->atoms->atomname[index[i][0]];
            }
        }
        average_data(fr->x, as_rvec_array(xav.data()), mass, ngrps, isize, index);
        fr_av        = *fr;
        fr_av.natoms = ngrps;
        fr_av.atoms  = atoms;
        fr_av.x      = as_rvec_array(xav.data());
        write_trxframe(status, &fr_av, nullptr);
    }
    else
    {
        write_trxframe_indexed(status, fr, isize[0], index[0], nullptr);
    }
}

static void make_legend(FILE*                   fp,
                        int                     ngrps,
                        int                     isize,
                        int                     index[],
                        char**                  name,
                        gmx_bool                bCom,
                        gmx_bool                bMol,
                        const gmx_bool          bDim[],
                        const gmx_output_env_t* oenv)
{
    std::vector<std::string> leg;
    const char*              dimtxt[] = { " X", " Y", " Z", "" };

    const int n = bCom ? ngrps : isize;

    for (int i = 0; i < n; i++)
    {
        for (int d = 0; d <= DIM; d++)
        {
            if (bDim[d])
            {
                if (bMol)
                {
                    leg.emplace_back(gmx::formatString("mol %d%s", index[i] + 1, dimtxt[d]));
                }
                else if (bCom)
                {
                    leg.emplace_back(gmx::formatString("%s%s", name[i], dimtxt[d]));
                }
                else
                {
                    leg.emplace_back(gmx::formatString("atom %d%s", index[i] + 1, dimtxt[d]));
                }
            }
        }
    }
    xvgrLegend(fp, leg, oenv);
}

static real ekrot(rvec x[], rvec v[], const real mass[], int isize, const int index[])
{
    real   TCM[DIM][DIM], L[DIM][DIM];
    double tm, m0, lxx, lxy, lxz, lyy, lyz, lzz, ekrot;
    rvec   a0, ocm;
    dvec   dx, b0;
    dvec   xcm, vcm, acm;
    int    i, j, m, n;

    clear_dvec(xcm);
    clear_dvec(vcm);
    clear_dvec(acm);
    tm = 0.0;
    for (i = 0; i < isize; i++)
    {
        j  = index[i];
        m0 = mass[j];
        tm += m0;
        cprod(x[j], v[j], a0);
        for (m = 0; (m < DIM); m++)
        {
            xcm[m] += m0 * x[j][m]; /* c.o.m. position */
            vcm[m] += m0 * v[j][m]; /* c.o.m. velocity */
            acm[m] += m0 * a0[m];   /* rotational velocity around c.o.m. */
        }
    }
    dcprod(xcm, vcm, b0);
    for (m = 0; (m < DIM); m++)
    {
        xcm[m] /= tm;
        vcm[m] /= tm;
        acm[m] -= b0[m] / tm;
    }

    lxx = lxy = lxz = lyy = lyz = lzz = 0.0;
    for (i = 0; i < isize; i++)
    {
        j  = index[i];
        m0 = mass[j];
        for (m = 0; m < DIM; m++)
        {
            dx[m] = x[j][m] - xcm[m];
        }
        lxx += dx[XX] * dx[XX] * m0;
        lxy += dx[XX] * dx[YY] * m0;
        lxz += dx[XX] * dx[ZZ] * m0;
        lyy += dx[YY] * dx[YY] * m0;
        lyz += dx[YY] * dx[ZZ] * m0;
        lzz += dx[ZZ] * dx[ZZ] * m0;
    }

    L[XX][XX] = lyy + lzz;
    L[YY][XX] = -lxy;
    L[ZZ][XX] = -lxz;
    L[XX][YY] = -lxy;
    L[YY][YY] = lxx + lzz;
    L[ZZ][YY] = -lyz;
    L[XX][ZZ] = -lxz;
    L[YY][ZZ] = -lyz;
    L[ZZ][ZZ] = lxx + lyy;

    m_inv_gen(&L[0][0], DIM, &TCM[0][0]);

    /* Compute omega (hoeksnelheid) */
    clear_rvec(ocm);
    ekrot = 0;
    for (m = 0; m < DIM; m++)
    {
        for (n = 0; n < DIM; n++)
        {
            ocm[m] += TCM[m][n] * acm[n];
        }
        ekrot += 0.5 * ocm[m] * acm[m];
    }

    return ekrot;
}

static real ektrans(rvec v[], const real mass[], int isize, const int index[])
{
    dvec   mvcom;
    double mtot;
    int    i, j, d;

    clear_dvec(mvcom);
    mtot = 0;
    for (i = 0; i < isize; i++)
    {
        j = index[i];
        for (d = 0; d < DIM; d++)
        {
            mvcom[d] += mass[j] * v[j][d];
        }
        mtot += mass[j];
    }

    return dnorm2(mvcom) / (mtot * 2);
}

static real temp(rvec v[], const real mass[], int isize, const int index[])
{
    double ekin2;
    int    i, j;

    ekin2 = 0;
    for (i = 0; i < isize; i++)
    {
        j = index[i];
        ekin2 += mass[j] * norm2(v[j]);
    }

    return ekin2 / (3 * isize * gmx::c_boltz);
}

static void remove_jump(matrix box, int natoms, rvec xp[], rvec x[])
{
    rvec hbox;
    int  d, i, m;

    for (d = 0; d < DIM; d++)
    {
        hbox[d] = 0.5 * box[d][d];
    }
    for (i = 0; i < natoms; i++)
    {
        for (m = DIM - 1; m >= 0; m--)
        {
            while (x[i][m] - xp[i][m] <= -hbox[m])
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] += box[m][d];
                }
            }
            while (x[i][m] - xp[i][m] > hbox[m])
            {
                for (d = 0; d <= m; d++)
                {
                    x[i][d] -= box[m][d];
                }
            }
        }
    }
}

static void write_pdb_bfac(const char*             fname,
                           const char*             xname,
                           const char*             title,
                           t_atoms*                atoms,
                           PbcType                 pbcType,
                           matrix                  box,
                           int                     isize,
                           int*                    index,
                           int                     nfr_x,
                           rvec*                   x,
                           int                     nfr_v,
                           rvec*                   sum,
                           const gmx_bool          bDim[],
                           real                    scale_factor,
                           const gmx_output_env_t* oenv)
{
    FILE* fp;
    real  max, len2, scale;
    int   maxi;
    int   i, m, onedim;

    if ((nfr_x == 0) || (nfr_v == 0))
    {
        fprintf(stderr, "No frames found for %s, will not write %s\n", title, fname);
    }
    else
    {
        fprintf(stderr, "Used %d frames for %s\n", nfr_x, "coordinates");
        fprintf(stderr, "Used %d frames for %s\n", nfr_v, title);
        onedim = -1;
        if (!bDim[DIM])
        {
            m = 0;
            for (i = 0; i < DIM; i++)
            {
                if (bDim[i])
                {
                    onedim = i;
                    m++;
                }
            }
            if (m != 1)
            {
                onedim = -1;
            }
        }
        scale = 1.0 / nfr_v;
        for (i = 0; i < isize; i++)
        {
            svmul(scale, sum[index[i]], sum[index[i]]);
        }

        fp = xvgropen(xname, title, "Atom", "Spatial component", oenv);
        for (i = 0; i < isize; i++)
        {
            fprintf(fp,
                    "%-5d  %10.3f  %10.3f  %10.3f\n",
                    1 + i,
                    sum[index[i]][XX],
                    sum[index[i]][YY],
                    sum[index[i]][ZZ]);
        }
        xvgrclose(fp);
        max  = 0;
        maxi = 0;
        for (i = 0; i < isize; i++)
        {
            len2 = 0;
            for (m = 0; m < DIM; m++)
            {
                if (bDim[m] || bDim[DIM])
                {
                    len2 += gmx::square(sum[index[i]][m]);
                }
            }
            if (len2 > max)
            {
                max  = len2;
                maxi = index[i];
            }
        }
        if (scale_factor != 0)
        {
            scale = scale_factor;
        }
        else
        {
            if (max == 0)
            {
                scale = 1;
            }
            else
            {
                scale = 10.0 / std::sqrt(max);
            }
        }

        printf("Maximum %s is %g on atom %d %s, res. %s %d\n",
               title,
               std::sqrt(max),
               maxi + 1,
               *(atoms->atomname[maxi]),
               *(atoms->resinfo[atoms->atom[maxi].resind].name),
               atoms->resinfo[atoms->atom[maxi].resind].nr);

        if (atoms->pdbinfo == nullptr)
        {
            snew(atoms->pdbinfo, atoms->nr);
        }
        atoms->havePdbInfo = TRUE;

        if (onedim == -1)
        {
            for (i = 0; i < isize; i++)
            {
                len2 = 0;
                for (m = 0; m < DIM; m++)
                {
                    if (bDim[m] || bDim[DIM])
                    {
                        len2 += gmx::square(sum[index[i]][m]);
                    }
                }
                atoms->pdbinfo[index[i]].bfac = std::sqrt(len2) * scale;
            }
        }
        else
        {
            for (i = 0; i < isize; i++)
            {
                atoms->pdbinfo[index[i]].bfac = sum[index[i]][onedim] * scale;
            }
        }
        write_sto_conf_indexed(fname, title, atoms, x, nullptr, pbcType, box, isize, index);
    }
}

static void update_histo(int gnx, const int index[], rvec v[], int* nhisto, int** histo, real binwidth)
{
    int  i, m, in, nnn;
    real vn, vnmax;

    if (*histo == nullptr)
    {
        vnmax = 0;
        for (i = 0; (i < gnx); i++)
        {
            vn    = norm(v[index[i]]);
            vnmax = std::max(vn, vnmax);
        }
        vnmax *= 2;
        *nhisto = static_cast<int>(1 + (vnmax / binwidth));
        snew(*histo, *nhisto);
    }
    for (i = 0; (i < gnx); i++)
    {
        vn = norm(v[index[i]]);
        in = static_cast<int>(vn / binwidth);
        if (in >= *nhisto)
        {
            nnn = in + 100;
            fprintf(stderr, "Extending histogram from %d to %d\n", *nhisto, nnn);

            srenew(*histo, nnn);
            for (m = *nhisto; (m < nnn); m++)
            {
                (*histo)[m] = 0;
            }
            *nhisto = nnn;
        }
        (*histo)[in]++;
    }
}

static void print_histo(const char* fn, int nhisto, int histo[], real binwidth, const gmx_output_env_t* oenv)
{
    FILE* fp;
    int   i;

    fp = xvgropen(fn, "Velocity distribution", "V (nm/ps)", "arbitrary units", oenv);
    for (i = 0; (i < nhisto); i++)
    {
        fprintf(fp, "%10.3e  %10d\n", i * binwidth, histo[i]);
    }
    xvgrclose(fp);
}

int gmx_traj(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] plots coordinates, velocities, forces and/or the box.",
        "With [TT]-com[tt] the coordinates, velocities and forces are",
        "calculated for the center of mass of each group.",
        "When [TT]-mol[tt] is set, the numbers in the index file are",
        "interpreted as molecule numbers and the same procedure as with",
        "[TT]-com[tt] is used for each molecule.[PAR]",
        "Option [TT]-ot[tt] plots the temperature of each group,",
        "provided velocities are present in the trajectory file.",
        "No corrections are made for constrained degrees of freedom!",
        "This implies [TT]-com[tt].[PAR]",
        "Options [TT]-ekt[tt] and [TT]-ekr[tt] plot the translational and",
        "rotational kinetic energy of each group,",
        "provided velocities are present in the trajectory file.",
        "This implies [TT]-com[tt].[PAR]",
        "Options [TT]-cv[tt] and [TT]-cf[tt] write the average velocities",
        "and average forces as temperature factors to a [REF].pdb[ref] file with",
        "the average coordinates or the coordinates at [TT]-ctime[tt].",
        "The temperature factors are scaled such that the maximum is 10.",
        "The scaling can be changed with the option [TT]-scale[tt].",
        "To get the velocities or forces of one",
        "frame set both [TT]-b[tt] and [TT]-e[tt] to the time of",
        "desired frame. When averaging over frames you might need to use",
        "the [TT]-nojump[tt] option to obtain the correct average coordinates.",
        "If you select either of these option the average force and velocity",
        "for each atom are written to an [REF].xvg[ref] file as well",
        "(specified with [TT]-av[tt] or [TT]-af[tt]).[PAR]",
        "Option [TT]-vd[tt] computes a velocity distribution, i.e. the",
        "norm of the vector is plotted. In addition in the same graph",
        "the kinetic energy distribution is given.",
        "",
        "See [gmx-trajectory] for plotting similar data for selections."
    };
    static gmx_bool bMol = FALSE, bCom = FALSE, bPBC = TRUE, bNoJump = FALSE;
    static gmx_bool bX = TRUE, bY = TRUE, bZ = TRUE, bNorm = FALSE, bFP = FALSE;
    static int      ngroups = 1;
    static real     ctime = -1, scale = 0, binwidth = 1;
    t_pargs         pa[] = {
        { "-com", FALSE, etBOOL, { &bCom }, "Plot data for the com of each group" },
        { "-pbc", FALSE, etBOOL, { &bPBC }, "Make molecules whole for COM" },
        { "-mol",
          FALSE,
          etBOOL,
          { &bMol },
          "Index contains molecule numbers instead of atom numbers" },
        { "-nojump", FALSE, etBOOL, { &bNoJump }, "Remove jumps of atoms across the box" },
        { "-x", FALSE, etBOOL, { &bX }, "Plot X-component" },
        { "-y", FALSE, etBOOL, { &bY }, "Plot Y-component" },
        { "-z", FALSE, etBOOL, { &bZ }, "Plot Z-component" },
        { "-ng", FALSE, etINT, { &ngroups }, "Number of groups to consider" },
        { "-len", FALSE, etBOOL, { &bNorm }, "Plot vector length" },
        { "-fp", FALSE, etBOOL, { &bFP }, "Full precision output" },
        { "-bin", FALSE, etREAL, { &binwidth }, "Binwidth for velocity histogram (nm/ps)" },
        { "-ctime",
          FALSE,
          etREAL,
          { &ctime },
          "Use frame at this time for x in [TT]-cv[tt] and [TT]-cf[tt] instead of the average x" },
        { "-scale",
          FALSE,
          etREAL,
          { &scale },
          "Scale factor for [REF].pdb[ref] output, 0 is autoscale" }
    };
    FILE *       outx = nullptr, *outv = nullptr, *outf = nullptr, *outb = nullptr, *outt = nullptr;
    FILE *       outekt = nullptr, *outekr = nullptr;
    t_topology   top;
    PbcType      pbcType;
    real *       mass, time;
    const char*  indexfn;
    t_trxframe   fr;
    int          flags, nvhisto = 0, *vhisto = nullptr;
    rvec *       xtop, *xp = nullptr;
    rvec *       sumx = nullptr, *sumv = nullptr, *sumf = nullptr;
    matrix       topbox;
    t_trxstatus* status;
    t_trxstatus* status_out = nullptr;
    gmx_rmpbc_t  gpbc       = nullptr;
    int          i, j;
    int          nr_xfr, nr_vfr, nr_ffr;
    char**       grpname;
    int *        isize0, *isize;
    int **       index0, **index;
    int*         atndx;
    t_block*     mols;
    gmx_bool     bTop, bOX, bOXT, bOV, bOF, bOB, bOT, bEKT, bEKR, bCV, bCF;
    gmx_bool     bDim[4], bDum[4], bVD;
    char         sffmt[STRLEN];
    std::array<std::string, 6> box_leg = { "XX", "YY", "ZZ", "YX", "ZX", "ZY" };
    gmx_output_env_t*          oenv;

    t_filenm fnm[] = {
        { efTRX, "-f", nullptr, ffREAD },       { efTPS, nullptr, nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffOPTRD },   { efXVG, "-ox", "coord", ffOPTWR },
        { efTRX, "-oxt", "coord", ffOPTWR },    { efXVG, "-ov", "veloc", ffOPTWR },
        { efXVG, "-of", "force", ffOPTWR },     { efXVG, "-ob", "box", ffOPTWR },
        { efXVG, "-ot", "temp", ffOPTWR },      { efXVG, "-ekt", "ektrans", ffOPTWR },
        { efXVG, "-ekr", "ekrot", ffOPTWR },    { efXVG, "-vd", "veldist", ffOPTWR },
        { efPDB, "-cv", "veloc", ffOPTWR },     { efPDB, "-cf", "force", ffOPTWR },
        { efXVG, "-av", "all_veloc", ffOPTWR }, { efXVG, "-af", "all_force", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc,
                           argv,
                           PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW,
                           NFILE,
                           fnm,
                           asize(pa),
                           pa,
                           asize(desc),
                           desc,
                           0,
                           nullptr,
                           &oenv))
    {
        return 0;
    }

    if (bMol)
    {
        fprintf(stderr,
                "Interpreting indexfile entries as molecules.\n"
                "Using center of mass.\n");
    }

    bOX  = opt2bSet("-ox", NFILE, fnm);
    bOXT = opt2bSet("-oxt", NFILE, fnm);
    bOV  = opt2bSet("-ov", NFILE, fnm);
    bOF  = opt2bSet("-of", NFILE, fnm);
    bOB  = opt2bSet("-ob", NFILE, fnm);
    bOT  = opt2bSet("-ot", NFILE, fnm);
    bEKT = opt2bSet("-ekt", NFILE, fnm);
    bEKR = opt2bSet("-ekr", NFILE, fnm);
    bCV  = opt2bSet("-cv", NFILE, fnm) || opt2bSet("-av", NFILE, fnm);
    bCF  = opt2bSet("-cf", NFILE, fnm) || opt2bSet("-af", NFILE, fnm);
    bVD  = opt2bSet("-vd", NFILE, fnm) || opt2parg_bSet("-bin", asize(pa), pa);
    if (bMol || bOT || bEKT || bEKR)
    {
        bCom = TRUE;
    }

    bDim[XX]  = bX;
    bDim[YY]  = bY;
    bDim[ZZ]  = bZ;
    bDim[DIM] = bNorm;

    if (bFP)
    {
        sprintf(sffmt, "\t%s", gmx_real_fullprecision_pfmt);
    }
    else
    {
        sprintf(sffmt, "\t%%g");
    }
    std::string sffmt6 = gmx::formatString("%s%s%s%s%s%s", sffmt, sffmt, sffmt, sffmt, sffmt, sffmt);

    bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm),
                         &top,
                         &pbcType,
                         &xtop,
                         nullptr,
                         topbox,
                         bCom && (bOX || bOXT || bOV || bOT || bEKT || bEKR));
    sfree(xtop);
    if ((bMol || bCV || bCF) && !bTop)
    {
        gmx_fatal(FARGS, "Need a run input file for option -mol, -cv or -cf");
    }

    if (bMol)
    {
        indexfn = ftp2fn(efNDX, NFILE, fnm);
    }
    else
    {
        indexfn = ftp2fn_null(efNDX, NFILE, fnm);
    }

    if (!(bCom && !bMol))
    {
        ngroups = 1;
    }
    snew(grpname, ngroups);
    snew(isize0, ngroups);
    snew(index0, ngroups);
    get_index(&(top.atoms), indexfn, ngroups, isize0, index0, grpname);

    if (bMol)
    {
        mols    = &(top.mols);
        atndx   = mols->index;
        ngroups = isize0[0];
        snew(isize, ngroups);
        snew(index, ngroups);
        for (i = 0; i < ngroups; i++)
        {
            if (index0[0][i] < 0 || index0[0][i] >= mols->nr)
            {
                gmx_fatal(FARGS, "Molecule index (%d) is out of range (%d-%d)", index0[0][i] + 1, 1, mols->nr);
            }
            isize[i] = atndx[index0[0][i] + 1] - atndx[index0[0][i]];
            snew(index[i], isize[i]);
            for (j = 0; j < isize[i]; j++)
            {
                index[i][j] = atndx[index0[0][i]] + j;
            }
        }
    }
    else
    {
        isize = isize0;
        index = index0;
    }
    if (bCom)
    {
        snew(mass, top.atoms.nr);
        for (i = 0; i < top.atoms.nr; i++)
        {
            mass[i] = top.atoms.atom[i].m;
        }
    }
    else
    {
        mass = nullptr;
    }

    flags = 0;
    std::string label(output_env_get_xvgr_tlabel(oenv));
    if (bOX)
    {
        flags = flags | TRX_READ_X;
        outx  = xvgropen(opt2fn("-ox", NFILE, fnm),
                        bCom ? "Center of mass" : "Coordinate",
                        label,
                        "Coordinate (nm)",
                        oenv);
        make_legend(outx, ngroups, isize0[0], index0[0], grpname, bCom, bMol, bDim, oenv);
    }
    if (bOXT)
    {
        flags      = flags | TRX_READ_X;
        status_out = open_trx(opt2fn("-oxt", NFILE, fnm), "w");
    }
    if (bOV)
    {
        flags = flags | TRX_READ_V;
        outv  = xvgropen(opt2fn("-ov", NFILE, fnm),
                        bCom ? "Center of mass velocity" : "Velocity",
                        label,
                        "Velocity (nm/ps)",
                        oenv);
        make_legend(outv, ngroups, isize0[0], index0[0], grpname, bCom, bMol, bDim, oenv);
    }
    if (bOF)
    {
        flags = flags | TRX_READ_F;
        outf  = xvgropen(
                opt2fn("-of", NFILE, fnm), "Force", label, "Force (kJ mol\\S-1\\N nm\\S-1\\N)", oenv);
        make_legend(outf, ngroups, isize0[0], index0[0], grpname, bCom, bMol, bDim, oenv);
    }
    if (bOB)
    {
        outb = xvgropen(opt2fn("-ob", NFILE, fnm), "Box vector elements", label, "(nm)", oenv);

        xvgrLegend(outb, box_leg, oenv);
    }
    if (bOT)
    {
        bDum[XX]  = FALSE;
        bDum[YY]  = FALSE;
        bDum[ZZ]  = FALSE;
        bDum[DIM] = TRUE;
        flags     = flags | TRX_READ_V;
        outt      = xvgropen(opt2fn("-ot", NFILE, fnm), "Temperature", label, "(K)", oenv);
        make_legend(outt, ngroups, isize[0], index[0], grpname, bCom, bMol, bDum, oenv);
    }
    if (bEKT)
    {
        bDum[XX]  = FALSE;
        bDum[YY]  = FALSE;
        bDum[ZZ]  = FALSE;
        bDum[DIM] = TRUE;
        flags     = flags | TRX_READ_V;
        outekt    = xvgropen(opt2fn("-ekt", NFILE, fnm),
                          "Center of mass translation",
                          label,
                          "Energy (kJ mol\\S-1\\N)",
                          oenv);
        make_legend(outekt, ngroups, isize[0], index[0], grpname, bCom, bMol, bDum, oenv);
    }
    if (bEKR)
    {
        bDum[XX]  = FALSE;
        bDum[YY]  = FALSE;
        bDum[ZZ]  = FALSE;
        bDum[DIM] = TRUE;
        flags     = flags | TRX_READ_X | TRX_READ_V;
        outekr    = xvgropen(opt2fn("-ekr", NFILE, fnm),
                          "Center of mass rotation",
                          label,
                          "Energy (kJ mol\\S-1\\N)",
                          oenv);
        make_legend(outekr, ngroups, isize[0], index[0], grpname, bCom, bMol, bDum, oenv);
    }
    if (bVD)
    {
        flags = flags | TRX_READ_V;
    }
    if (bCV)
    {
        flags = flags | TRX_READ_X | TRX_READ_V;
    }
    if (bCF)
    {
        flags = flags | TRX_READ_X | TRX_READ_F;
    }
    if ((flags == 0) && !bOB)
    {
        fprintf(stderr, "Please select one or more output file options\n");
        exit(0);
    }

    read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);


    if ((bOV || bOF) && fn2ftp(ftp2fn(efTRX, NFILE, fnm)) == efXTC)
    {
        gmx_fatal(FARGS,
                  "Cannot extract velocities or forces since your input XTC file does not contain "
                  "them.");
    }

    if (bCV || bCF)
    {
        snew(sumx, fr.natoms);
    }
    if (bCV)
    {
        snew(sumv, fr.natoms);
    }
    if (bCF)
    {
        snew(sumf, fr.natoms);
    }
    nr_xfr = 0;
    nr_vfr = 0;
    nr_ffr = 0;

    if (bCom && bPBC)
    {
        gpbc = gmx_rmpbc_init(&top.idef, pbcType, fr.natoms);
    }

    do
    {
        time = output_env_conv_time(oenv, fr.time);

        if (fr.bX && bNoJump && fr.bBox)
        {
            if (xp)
            {
                remove_jump(fr.box, fr.natoms, xp, fr.x);
            }
            else
            {
                snew(xp, fr.natoms);
            }
            for (i = 0; i < fr.natoms; i++)
            {
                copy_rvec(fr.x[i], xp[i]);
            }
        }

        if (fr.bX && bCom && bPBC)
        {
            gmx_rmpbc_trxfr(gpbc, &fr);
        }

        if (bVD && fr.bV)
        {
            update_histo(isize[0], index[0], fr.v, &nvhisto, &vhisto, binwidth);
        }

        if (bOX && fr.bX)
        {
            print_data(outx, time, fr.x, mass, bCom, ngroups, isize, index, bDim, sffmt);
        }
        if (bOXT && fr.bX)
        {
            t_trxframe frout = fr;
            if (!frout.bAtoms)
            {
                frout.atoms  = &top.atoms;
                frout.bAtoms = TRUE;
            }
            frout.bV = FALSE;
            frout.bF = FALSE;
            write_trx_x(status_out, &frout, mass, bCom, ngroups, isize, index);
        }
        if (bOV && fr.bV)
        {
            print_data(outv, time, fr.v, mass, bCom, ngroups, isize, index, bDim, sffmt);
        }
        if (bOF && fr.bF)
        {
            print_data(outf, time, fr.f, nullptr, bCom, ngroups, isize, index, bDim, sffmt);
        }
        if (bOB && fr.bBox)
        {
            fprintf(outb, "\t%g", fr.time);
            fprintf(outb,
                    sffmt6.c_str(),
                    fr.box[XX][XX],
                    fr.box[YY][YY],
                    fr.box[ZZ][ZZ],
                    fr.box[YY][XX],
                    fr.box[ZZ][XX],
                    fr.box[ZZ][YY]);
            fprintf(outb, "\n");
        }
        if (bOT && fr.bV)
        {
            fprintf(outt, " %g", time);
            for (i = 0; i < ngroups; i++)
            {
                fprintf(outt, sffmt, temp(fr.v, mass, isize[i], index[i]));
            }
            fprintf(outt, "\n");
        }
        if (bEKT && fr.bV)
        {
            fprintf(outekt, " %g", time);
            for (i = 0; i < ngroups; i++)
            {
                fprintf(outekt, sffmt, ektrans(fr.v, mass, isize[i], index[i]));
            }
            fprintf(outekt, "\n");
        }
        if (bEKR && fr.bX && fr.bV)
        {
            fprintf(outekr, " %g", time);
            for (i = 0; i < ngroups; i++)
            {
                fprintf(outekr, sffmt, ekrot(fr.x, fr.v, mass, isize[i], index[i]));
            }
            fprintf(outekr, "\n");
        }
        if ((bCV || bCF) && fr.bX
            && (ctime < 0 || (fr.time >= ctime * 0.999999 && fr.time <= ctime * 1.000001)))
        {
            for (i = 0; i < fr.natoms; i++)
            {
                rvec_inc(sumx[i], fr.x[i]);
            }
            nr_xfr++;
        }
        if (bCV && fr.bV)
        {
            for (i = 0; i < fr.natoms; i++)
            {
                rvec_inc(sumv[i], fr.v[i]);
            }
            nr_vfr++;
        }
        if (bCF && fr.bF)
        {
            for (i = 0; i < fr.natoms; i++)
            {
                rvec_inc(sumf[i], fr.f[i]);
            }
            nr_ffr++;
        }

    } while (read_next_frame(oenv, status, &fr));

    if (gpbc != nullptr)
    {
        gmx_rmpbc_done(gpbc);
    }

    /* clean up a bit */
    close_trx(status);

    if (bOX)
    {
        xvgrclose(outx);
    }
    if (bOXT)
    {
        close_trx(status_out);
    }
    if (bOV)
    {
        xvgrclose(outv);
    }
    if (bOF)
    {
        xvgrclose(outf);
    }
    if (bOB)
    {
        xvgrclose(outb);
    }
    if (bOT)
    {
        xvgrclose(outt);
    }
    if (bEKT)
    {
        xvgrclose(outekt);
    }
    if (bEKR)
    {
        xvgrclose(outekr);
    }

    if (bVD)
    {
        print_histo(opt2fn("-vd", NFILE, fnm), nvhisto, vhisto, binwidth, oenv);
    }

    if (bCV || bCF)
    {
        if (nr_xfr > 1)
        {
            if (pbcType != PbcType::No && !bNoJump)
            {
                fprintf(stderr,
                        "\nWARNING: More than one frame was used for option -cv or -cf\n"
                        "If atoms jump across the box you should use the -nojump or -ctime "
                        "option\n\n");
            }
            for (i = 0; i < isize[0]; i++)
            {
                svmul(1.0 / nr_xfr, sumx[index[0][i]], sumx[index[0][i]]);
            }
        }
        else if (nr_xfr == 0)
        {
            fprintf(stderr, "\nWARNING: No coordinate frames found for option -cv or -cf\n\n");
        }
    }
    if (bCV)
    {
        write_pdb_bfac(opt2fn("-cv", NFILE, fnm),
                       opt2fn("-av", NFILE, fnm),
                       "average velocity",
                       &(top.atoms),
                       pbcType,
                       topbox,
                       isize[0],
                       index[0],
                       nr_xfr,
                       sumx,
                       nr_vfr,
                       sumv,
                       bDim,
                       scale,
                       oenv);
    }
    if (bCF)
    {
        write_pdb_bfac(opt2fn("-cf", NFILE, fnm),
                       opt2fn("-af", NFILE, fnm),
                       "average force",
                       &(top.atoms),
                       pbcType,
                       topbox,
                       isize[0],
                       index[0],
                       nr_xfr,
                       sumx,
                       nr_ffr,
                       sumf,
                       bDim,
                       scale,
                       oenv);
    }

    /* view it */
    view_all(oenv, NFILE, fnm);

    done_top(&top);
    // Free index and isize only if they are distinct from index0 and isize0
    if (bMol)
    {
        for (int i = 0; i < ngroups; i++)
        {
            sfree(index[i]);
        }
        sfree(index);
        sfree(isize);
    }
    for (int i = 0; i < ngroups; i++)
    {
        sfree(index0[i]);
        sfree(grpname[i]);
    }
    sfree(index0);
    sfree(isize0);
    sfree(grpname);
    done_frame(&fr);
    output_env_done(oenv);

    return 0;
}
