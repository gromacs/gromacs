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
#include <cstring>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/fitahx.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/hxprops.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;
struct gmx_output_env_t;

int gmx_helix(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes all kinds of helix properties. First, the peptide",
        "is checked to find the longest helical part, as determined by",
        "hydrogen bonds and [GRK]phi[grk]/[GRK]psi[grk] angles.",
        "That bit is fitted",
        "to an ideal helix around the [IT]z[it]-axis and centered around the origin.",
        "Then the following properties are computed:",
        "",
        " * Helix radius (file [TT]radius.xvg[tt]). This is merely the",
        "   RMS deviation in two dimensions for all C[GRK]alpha[grk] atoms.",
        "   it is calculated as [SQRT]([SUM][sum][SUB]i[sub] (x^2(i)+y^2(i)))/N[sqrt] where N is ",
        "   the number of backbone atoms. For an ideal helix the radius is 0.23 nm.",
        " * Twist (file [TT]twist.xvg[tt]). The average helical angle per",
        "   residue is calculated. For an [GRK]alpha[grk]-helix it is 100 degrees,",
        "   for 3-10 helices it will be smaller, and ",
        "   for 5-helices it will be larger.",
        " * Rise per residue (file [TT]rise.xvg[tt]). The helical rise per",
        "   residue is plotted as the difference in [IT]z[it]-coordinate between C[GRK]alpha[grk]",
        "   atoms. For an ideal helix, this is 0.15 nm.",
        " * Total helix length (file [TT]len-ahx.xvg[tt]). The total length",
        "   of the",
        "   helix in nm. This is simply the average rise (see above) times the",
        "   number of helical residues (see below).",
        " * Helix dipole, backbone only (file [TT]dip-ahx.xvg[tt]).",
        " * RMS deviation from ideal helix, calculated for the C[GRK]alpha[grk]",
        "   atoms only (file [TT]rms-ahx.xvg[tt]).",
        " * Average C[GRK]alpha[grk] - C[GRK]alpha[grk] dihedral angle (file [TT]phi-ahx.xvg[tt]).",
        " * Average [GRK]phi[grk] and [GRK]psi[grk] angles (file [TT]phipsi.xvg[tt]).",
        " * Ellipticity at 222 nm according to Hirst and Brooks."
    };
    static gmx_bool bCheck = FALSE, bFit = TRUE, bDBG = FALSE, bEV = FALSE;
    static int      rStart = 0, rEnd = 0, r0 = 1;
    t_pargs pa[] = { { "-r0", FALSE, etINT, { &r0 }, "The first residue number in the sequence" },
                     { "-q",
                       FALSE,
                       etBOOL,
                       { &bCheck },
                       "Check at every step which part of the sequence is helical" },
                     { "-F", FALSE, etBOOL, { &bFit }, "Toggle fit to a perfect helix" },
                     { "-db", FALSE, etBOOL, { &bDBG }, "Print debug info" },
                     { "-ev", FALSE, etBOOL, { &bEV }, "Write a new 'trajectory' file for ED" },
                     { "-ahxstart", FALSE, etINT, { &rStart }, "First residue in helix" },
                     { "-ahxend", FALSE, etINT, { &rEnd }, "Last residue in helix" } };

    typedef struct
    { //NOLINT(clang-analyzer-optin.performance.Padding)
        FILE *      fp, *fp2;
        gmx_bool    bfp2;
        const char* filenm;
        const char* title;
        const char* xaxis;
        const char* yaxis;
        real        val;
    } t_xvgrfile;

    t_xvgrfile xf[efhNR] = {
        { nullptr, nullptr, TRUE, "radius", "Helix radius", nullptr, "r (nm)", 0.0 },
        { nullptr, nullptr, TRUE, "twist", "Twist per residue", nullptr, "Angle (deg)", 0.0 },
        { nullptr, nullptr, TRUE, "rise", "Rise per residue", nullptr, "Rise (nm)", 0.0 },
        { nullptr, nullptr, FALSE, "len-ahx", "Length of the Helix", nullptr, "Length (nm)", 0.0 },
        { nullptr, nullptr, FALSE, "dip-ahx", "Helix Backbone Dipole", nullptr, "rq (nm e)", 0.0 },
        { nullptr, nullptr, TRUE, "rms-ahx", "RMS Deviation from Ideal Helix", nullptr, "RMS (nm)", 0.0 },
        { nullptr, nullptr, FALSE, "rmsa-ahx", "Average RMSD per Residue", "Residue", "RMS (nm)", 0.0 },
        { nullptr, nullptr, FALSE, "cd222", "Ellipticity at 222 nm", nullptr, "nm", 0.0 },
        { nullptr, nullptr, TRUE, "pprms", "RMS Distance from \\8a\\4-helix", nullptr, "deg", 0.0 },
        { nullptr, nullptr, TRUE, "caphi", "Average Ca-Ca Dihedral", nullptr, "\\8F\\4(deg)", 0.0 },
        { nullptr, nullptr, TRUE, "phi", "Average \\8F\\4 angles", nullptr, "deg", 0.0 },
        { nullptr, nullptr, TRUE, "psi", "Average \\8Y\\4 angles", nullptr, "deg", 0.0 },
        { nullptr, nullptr, TRUE, "hb3", "Average n-n+3 hbond length", nullptr, "nm", 0.0 },
        { nullptr, nullptr, TRUE, "hb4", "Average n-n+4 hbond length", nullptr, "nm", 0.0 },
        { nullptr, nullptr, TRUE, "hb5", "Average n-n+5 hbond length", nullptr, "nm", 0.0 },
        { nullptr, nullptr, FALSE, "JCaHa", "J-Coupling Values", "Residue", "Hz", 0.0 },
        { nullptr, nullptr, FALSE, "helicity", "Helicity per Residue", "Residue", "% of time", 0.0 }
    };

    gmx_output_env_t* oenv;
    char              buf[54];
    t_trxstatus*      status;
    int               natoms, nres;
    t_bb*             bb;
    int               i, j, nall, nbb, nca, teller;
    int *             bbindex, *caindex, *allindex;
    t_topology*       top;
    PbcType           pbcType;
    rvec *            x, *xref;
    real              t;
    real              rms;
    matrix            box;
    gmx_rmpbc_t       gpbc = nullptr;
    gmx_bool          bRange;
    t_filenm          fnm[] = {
        { efTPR, nullptr, nullptr, ffREAD },
        { efNDX, nullptr, nullptr, ffREAD },
        { efTRX, "-f", nullptr, ffREAD },
        { efSTO, "-cz", "zconf", ffWRITE },
    };
#define NFILE asize(fnm)

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    bRange = (opt2parg_bSet("-ahxstart", asize(pa), pa) && opt2parg_bSet("-ahxend", asize(pa), pa));

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &pbcType);

    natoms = read_first_x(oenv, &status, opt2fn("-f", NFILE, fnm), &t, &x, box);

    if (natoms != top->atoms.nr)
    {
        gmx_fatal(FARGS,
                  "Sorry can only run when the number of atoms in the run input file (%d) is equal "
                  "to the number in the trajectory (%d)",
                  top->atoms.nr,
                  natoms);
    }

    bb = mkbbind(ftp2fn(efNDX, NFILE, fnm),
                 &nres,
                 &nbb,
                 r0,
                 &nall,
                 &allindex,
                 top->atoms.atomname,
                 top->atoms.atom,
                 top->atoms.resinfo);
    snew(bbindex, natoms);
    snew(caindex, nres);

    fprintf(stderr, "nall=%d\n", nall);

    /* Open output files, default x-axis is time */
    for (i = 0; (i < efhNR); i++)
    {
        sprintf(buf, "%s.xvg", xf[i].filenm);
        remove(buf);
        xf[i].fp = xvgropen(buf, xf[i].title, xf[i].xaxis ? xf[i].xaxis : "Time (ps)", xf[i].yaxis, oenv);
        if (xf[i].bfp2)
        {
            sprintf(buf, "%s.out", xf[i].filenm);
            remove(buf);
            xf[i].fp2 = gmx_ffopen(buf, "w");
        }
    }

    /* Read reference frame from tpx file to compute helix length */
    snew(xref, top->atoms.nr);
    read_tpx(ftp2fn(efTPR, NFILE, fnm), nullptr, nullptr, nullptr, xref, nullptr, nullptr);
    calc_hxprops(nres, bb, xref);
    do_start_end(nres, bb, &nbb, bbindex, &nca, caindex, bRange, rStart, rEnd);
    sfree(xref);
    if (bDBG)
    {
        fprintf(stderr, "nca=%d, nbb=%d\n", nca, nbb);
        pr_bb(stdout, nres, bb);
    }

    gpbc = gmx_rmpbc_init(&top->idef, pbcType, natoms);

    teller = 0;
    do
    {
        if ((teller++ % 10) == 0)
        {
            fprintf(stderr, "\rt=%.2f", t);
            fflush(stderr);
        }
        gmx_rmpbc_apply(gpbc, natoms, box, x);


        calc_hxprops(nres, bb, x);
        if (bCheck)
        {
            do_start_end(nres, bb, &nbb, bbindex, &nca, caindex, FALSE, 0, 0);
        }

        if (nca >= 5)
        {
            rms = fit_ahx(nres, bb, natoms, nall, allindex, x, nca, caindex, bFit);

            if (teller == 1)
            {
                write_sto_conf(
                        opt2fn("-cz", NFILE, fnm), "Helix fitted to Z-Axis", &(top->atoms), x, nullptr, pbcType, box);
            }

            xf[efhRAD].val   = radius(xf[efhRAD].fp2, nca, caindex, x);
            xf[efhTWIST].val = twist(nca, caindex, x);
            xf[efhRISE].val  = rise(nca, caindex, x);
            xf[efhLEN].val   = ahx_len(nca, caindex, x);
            xf[efhCD222].val = ellipticity(nres, bb);
            xf[efhDIP].val   = dip(nbb, bbindex, x, top->atoms.atom);
            xf[efhRMS].val   = rms;
            xf[efhCPHI].val  = ca_phi(nca, caindex, x);
            xf[efhPPRMS].val = pprms(xf[efhPPRMS].fp2, nres, bb);

            for (j = 0; (j <= efhCPHI); j++)
            {
                fprintf(xf[j].fp, "%10g  %10g\n", t, xf[j].val);
            }

            av_phipsi(xf[efhPHI].fp, xf[efhPSI].fp, xf[efhPHI].fp2, xf[efhPSI].fp2, t, nres, bb);
            av_hblen(xf[efhHB3].fp,
                     xf[efhHB3].fp2,
                     xf[efhHB4].fp,
                     xf[efhHB4].fp2,
                     xf[efhHB5].fp,
                     xf[efhHB5].fp2,
                     t,
                     nres,
                     bb);
        }
    } while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");

    gmx_rmpbc_done(gpbc);

    close_trx(status);

    for (i = 0; (i < nres); i++)
    {
        if (bb[i].nrms > 0)
        {
            fprintf(xf[efhRMSA].fp, "%10d  %10g\n", r0 + i, bb[i].rmsa / bb[i].nrms);
        }
        fprintf(xf[efhAHX].fp, "%10d  %10g\n", r0 + i, (bb[i].nhx * 100.0) / static_cast<real>(teller));
        fprintf(xf[efhJCA].fp, "%10d  %10g\n", r0 + i, 140.3 + (bb[i].jcaha / static_cast<double>(teller)));
    }

    for (i = 0; (i < efhNR); i++)
    {
        xvgrclose(xf[i].fp);
        if (xf[i].bfp2)
        {
            xvgrclose(xf[i].fp2);
        }
        do_view(oenv, xf[i].filenm, "-nxy");
    }

    return 0;
}
