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

#include "filetypes.h"

#include <cstring>

#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"

enum
{
    eftASC, eftXDR, eftTNG, eftGEN, eftNR
};

/* To support multiple file types with one general (eg TRX) we have
 * these arrays.
 */
static const int trxs[] =
{
    efXTC, efTRR, efCPT,
    efGRO, efG96, efPDB, efTNG
};
#define NTRXS asize(trxs)

static const int trcompressed[] =
{
    efXTC,
    efTNG
};
#define NTRCOMPRESSED asize(trcompressed)

static const int tros[] =
{
    efXTC, efTRR,
    efGRO, efG96, efPDB, efTNG
};
#define NTROS asize(tros)

static const int trns[] =
{
    efTRR, efCPT,
    efTNG
};
#define NTRNS asize(trns)

static const int stos[] =
{ efGRO, efG96, efPDB, efBRK, efENT, efESP };
#define NSTOS asize(stos)

static const int stxs[] =
{
    efGRO, efG96, efPDB, efBRK, efENT, efESP,
    efTPR
};
#define NSTXS asize(stxs)

static const int tpss[] =
{
    efTPR,
    efGRO, efG96, efPDB, efBRK, efENT
};
#define NTPSS asize(tpss)

#ifdef BUILD_WITH_FDA
// FDA scalar stress type
static const int tstrs[] =
{
    efPSA, efPSR, efVMA
};
#define NTSTRS asize(tstrs)

// FDA atomic or residue based force pairs
static const int tpfx[] =
{
    efPFA, efPFR
};
#define NTPFX asize(tpfx)

// FDA resulting networks (pdb or dimacs)
static const int tgrx[] =
{
    efPDB, efDIM
};
#define NTGRX asize(tgrx)

// FDA atomic or residue based punctual stress
static const int tpsx[] =
{
    efPSA, efPSR
};
#define NTPSX asize(tpsx)

// FDA view stress format (pdb or xpm)
static const int tvst[] =
{
    efPDB, efXPM
};
#define NTVST asize(tvst)
#endif

typedef struct
{
    int         ftype;
    const char *ext;
    const char *defnm;
    const char *defopt;
    const char *descr;
    int         ntps;
    const int  *tps;
} t_deffile;

/* this array should correspond to the enum in filetypes.h */
static const t_deffile deffile[efNR] =
{
    { eftASC, ".mdp", "grompp", "-f", "grompp input file with MD parameters" },
    { eftGEN, ".???", "traj", "-f", "Trajectory", NTRXS, trxs },
    { eftGEN, ".???", "trajout", "-f", "Trajectory", NTROS, tros },
    { eftGEN, ".???", "traj", nullptr,
      "Full precision trajectory", NTRNS, trns },
    { eftXDR, ".trr", "traj", nullptr, "Trajectory in portable xdr format" },
    { eftGEN, ".???", "traj_comp", nullptr,
      "Compressed trajectory (tng format or portable xdr format)", NTRCOMPRESSED, trcompressed},
    { eftXDR, ".xtc", "traj", nullptr,
      "Compressed trajectory (portable xdr format): xtc" },
    { eftTNG, ".tng", "traj", nullptr,
      "Trajectory file (tng format)" },
    { eftXDR, ".edr", "ener",   nullptr, "Energy file"},
    { eftGEN, ".???", "conf", "-c", "Structure file", NSTXS, stxs },
    { eftGEN, ".???", "out", "-o", "Structure file", NSTOS, stos },
    { eftASC, ".gro", "conf", "-c", "Coordinate file in Gromos-87 format" },
    { eftASC, ".g96", "conf", "-c", "Coordinate file in Gromos-96 format" },
    { eftASC, ".pdb", "eiwit",  "-f", "Protein data bank file"},
    { eftASC, ".brk", "eiwit",  "-f", "Brookhaven data bank file"},
    { eftASC, ".ent", "eiwit", "-f", "Entry in the protein date bank" },
    { eftASC, ".esp", "conf", "-f", "Coordinate file in Espresso format" },
    { eftASC, ".pqr", "state",  "-o", "Coordinate file for MEAD"},
    { eftXDR, ".cpt", "state",  "-cp", "Checkpoint file"},
    { eftASC, ".log", "run",    "-l", "Log file"},
    { eftASC, ".xvg", "graph",  "-o", "xvgr/xmgr file"},
    { eftASC, ".out", "hello",  "-o", "Generic output file"},
    { eftASC, ".ndx", "index",  "-n", "Index file", },
    { eftASC, ".top", "topol",  "-p", "Topology file"},
    { eftASC, ".itp", "topinc", nullptr, "Include file for topology"},
    { eftGEN, ".???", "topol", "-s", "Structure+mass(db)", NTPSS, tpss },
    { eftXDR, ".tpr", "topol",  "-s", "Portable xdr run input file"},
    { eftASC, ".tex", "doc",    "-o", "LaTeX file"},
    { eftASC, ".rtp", "residue", nullptr, "Residue Type file used by pdb2gmx" },
    { eftASC, ".atp", "atomtp", nullptr, "Atomtype file used by pdb2gmx" },
    { eftASC, ".hdb", "polar",  nullptr, "Hydrogen data base"},
    { eftASC, ".dat", "nnnice", nullptr, "Generic data file"},
    { eftASC, ".dlg", "user",   nullptr, "Dialog Box data for ngmx"},
    { eftASC, ".map", "ss", nullptr, "File that maps matrix data to colors" },
    { eftASC, ".eps", "plot", nullptr, "Encapsulated PostScript (tm) file" },
    { eftASC, ".mat", "ss",     nullptr, "Matrix Data file"},
    { eftASC, ".m2p", "ps",     nullptr, "Input file for mat2ps"},
    { eftXDR, ".mtx", "hessian", "-m", "Hessian matrix"},
    { eftASC, ".edi", "sam",    nullptr, "ED sampling input"},
    { eftASC, ".cub", "pot",  nullptr, "Gaussian cube file" },
    { eftASC, ".xpm", "root", nullptr, "X PixMap compatible matrix file" },
#ifdef BUILD_WITH_FDA
    { eftASC, ".pfi", "fda", NULL, "FDA paired force input file" },
    { eftASC, ".pfa", "fda", NULL, "FDA paired force output atom-based matrix file" },
    { eftASC, ".pfr", "fda", NULL, "FDA paired force output residue-based matrix file" },
    { eftASC, ".psa", "fda", NULL, "FDA punctual stress output atom-based matrix file" },
    { eftASC, ".psr", "fda", NULL, "FDA punctual stress output residue-based matrix file" },
    { eftASC, ".vsa", "fda", NULL, "FDA virial stress output atom-based matrix file" },
    { eftASC, ".vma", "fda", NULL, "FDA von Mises virial stress output atom-based matrix file" },
    { eftGEN, ".???", "fda", "-f", "FDA punctual or virial stress", NTSTRS, tstrs },
    { eftGEN, ".???", "fda", "-f", "FDA atomic or residue based force pairs", NTPFX, tpfx },
    { eftGEN, ".???", "fda", "-o", "FDA resulting networks (pdb or dimacs)", NTGRX, tgrx },
    { eftASC, ".dmc", "fda", NULL, "FDA network graph DIMACS" },
    { eftGEN, ".???", "fda", NULL, "FDA atomic or residue based punctual stress", NTPSX, tpsx },
    { eftGEN, ".???", "fda", NULL, "FDA view stress format (pdb or xpm)", NTVST, tvst },
    { eftASC, ".pml", "fda", NULL, "FDA pymol script for pdb-trajectory" },
#endif
    { eftASC, "", "rundir", nullptr, "Run directory" }
};

const char *ftp2ext(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ext[0] != '\0' ? deffile[ftp].ext + 1 : "";
    }
    else
    {
        return "unknown";
    }
}

const char *ftp2ext_generic(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        switch (ftp)
        {
            case efTRX:
                return "trx";
            case efTRN:
                return "trn";
            case efSTO:
                return "sto";
            case efSTX:
                return "stx";
            case efTPS:
                return "tps";
            default:
                return ftp2ext(ftp);
        }
    }
    else
    {
        return "unknown";
    }
}

const char *ftp2ext_with_dot(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ext;
    }
    else
    {
        return "unknown";
    }
}

int ftp2generic_count(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].ntps;
    }
    else
    {
        return 0;
    }
}

const int *ftp2generic_list(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].tps;
    }
    else
    {
        return nullptr;
    }
}

const char *ftp2desc(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].descr;
    }
    else
    {
        return "unknown filetype";
    }
}

gmx_bool ftp_is_text(int ftp)
{
    if ((ftp >= 0) && (ftp < efNR))
    {
        return deffile[ftp].ftype == eftASC;
    }
    return FALSE;
}

gmx_bool ftp_is_xdr(int ftp)
{
    if ((ftp >= 0) && (ftp < efNR))
    {
        return deffile[ftp].ftype == eftXDR;
    }
    return FALSE;
}

const char *ftp2defnm(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].defnm;
    }
    else
    {
        return nullptr;
    }
}

const char *ftp2defopt(int ftp)
{
    if ((0 <= ftp) && (ftp < efNR))
    {
        return deffile[ftp].defopt;
    }
    else
    {
        return nullptr;
    }
}

int fn2ftp(const char *fn)
{
    int         i, len;
    const char *feptr;
    const char *eptr;

    if (!fn)
    {
        return efNR;
    }

    len = std::strlen(fn);
    if ((len >= 4) && (fn[len - 4] == '.'))
    {
        feptr = &(fn[len - 4]);
    }
    else
    {
        return efNR;
    }

    for (i = 0; (i < efNR); i++)
    {
        if ((eptr = deffile[i].ext) != nullptr)
        {
            if (gmx_strcasecmp(feptr, eptr) == 0)
            {
                break;
            }
        }
    }

    return i;
}
