/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <ctype.h>
#include <stdlib.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "babel_io.h"
#include "fill_inputrec.h"
#include "getmdlogger.h"
#include "molprop_xml.h"
#include "mymol.h"
#include "poldata_xml.h"

int alex_gentop(int argc, char *argv[])
{
    static const char               *desc[] = {
        "gentop generates a topology from molecular coordinates",
        "either from a file, from a database, or from a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.[PAR]",
        "If the [TT]-oi[tt] option is set an [TT]itp[tt] file will be generated",
        "instead of a [TT]top[tt] file.",
        "When [TT]-param[tt] is set, equilibrium distances and angles",
        "and force constants will be printed in the topology for all",
        "interactions. The equilibrium distances and angles are taken",
        "from the input coordinates, the force constant are set with",
        "command line options.",
        "With the [TT]-db molecule[tt] option a file is extracted from the",
        "database from one of the specified QM calculations (given with [TT]-lot[tt]).",
        "An alternative to the system-wide database [TT]molprops.dat[tt]",
        "can be passed along using the [TT]-mpdb[tt] flag.[PAR]",
        "If the flag [TT]-qgen[tt] is given, charges will be generated using the",
        "specified algorithm. Without the flag the charges from the QM calculation",
        "will be used.",
        "The only supported force field for this tool is Alexandria.[PAR]",
        "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
        "The corresponding data files can be found in the library directory",
        "with names like ffXXXX.YYY. Check chapter 5 of the manual for more",
        "information about file formats. The default forcefield is Alexandria",
        "but selection can be made interactive, using the [TT]-ff select[tt] option.",
        "one of the short names above on the command line instead.[PAR]"
    };
    const char                      *bugs[] = {
        "No force constants for impropers are generated"
    };
    gmx_output_env_t                *oenv;
    gmx_atomprop_t                   aps;
    char                             forcefield[STRLEN];
    char                             ffdir[STRLEN];
    char                             ffname[STRLEN];
    std::vector<alexandria::MolProp> mps;
    alexandria::MolPropIterator      mpi;
    alexandria::MyMol                mymol;
    immStatus                        imm;

    t_filenm                         fnm[] = {
        { efSTX, "-f",        "conf",      ffOPTRD },
        { efTOP, "-p",        "out",       ffOPTWR },
        { efITP, "-oi",       "out",       ffOPTWR },
        { efSTO, "-c",        "out",       ffWRITE },
        { efLOG, "-g03",      "gauss",     ffOPTRDMULT},
        { efNDX, "-n",        "renum",     ffOPTWR },
        { efDAT, "-q",        "qout",      ffOPTWR },
        { efDAT, "-mpdb",     "molprops",  ffOPTRD },
        { efDAT, "-d",        "gentop",    ffOPTRD },
        { efXVG, "-table",    "table",     ffOPTRD },
        { efCUB, "-pot",      "potential", ffOPTWR },
        { efCUB, "-ref",      "refpot",    ffOPTRD },
        { efCUB, "-diff",     "diffpot",   ffOPTWR },
        { efCUB, "-rho",      "density",   ffOPTWR },
        { efXVG, "-diffhist", "diffpot",   ffOPTWR },
        { efXVG, "-his",      "pot-histo", ffOPTWR },
        { efXVG, "-pc",       "pot-comp",  ffOPTWR },
        { efPDB, "-pdbdiff",  "pdbdiff",   ffOPTWR },
        { efXVG, "-plotESP",  "ESPcorr",   ffOPTWR }
    };

    const  int                       NFILE          = asize(fnm);
    
    static int                       maxpot         = 100;
    static int                       nsymm          = 0;
    static int                       qcycle         = 1000;
    static int                       nexcl          = 2;
    static real                      qtol           = 1e-6;
    static real                      qtot           = 0;
    static real                      hfac           = 0;
    static real                      watoms         = 0;
    static real                      spacing        = 0.1;
    static real                      efield         = 0;
    static char                     *molnm          = (char *)"";
    static char                     *iupac          = (char *)"";
    static char                     *dbname         = (char *)"";
    static char                     *symm_string    = (char *)"";
    static char                     *conf           = (char *)"minimum";
    static char                     *basis          = (char *)"";
    static char                     *jobtype        = (char *)"unknown";
    static gmx_bool                  bQsym          = false;
    static gmx_bool                  bITP           = false;
    static gmx_bool                  bPairs         = false;
    static gmx_bool                  bUsePDBcharge  = false;
    static gmx_bool                  bGenVSites     = false;
    static gmx_bool                  bDihedral      = false;
    static gmx_bool                  bLOG           = false;
    static gmx_bool                  bCUBE          = false;
    static gmx_bool                  bH14           = true;
    static gmx_bool                  bVerbose       = true;

    static const char               *cqdist[]       = {nullptr, "AXp", "AXg", "AXs", "AXpp", "AXpg", "AXps", "Yang", "Bultinck", "Rappe", nullptr};
    static const char               *cqgen[]        = {nullptr, "None", "ACM", "ESP", "RESP", nullptr};
    static const char               *cgopt[]        = {nullptr, "Atom", "Group", "Neutral", nullptr};
    static const char               *lot            = "AFF/ACM";
    static const char               *ff             = "alexandria";

    t_pargs                          pa[]     = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose output in the top file and on terminal." },
        { "-cube",   FALSE, etBOOL, {&bCUBE},
          "Generate cube." },
        { "-ff",     FALSE, etSTR,  {&ff},
          "Force field, interactive by default. Use -h for information." },
        { "-db",     FALSE, etSTR,  {&dbname},
          "Read a molecule from the database rather than from a file" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges" },
        { "-dih",    FALSE, etBOOL, {&bDihedral},
          "Add dihedrals to the topology" },
        { "-H14",    FALSE, etBOOL, {&bH14},
          "HIDDENUse 3rd neighbour interactions for hydrogen atoms" },
        { "-pairs",  FALSE, etBOOL, {&bPairs},
          "HIDDENOutput 1-4 interactions (pairs) in topology file. Check consistency of your option with the [TT]-nexcl[tt] flag." },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-iupac",   FALSE, etSTR,  {&iupac},
          "IUPAC Name of your molecule" },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-basis",  FALSE, etSTR, {&basis},
          "Basis-set used in this calculation for those case where it is difficult to extract from a Gaussian file" },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Fraction of potential points to read from the gaussian file (percent). If 100 all points are registered, else a selection of points evenly spread over the range of values is taken" },
        { "-nsymm", FALSE, etINT, {&nsymm},
          "Symmetry number of the molecule can be supplied here if you know there is an error in the input file" },
        { "-genvsites", FALSE, etBOOL, {&bGenVSites},
          "Generate virtual sites. Check and double check." },
        { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
          "HIDDENUse the B-factor supplied in a pdb file for the atomic charges" },
        { "-spacing", FALSE, etREAL, {&spacing},
          "Spacing of grid points for computing the potential (not used when a reference file is read)." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use 0. For point+smeared charges 1 is recommended." },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Charge distribution used" },
        { "-qtol",   FALSE, etREAL, {&qtol},
          "Tolerance for assigning charge generation algorithm" },
        { "-qtot",   FALSE, etREAL, {&qtot},
          "Total charge of the molecule. If the input file is a Gaussian log file, qtot will be taken from the log file." },
        { "-qcycle", FALSE, etINT, {&qcycle},
          "Max number of tries for optimizing the charges. The trial with lowest chi2 will be used for generating a topology. Will be turned off if randzeta is No." },
        { "-hfac",    FALSE, etREAL, {&hfac},
          "HIDDENFudge factor for AXx algorithms that modulates J00 for hydrogen atoms by multiplying it by (1 + hfac*qH). This hack is originally due to Rappe & Goddard." },
        { "-qsymm",  FALSE, etBOOL, {&bQsym},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2." },                             
        { "-symm",   FALSE, etSTR, {&symm_string},
          "Use the order given here for symmetrizing, e.g. when specifying [TT]-symm '0 1 0'[tt] for a water molecule (H-O-H) the hydrogens will have obtain the same charge. For simple groups, like methyl (or water) this is done automatically, but higher symmetry is not detected by the program. The numbers should correspond to atom numbers minus 1, and point to either the atom itself or to a previous atom." },
        { "-cgsort", FALSE, etSTR, {cgopt},
          "HIDDENOption for assembling charge groups: based on Atom (default, does not change the atom order), Group (e.g. CH3 groups are kept together), or Neutral sections (try to find groups that together are neutral). If the order of atoms is changed an index file is written in order to facilitate changing the order in old files." },
        { "-nexcl",    FALSE, etINT, {&nexcl},
          "HIDDENNumber of exclusion" },
        { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." },
        { "-efield",  FALSE, etREAL, {&efield},
          "The magnitude of the external electeric field to calculate polarizability tensor." },
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }

    alexandria::Poldata       pd;
    t_inputrec               *inputrec                   = new t_inputrec();
    t_commrec                *cr                         = init_commrec();
    const char               *tabfn                      = opt2fn_null("-table", NFILE, fnm);
    eChargeGroup              ecg                        = (eChargeGroup) get_option(cgopt);
    ChargeGenerationAlgorithm iChargeGenerationAlgorithm = (ChargeGenerationAlgorithm) get_option(cqgen);
    ChargeDistributionModel   iChargeDistributionModel;
    gmx::MDLogger             mdlog                      = getMdLogger(cr, stdout);

    /* Force field selection, interactive or direct */
    choose_ff(strcmp(ff, "select") == 0 ? nullptr : ff,
              forcefield, sizeof(forcefield), ffdir, sizeof(ffdir));

    if (strlen(forcefield) > 0)
    {
        strcpy(ffname, forcefield);
        ffname[0] = toupper(ffname[0]);
    }
    else
    {
        gmx_fatal(FARGS, "Empty forcefield string");
    }

    /* Check the options */
    bITP = opt2bSet("-oi", NFILE, fnm);
    if ((qtol < 0) || (qtol > 1))
    {
        gmx_fatal(FARGS, "Charge tolerance should be between 0 and 1 (not %g)", qtol);
    }
    if ((iChargeDistributionModel = name2eemtype(cqdist[0])) == eqdNR)
    {
        gmx_fatal(FARGS, "Invalid Charge Distribution model %s.\n", cqdist[0]);
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();
    try
    {
        const char *gentop_fnm = opt2fn_null("-d", NFILE, fnm);
        alexandria::readPoldata(gentop_fnm ? gentop_fnm : "", pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    if (pd.getNexcl() != nexcl)
    {
        printf("Exclusion number changed from %d to %d.\n", pd.getNexcl(), nexcl);
        pd.setNexcl(nexcl);
    }
    if (bVerbose)
    {
        printf("Reading force field information. There are %d atomtypes.\n",
               static_cast<int>(pd.getNatypes()));
    }
    if (strlen(dbname) > 0)
    {
        if (bVerbose)
        {
            printf("Reading molecule database.\n");
        }
        MolPropRead(opt2fn_null("-mpdb", NFILE, fnm), mps);
        for (mpi = mps.begin(); (mpi < mps.end()); mpi++)
        {
            if (strcasecmp(dbname, mpi->getMolname().c_str()) == 0)
            {
                break;
            }
        }
        if (mpi == mps.end())
        {
            gmx_fatal(FARGS, "Molecule %s not found in database", dbname);
        }
    }
    else
    {
        gmx::ArrayRef<const std::string> fns;
        if (strlen(molnm) == 0)
        {
            molnm = (char *)"MOL";
        }
        bLOG = opt2bSet("-g03", NFILE, fnm);
        if (bLOG)
        {
            fns = ftp2fns(efLOG, NFILE, fnm);
        }
        else if (opt2bSet("-f", NFILE, fnm))
        {
            fns = ftp2fns(efSTX, NFILE, fnm);
        }
        if (fns.size() > 0)
        {
            for (auto &i : fns)
            {
                alexandria::MolProp  mp;
                readBabel(i.c_str(),
                          mp,
                          molnm,
                          iupac,
                          conf,
                          basis,
                          maxpot,
                          nsymm,
                          pd.getForceField().c_str(),
                          jobtype,
                          qtot);
                mps.push_back(mp);
            }
        }
        else
        {
            gmx_fatal(FARGS, "No input file has been specified.");
        }
    }
    for (auto mpi = mps.begin(); mpi < mps.end(); mpi++)
    {
        mymol.molProp()->Merge(mpi);
    }
    mymol.SetForceField(forcefield);
    fill_inputrec(inputrec);
    mymol.setInputrec(inputrec);
    imm = mymol.GenerateTopology(aps,
                                 pd,
                                 lot,
                                 iChargeDistributionModel,
                                 bGenVSites,
                                 bPairs,
                                 bDihedral,
                                 false,
                                 tabfn);

    if (immOK == imm)
    {
        maxpot = 100; //Use 100 percent of the ESP read from Gaussian file.
        imm    = mymol.GenerateCharges(pd,
                                       mdlog,
                                       aps,
                                       iChargeDistributionModel,
                                       iChargeGenerationAlgorithm,
                                       watoms,
                                       hfac,
                                       lot,
                                       bQsym,
                                       symm_string,
                                       cr,
                                       tabfn,
                                       nullptr,
                                       qcycle,
                                       maxpot,
                                       qtol,
                                       oenv,
                                       opt2fn_null("-plotESP", NFILE, fnm));
    }
    if (bCUBE && immOK == imm)
    {
        fprintf(stderr, "Fix me: GenerateCube is broken\n");
        mymol.GenerateCube(iChargeDistributionModel,
                           pd,
                           spacing,
                           opt2fn_null("-ref",      NFILE, fnm),
                           opt2fn_null("-pc",       NFILE, fnm),
                           opt2fn_null("-pdbdiff",  NFILE, fnm),
                           opt2fn_null("-pot",      NFILE, fnm),
                           opt2fn_null("-rho",      NFILE, fnm),
                           opt2fn_null("-his",      NFILE, fnm),
                           opt2fn_null("-diff",     NFILE, fnm),
                           opt2fn_null("-diffhist", NFILE, fnm),
                           oenv);
    }

    if (immOK == imm)
    {
        imm = mymol.GenerateChargeGroups(ecg, bUsePDBcharge);
    }

    if (immOK == imm)
    {
        mymol.PrintConformation(opt2fn("-c", NFILE, fnm));
        mymol.PrintTopology(bITP ? ftp2fn(efITP, NFILE, fnm) : ftp2fn(efTOP, NFILE, fnm),
                            iChargeDistributionModel,
                            bVerbose,
                            pd,
                            aps,
                            cr,
                            efield,
                            lot);
    }
    else
    {
        gmx_fatal(FARGS,"Alexandria gentop ended prematurely due to \"%s\"\n", alexandria::immsg(imm));
    }
    return 0;
}
