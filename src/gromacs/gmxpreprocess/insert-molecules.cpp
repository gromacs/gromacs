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

#include "insert-molecules.h"

#include <algorithm>
#include <string>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/gmxpreprocess/read-conformation.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* enum for random rotations of inserted solutes */
enum {
    en_rotXYZ, en_rotZ, en_rotNone
};

static void center_molecule(int atomCount, rvec x[])
{
    rvec center;
    clear_rvec(center);
    for (int i = 0; i < atomCount; ++i)
    {
        rvec_inc(center, x[i]);
    }
    svmul(1.0/atomCount, center, center);
    for (int i = 0; i < atomCount; ++i)
    {
        rvec_dec(x[i], center);
    }
}

static void generate_trial_conf(int atomCount, const rvec xin[],
                                const rvec offset, int enum_rot, gmx_rng_t rng,
                                rvec xout[])
{
    for (int i = 0; i < atomCount; ++i)
    {
        copy_rvec(xin[i], xout[i]);
    }
    real alfa = 0.0, beta = 0.0, gamma = 0.0;
    switch (enum_rot)
    {
        case en_rotXYZ:
            alfa  = 2*M_PI * gmx_rng_uniform_real(rng);
            beta  = 2*M_PI * gmx_rng_uniform_real(rng);
            gamma = 2*M_PI * gmx_rng_uniform_real(rng);
            break;
        case en_rotZ:
            alfa  = beta = 0.;
            gamma = 2*M_PI * gmx_rng_uniform_real(rng);
            break;
        case en_rotNone:
            alfa = beta = gamma = 0.;
            break;
    }
    if (enum_rot == en_rotXYZ || (enum_rot == en_rotZ))
    {
        rotate_conf(atomCount, xout, NULL, alfa, beta, gamma);
    }
    for (int i = 0; i < atomCount; ++i)
    {
        rvec_inc(xout[i], offset);
    }
}

static bool is_insertion_allowed(gmx::AnalysisNeighborhoodSearch *search,
                                 const real *exclusionDistances,
                                 int atomCount, const rvec *x,
                                 const real *exclusionDistances_insrt)
{
    gmx::AnalysisNeighborhoodPositions  pos(x, atomCount);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search->startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const real r1 = exclusionDistances[pair.refIndex()];
        const real r2 = exclusionDistances_insrt[pair.testIndex()];
        if (pair.distance2() < sqr(r1 + r2))
        {
            return false;
        }
    }
    return true;
}

static void merge_atoms_noalloc(t_atoms *atoms, const t_atoms *atoms_add)
{
    int resnr = 0;
    if (atoms->nr > 0)
    {
        resnr = atoms->resinfo[atoms->atom[atoms->nr-1].resind].nr;
    }
    int prevResInd = -1;
    for (int i = 0; i < atoms_add->nr; ++i)
    {
        if (atoms_add->atom[i].resind != prevResInd)
        {
            prevResInd = atoms_add->atom[i].resind;
            ++resnr;
            atoms->resinfo[atoms->nres]    = atoms_add->resinfo[prevResInd];
            atoms->resinfo[atoms->nres].nr = resnr;
            ++atoms->nres;
        }
        atoms->atom[atoms->nr]        = atoms_add->atom[i];
        atoms->atomname[atoms->nr]    = atoms_add->atomname[i];
        atoms->atom[atoms->nr].resind = atoms->nres-1;
        ++atoms->nr;
    }
}

static void insert_mols(int nmol_insrt, int ntry, int seed,
                        real defaultDistance, real scaleFactor,
                        t_atoms *atoms, rvec **x,
                        const t_atoms *atoms_insrt, const rvec *x_insrt,
                        int ePBC, matrix box,
                        const std::string &posfn, const rvec deltaR, int enum_rot)
{
    t_pbc            pbc;
    rvec            *x_n;

    fprintf(stderr, "Initialising inter-atomic distances...\n");
    gmx_atomprop_t   aps = gmx_atomprop_init();
    real            *exclusionDistances
        = makeExclusionDistances(atoms, aps, defaultDistance, scaleFactor);
    real            *exclusionDistances_insrt
        = makeExclusionDistances(atoms_insrt, aps, defaultDistance, scaleFactor);
    gmx_atomprop_destroy(aps);

    const real       maxInsertRadius
        = *std::max_element(exclusionDistances_insrt,
                            exclusionDistances_insrt + atoms_insrt->nr);
    real             maxRadius = maxInsertRadius;
    if (atoms->nr > 0)
    {
        const real maxExistingRadius
            = *std::max_element(exclusionDistances,
                                exclusionDistances + atoms->nr);
        maxRadius = std::max(maxInsertRadius, maxExistingRadius);
    }

    // TODO: Make all of this exception-safe.
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(maxInsertRadius + maxRadius);

    gmx_rng_t        rng = gmx_rng_init(seed);
    set_pbc(&pbc, ePBC, box);

    snew(x_n, atoms_insrt->nr);

    /* With -ip, take nmol_insrt from file posfn */
    double         **rpos = NULL;
    if (!posfn.empty())
    {
        int ncol;
        nmol_insrt = read_xvg(posfn.c_str(), &rpos, &ncol);
        if (ncol != 3)
        {
            gmx_fatal(FARGS, "Expected 3 columns (x/y/z coordinates) in file %s\n",
                      posfn.c_str());
        }
        fprintf(stderr, "Read %d positions from file %s\n\n",
                nmol_insrt, posfn.c_str());
    }

    {
        const int finalAtomCount    = atoms->nr + nmol_insrt * atoms_insrt->nr;
        const int finalResidueCount = atoms->nres + nmol_insrt * atoms_insrt->nres;
        srenew(atoms->resinfo,      finalResidueCount);
        srenew(atoms->atomname,     finalAtomCount);
        srenew(atoms->atom,         finalAtomCount);
        srenew(*x,                  finalAtomCount);
        srenew(exclusionDistances,  finalAtomCount);
    }

    int mol        = 0;
    int trial      = 0;
    int firstTrial = 0;
    int failed     = 0;
    while ((mol < nmol_insrt) && (trial < ntry*nmol_insrt))
    {
        rvec offset_x;
        if (posfn.empty())
        {
            // Insert at random positions.
            offset_x[XX] = box[XX][XX] * gmx_rng_uniform_real(rng);
            offset_x[YY] = box[YY][YY] * gmx_rng_uniform_real(rng);
            offset_x[ZZ] = box[ZZ][ZZ] * gmx_rng_uniform_real(rng);
        }
        else
        {
            // Skip a position if ntry trials were not successful.
            if (trial >= firstTrial + ntry)
            {
                fprintf(stderr, " skipped position (%.3f, %.3f, %.3f)\n",
                        rpos[XX][mol], rpos[YY][mol], rpos[ZZ][mol]);
                ++mol;
                ++failed;
            }
            // Insert at positions taken from option -ip file.
            offset_x[XX] = rpos[XX][mol] + deltaR[XX]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[YY] = rpos[YY][mol] + deltaR[YY]*(2 * gmx_rng_uniform_real(rng)-1);
            offset_x[ZZ] = rpos[ZZ][mol] + deltaR[ZZ]*(2 * gmx_rng_uniform_real(rng)-1);
        }
        fprintf(stderr, "\rTry %d", ++trial);
        generate_trial_conf(atoms_insrt->nr, x_insrt, offset_x, enum_rot, rng,
                            x_n);
        gmx::AnalysisNeighborhoodPositions pos(*x, atoms->nr);
        gmx::AnalysisNeighborhoodSearch    search = nb.initSearch(&pbc, pos);
        if (is_insertion_allowed(&search, exclusionDistances, atoms_insrt->nr,
                                 x_n, exclusionDistances_insrt))
        {
            const int firstIndex = atoms->nr;
            for (int i = 0; i < atoms_insrt->nr; ++i)
            {
                copy_rvec(x_n[i], (*x)[firstIndex + i]);
                exclusionDistances[firstIndex + i] = exclusionDistances_insrt[i];
            }
            merge_atoms_noalloc(atoms, atoms_insrt);
            ++mol;
            firstTrial = trial;
            fprintf(stderr, " success (now %d atoms)!\n", atoms->nr);
        }
    }
    gmx_rng_destroy(rng);
    srenew(atoms->resinfo,  atoms->nres);
    srenew(atoms->atomname, atoms->nr);
    srenew(atoms->atom,     atoms->nr);
    srenew(*x,              atoms->nr);

    fprintf(stderr, "\n");
    /* print number of molecules added */
    fprintf(stderr, "Added %d molecules (out of %d requested)\n",
            mol - failed, nmol_insrt);

    sfree(x_n);
    sfree(exclusionDistances);
    sfree(exclusionDistances_insrt);
    if (rpos != NULL)
    {
        for (int i = 0; i < DIM; ++i)
        {
            sfree(rpos[i]);
        }
        sfree(rpos);
    }
}

namespace gmx
{

namespace
{

class InsertMolecules : public CommandLineOptionsModuleInterface
{
    public:
        InsertMolecules()
            : bBox_(false), nmolIns_(0), nmolTry_(10), seed_(1997),
              defaultDistance_(0.105), scaleFactor_(0.57), enumRot_(en_rotXYZ)
        {
            clear_rvec(newBox_);
            clear_rvec(deltaR_);
        }

        virtual void init(CommandLineModuleSettings * /*settings*/)
        {
        }

        virtual void initOptions(Options *options);
        virtual void optionsFinished(Options *options);

        virtual int run();

    private:
        std::string inputConfFile_;
        std::string insertConfFile_;
        std::string positionFile_;
        std::string outputConfFile_;
        rvec        newBox_;
        bool        bBox_;
        int         nmolIns_;
        int         nmolTry_;
        int         seed_;
        real        defaultDistance_;
        real        scaleFactor_;
        rvec        deltaR_;
        int         enumRot_;
};

void InsertMolecules::initOptions(Options *options)
{
    const char *const desc[] = {
        "[THISMODULE] inserts [TT]-nmol[tt] copies of the system specified in",
        "the [TT]-ci[tt] input file. The insertions take place either into",
        "vacant space in the solute conformation given with [TT]-f[tt], or",
        "into an empty box given by [TT]-box[tt]. Specifying both [TT]-f[tt]",
        "and [TT]-box[tt] behaves like [TT]-f[tt], but places a new box",
        "around the solute before insertions. Any velocities present are",
        "discarded.[PAR]",

        "By default, the insertion positions are random (with initial seed",
        "specified by [TT]-seed[tt]). The program iterates until [TT]-nmol[tt]",
        "molecules have been inserted in the box. Molecules are not inserted",
        "where the distance between any existing atom and any atom of the",
        "inserted molecule is less than the sum based on the van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those"
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].[PAR]",

        "A total of [TT]-nmol[tt] * [TT]-try[tt] insertion attempts are made",
        "before giving up. Increase [TT]-try[tt] if you have several small",
        "holes to fill. Option [TT]-rot[tt] specifies whether the insertion",
        "molecules are randomly oriented before insertion attempts.[PAR]",

        "Alternatively, the molecules can be inserted only at positions defined in",
        "positions.dat ([TT]-ip[tt]). That file should have 3 columns (x,y,z),",
        "that give the displacements compared to the input molecule position",
        "([TT]-ci[tt]). Hence, if that file should contain the absolute",
        "positions, the molecule must be centered on (0,0,0) before using",
        "[THISMODULE] (e.g. from [gmx-editconf] [TT]-center[tt]).",
        "Comments in that file starting with # are ignored. Option [TT]-dr[tt]",
        "defines the maximally allowed displacements during insertial trials.",
        "[TT]-try[tt] and [TT]-rot[tt] work as in the default mode (see above)."
    };

    options->setDescription(desc);

    // TODO: Replace use of legacyType.
    options->addOption(FileNameOption("f")
                           .legacyType(efSTX).inputFile()
                           .store(&inputConfFile_)
                           .defaultBasename("protein")
                           .description("Existing configuration to insert into"));
    options->addOption(FileNameOption("ci")
                           .legacyType(efSTX).inputFile().required()
                           .store(&insertConfFile_)
                           .defaultBasename("insert")
                           .description("Configuration to insert"));
    options->addOption(FileNameOption("ip")
                           .filetype(eftGenericData).inputFile()
                           .store(&positionFile_)
                           .defaultBasename("positions")
                           .description("Predefined insertion trial positions"));
    options->addOption(FileNameOption("o")
                           .legacyType(efSTO).outputFile().required()
                           .store(&outputConfFile_)
                           .defaultBasename("out")
                           .description("Output configuration after insertion"));

    options->addOption(RealOption("box").vector()
                           .store(newBox_)
                           .description("Box size (in nm)"));
    options->addOption(IntegerOption("nmol")
                           .store(&nmolIns_)
                           .description("Number of extra molecules to insert"));
    options->addOption(IntegerOption("try")
                           .store(&nmolTry_)
                           .description("Try inserting [TT]-nmol[tt] times [TT]-try[tt] times"));
    options->addOption(IntegerOption("seed")
                           .store(&seed_)
                           .description("Random generator seed"));
    options->addOption(RealOption("radius")
                           .store(&defaultDistance_)
                           .description("Default van der Waals distance"));
    options->addOption(RealOption("scale")
                           .store(&scaleFactor_)
                           .description("Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water."));
    options->addOption(RealOption("dr").vector()
                           .store(deltaR_)
                           .description("Allowed displacement in x/y/z from positions in [TT]-ip[tt] file"));
    const char *const cRotationEnum[] = {"xyz", "z", "none"};
    options->addOption(StringOption("rot").enumValue(cRotationEnum)
                           .storeEnumIndex(&enumRot_)
                           .description("Rotate inserted molecules randomly"));
}

void InsertMolecules::optionsFinished(Options *options)
{
    bBox_ = options->isSet("box");
}

int InsertMolecules::run()
{
    const bool bProt = !inputConfFile_.empty();

    /* check input */
    if (nmolIns_ <= 0 && positionFile_.empty())
    {
        gmx_fatal(FARGS, "Either -nmol must be larger than 0, "
                  "or positions must be given with -ip");
    }
    if (!bProt && !bBox_)
    {
        gmx_fatal(FARGS, "When no solute (-f) is specified, "
                  "a box size (-box) must be specified");
    }

    char    *title = NULL;
    t_atoms *atoms;
    rvec    *x = NULL;
    matrix   box;
    int      ePBC = -1;
    snew(atoms, 1);
    init_t_atoms(atoms, 0, FALSE);
    if (bProt)
    {
        /* Generate a solute configuration */
        title = readConformation(inputConfFile_.c_str(), atoms, &x, NULL,
                                 &ePBC, box, "solute");
        if (atoms->nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", inputConfFile_.c_str());
            sfree(title);
            title = NULL;
        }
    }
    if (bBox_)
    {
        ePBC = epbcXYZ;
        clear_mat(box);
        box[XX][XX] = newBox_[XX];
        box[YY][YY] = newBox_[YY];
        box[ZZ][ZZ] = newBox_[ZZ];
    }
    if (det(box) == 0)
    {
        gmx_fatal(FARGS, "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }

    t_atoms *atoms_insrt;
    rvec    *x_insrt = NULL;
    snew(atoms_insrt, 1);
    init_t_atoms(atoms_insrt, 0, FALSE);
    {
        int         ePBC_dummy;
        matrix      box_dummy;
        char       *title_ins
            = readConformation(insertConfFile_.c_str(), atoms_insrt, &x_insrt,
                               NULL, &ePBC_dummy, box_dummy, "molecule");
        if (atoms_insrt->nr == 0)
        {
            gmx_fatal(FARGS, "No molecule in %s, please check your input",
                      insertConfFile_.c_str());
        }
        if (title == NULL)
        {
            title = title_ins;
        }
        else
        {
            sfree(title_ins);
        }
        if (positionFile_.empty())
        {
            center_molecule(atoms_insrt->nr, x_insrt);
        }
    }

    /* add nmol_ins molecules of atoms_ins
       in random orientation at random place */
    insert_mols(nmolIns_, nmolTry_, seed_, defaultDistance_, scaleFactor_,
                atoms, &x, atoms_insrt, x_insrt,
                ePBC, box, positionFile_, deltaR_, enumRot_);

    /* write new configuration to file confout */
    fprintf(stderr, "Writing generated configuration to %s\n",
            outputConfFile_.c_str());
    write_sto_conf(outputConfFile_.c_str(), title, atoms, x, NULL, ePBC, box);

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            atoms->nr, atoms->nres);

    sfree(x);
    sfree(x_insrt);
    done_atom(atoms);
    done_atom(atoms_insrt);
    sfree(atoms);
    sfree(atoms_insrt);
    sfree(title);

    return 0;
}

}   // namespace

const char InsertMoleculesInfo::name[]             = "insert-molecules";
const char InsertMoleculesInfo::shortDescription[] =
    "Insert molecules into existing vacancies";
CommandLineOptionsModuleInterface *InsertMoleculesInfo::create()
{
    return new InsertMolecules();
}

} // namespace gmx
