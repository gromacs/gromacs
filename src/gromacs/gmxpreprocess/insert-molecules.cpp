/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/gmxpreprocess/read-conformation.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/atomsbuilder.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

using gmx::RVec;

/* enum for random rotations of inserted solutes */
enum RotationType {
    en_rotXYZ, en_rotZ, en_rotNone
};
const char *const cRotationEnum[] = {"xyz", "z", "none"};

static void center_molecule(std::vector<RVec> *x)
{
    const size_t atomCount = x->size();
    rvec         center;
    clear_rvec(center);
    for (size_t i = 0; i < atomCount; ++i)
    {
        rvec_inc(center, (*x)[i]);
    }
    svmul(1.0/atomCount, center, center);
    for (size_t i = 0; i < atomCount; ++i)
    {
        rvec_dec((*x)[i], center);
    }
}

static void generate_trial_conf(const std::vector<RVec> &xin,
                                const rvec offset, RotationType enum_rot,
                                gmx::DefaultRandomEngine * rng,
                                std::vector<RVec> *xout)
{
    gmx::UniformRealDistribution<real> dist(0, 2.0*M_PI);
    *xout = xin;

    real alfa = 0.0, beta = 0.0, gamma = 0.0;
    switch (enum_rot)
    {
        case en_rotXYZ:
            alfa  = dist(*rng);
            beta  = dist(*rng);
            gamma = dist(*rng);
            break;
        case en_rotZ:
            alfa  = beta = 0.;
            gamma = dist(*rng);
            break;
        case en_rotNone:
            alfa = beta = gamma = 0.;
            break;
    }
    if (enum_rot == en_rotXYZ || enum_rot == en_rotZ)
    {
        rotate_conf(xout->size(), as_rvec_array(xout->data()), NULL, alfa, beta, gamma);
    }
    for (size_t i = 0; i < xout->size(); ++i)
    {
        rvec_inc((*xout)[i], offset);
    }
}

static bool isInsertionAllowed(gmx::AnalysisNeighborhoodSearch *search,
                               const std::vector<real>         &exclusionDistances,
                               const std::vector<RVec>         &x,
                               const std::vector<real>         &exclusionDistances_insrt,
                               const t_atoms                   &atoms,
                               const std::set<int>             &removableAtoms,
                               gmx::AtomsRemover               *remover)
{
    gmx::AnalysisNeighborhoodPositions  pos(x);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search->startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const real r1 = exclusionDistances[pair.refIndex()];
        const real r2 = exclusionDistances_insrt[pair.testIndex()];
        if (pair.distance2() < gmx::square(r1 + r2))
        {
            if (removableAtoms.count(pair.refIndex()) == 0)
            {
                return false;
            }
            // TODO: If molecule information is available, this should ideally
            // use it to remove whole molecules.
            remover->markResidue(atoms, pair.refIndex(), true);
        }
    }
    return true;
}

static void insert_mols(int nmol_insrt, int ntry, int seed,
                        real defaultDistance, real scaleFactor,
                        t_topology *top, std::vector<RVec> *x,
                        const std::set<int> &removableAtoms,
                        const t_atoms &atoms_insrt, const std::vector<RVec> &x_insrt,
                        int ePBC, matrix box,
                        const std::string &posfn, const rvec deltaR,
                        RotationType enum_rot)
{
    fprintf(stderr, "Initialising inter-atomic distances...\n");
    gmx_atomprop_t          aps = gmx_atomprop_init();
    std::vector<real>       exclusionDistances(
            makeExclusionDistances(&top->atoms, aps, defaultDistance, scaleFactor));
    const std::vector<real> exclusionDistances_insrt(
            makeExclusionDistances(&atoms_insrt, aps, defaultDistance, scaleFactor));
    gmx_atomprop_destroy(aps);

    const real       maxInsertRadius
        = *std::max_element(exclusionDistances_insrt.begin(),
                            exclusionDistances_insrt.end());
    real             maxRadius = maxInsertRadius;
    if (!exclusionDistances.empty())
    {
        const real maxExistingRadius
            = *std::max_element(exclusionDistances.begin(),
                                exclusionDistances.end());
        maxRadius = std::max(maxInsertRadius, maxExistingRadius);
    }

    // TODO: Make all of this exception-safe.
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(maxInsertRadius + maxRadius);


    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
    }
    fprintf(stderr, "Using random seed %d\n", seed);

    gmx::DefaultRandomEngine rng(seed);

    t_pbc                    pbc;
    set_pbc(&pbc, ePBC, box);

    /* With -ip, take nmol_insrt from file posfn */
    double     **rpos              = NULL;
    const bool   insertAtPositions = !posfn.empty();
    if (insertAtPositions)
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

    gmx::AtomsBuilder builder(&top->atoms, &top->symtab);
    gmx::AtomsRemover remover(top->atoms);
    {
        const int finalAtomCount    = top->atoms.nr + nmol_insrt * atoms_insrt.nr;
        const int finalResidueCount = top->atoms.nres + nmol_insrt * atoms_insrt.nres;
        builder.reserve(finalAtomCount, finalResidueCount);
        x->reserve(finalAtomCount);
        exclusionDistances.reserve(finalAtomCount);
    }

    std::vector<RVec>                    x_n(x_insrt.size());

    int                                  mol        = 0;
    int                                  trial      = 0;
    int                                  firstTrial = 0;
    int                                  failed     = 0;
    gmx::UniformRealDistribution<real>   dist;

    while (mol < nmol_insrt && trial < ntry*nmol_insrt)
    {
        // cppcheck 1.72 complains about uninitialized variables in the
        // assignments below otherwise...
        rvec offset_x = {0};
        if (!insertAtPositions)
        {
            // Insert at random positions.
            offset_x[XX] = box[XX][XX] * dist(rng);
            offset_x[YY] = box[YY][YY] * dist(rng);
            offset_x[ZZ] = box[ZZ][ZZ] * dist(rng);
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
            offset_x[XX] = rpos[XX][mol] + deltaR[XX]*(2 * dist(rng)-1);
            offset_x[YY] = rpos[YY][mol] + deltaR[YY]*(2 * dist(rng)-1);
            offset_x[ZZ] = rpos[ZZ][mol] + deltaR[ZZ]*(2 * dist(rng)-1);
        }
        fprintf(stderr, "\rTry %d", ++trial);
        fflush(stderr);

        generate_trial_conf(x_insrt, offset_x, enum_rot, &rng, &x_n);
        gmx::AnalysisNeighborhoodPositions pos(*x);
        gmx::AnalysisNeighborhoodSearch    search = nb.initSearch(&pbc, pos);
        if (isInsertionAllowed(&search, exclusionDistances, x_n, exclusionDistances_insrt,
                               top->atoms, removableAtoms, &remover))
        {
            x->insert(x->end(), x_n.begin(), x_n.end());
            exclusionDistances.insert(exclusionDistances.end(),
                                      exclusionDistances_insrt.begin(),
                                      exclusionDistances_insrt.end());
            builder.mergeAtoms(atoms_insrt);
            ++mol;
            firstTrial = trial;
            fprintf(stderr, " success (now %d atoms)!\n", builder.currentAtomCount());
        }
    }

    fprintf(stderr, "\n");
    /* print number of molecules added */
    fprintf(stderr, "Added %d molecules (out of %d requested)\n",
            mol - failed, nmol_insrt);

    const int originalAtomCount    = top->atoms.nr;
    const int originalResidueCount = top->atoms.nres;
    remover.refreshAtomCount(top->atoms);
    remover.removeMarkedElements(x);
    remover.removeMarkedAtoms(&top->atoms);
    if (top->atoms.nr < originalAtomCount)
    {
        fprintf(stderr, "Replaced %d residues (%d atoms)\n",
                originalResidueCount - top->atoms.nres,
                originalAtomCount - top->atoms.nr);
    }

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

class InsertMolecules : public ICommandLineOptionsModule, public ITopologyProvider
{
    public:
        InsertMolecules()
            : bBox_(false), nmolIns_(0), nmolTry_(10), seed_(0),
              defaultDistance_(0.105), scaleFactor_(0.57), enumRot_(en_rotXYZ),
              top_(NULL), ePBC_(-1)
        {
            clear_rvec(newBox_);
            clear_rvec(deltaR_);
            clear_mat(box_);
        }
        virtual ~InsertMolecules()
        {
            if (top_ != NULL)
            {
                done_top(top_);
                sfree(top_);
            }
        }

        // From ITopologyProvider
        virtual t_topology *getTopology(bool /*required*/) { return top_; }
        virtual int getAtomCount() { return 0; }

        // From ICommandLineOptionsModule
        virtual void init(CommandLineModuleSettings * /*settings*/)
        {
        }
        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings);
        virtual void optionsFinished();
        virtual int run();

    private:
        void loadSolute();

        SelectionCollection selections_;

        std::string         inputConfFile_;
        std::string         insertConfFile_;
        std::string         positionFile_;
        std::string         outputConfFile_;
        rvec                newBox_;
        bool                bBox_;
        int                 nmolIns_;
        int                 nmolTry_;
        int                 seed_;
        real                defaultDistance_;
        real                scaleFactor_;
        rvec                deltaR_;
        RotationType        enumRot_;
        Selection           replaceSel_;

        t_topology         *top_;
        std::vector<RVec>   x_;
        matrix              box_;
        int                 ePBC_;
};

void InsertMolecules::initOptions(IOptionsContainer                 *options,
                                  ICommandLineOptionsModuleSettings *settings)
{
    const char *const desc[] = {
        "[THISMODULE] inserts [TT]-nmol[tt] copies of the system specified in",
        "the [TT]-ci[tt] input file. The insertions take place either into",
        "vacant space in the solute conformation given with [TT]-f[tt], or",
        "into an empty box given by [TT]-box[tt]. Specifying both [TT]-f[tt]",
        "and [TT]-box[tt] behaves like [TT]-f[tt], but places a new box",
        "around the solute before insertions. Any velocities present are",
        "discarded.",
        "",
        "It is possible to also insert into a solvated configuration and",
        "replace solvent atoms with the inserted atoms. To do this, use",
        "[TT]-replace[tt] to specify a selection that identifies the atoms",
        "that can be replaced. The tool assumes that all molecules in this",
        "selection consist of single residues: each residue from this",
        "selection that overlaps with the inserted molecules will be removed",
        "instead of preventing insertion.",
        "",
        "By default, the insertion positions are random (with initial seed",
        "specified by [TT]-seed[tt]). The program iterates until [TT]-nmol[tt]",
        "molecules have been inserted in the box. Molecules are not inserted",
        "where the distance between any existing atom and any atom of the",
        "inserted molecule is less than the sum based on the van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those"
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].",
        "",
        "A total of [TT]-nmol[tt] * [TT]-try[tt] insertion attempts are made",
        "before giving up. Increase [TT]-try[tt] if you have several small",
        "holes to fill. Option [TT]-rot[tt] specifies whether the insertion",
        "molecules are randomly oriented before insertion attempts.",
        "",
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

    settings->setHelpText(desc);

    std::shared_ptr<SelectionOptionBehavior> selectionOptionBehavior(
            new SelectionOptionBehavior(&selections_, this));
    settings->addOptionsBehavior(selectionOptionBehavior);

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

    options->addOption(SelectionOption("replace").onlyAtoms()
                           .store(&replaceSel_)
                           .description("Atoms that can be removed if overlapping"));
    selectionOptionBehavior->initOptions(options);

    options->addOption(RealOption("box").vector()
                           .store(newBox_).storeIsSet(&bBox_)
                           .description("Box size (in nm)"));
    options->addOption(IntegerOption("nmol")
                           .store(&nmolIns_)
                           .description("Number of extra molecules to insert"));
    options->addOption(IntegerOption("try")
                           .store(&nmolTry_)
                           .description("Try inserting [TT]-nmol[tt] times [TT]-try[tt] times"));
    options->addOption(IntegerOption("seed")
                           .store(&seed_)
                           .description("Random generator seed (0 means generate)"));
    options->addOption(RealOption("radius")
                           .store(&defaultDistance_)
                           .description("Default van der Waals distance"));
    options->addOption(RealOption("scale")
                           .store(&scaleFactor_)
                           .description("Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water."));
    options->addOption(RealOption("dr").vector()
                           .store(deltaR_)
                           .description("Allowed displacement in x/y/z from positions in [TT]-ip[tt] file"));
    options->addOption(EnumOption<RotationType>("rot").enumValue(cRotationEnum)
                           .store(&enumRot_)
                           .description("Rotate inserted molecules randomly"));
}

void InsertMolecules::optionsFinished()
{
    if (nmolIns_ <= 0 && positionFile_.empty())
    {
        GMX_THROW(InconsistentInputError("Either -nmol must be larger than 0, "
                                         "or positions must be given with -ip."));
    }
    if (inputConfFile_.empty() && !bBox_)
    {
        GMX_THROW(InconsistentInputError("When no solute (-f) is specified, "
                                         "a box size (-box) must be specified."));
    }
    if (replaceSel_.isValid() && inputConfFile_.empty())
    {
        GMX_THROW(InconsistentInputError("Replacement (-replace) only makes sense "
                                         "together with an existing configuration (-f)."));
    }

    snew(top_, 1);
    if (!inputConfFile_.empty())
    {
        readConformation(inputConfFile_.c_str(), top_, &x_, NULL,
                         &ePBC_, box_, "solute");
        if (top_->atoms.nr == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", inputConfFile_.c_str());
        }
    }
}

int InsertMolecules::run()
{
    std::set<int> removableAtoms;
    if (replaceSel_.isValid())
    {
        t_pbc       pbc;
        set_pbc(&pbc, ePBC_, box_);
        t_trxframe *fr;
        snew(fr, 1);
        fr->natoms = x_.size();
        fr->bX     = TRUE;
        fr->x      = as_rvec_array(x_.data());
        selections_.evaluate(fr, &pbc);
        sfree(fr);
        removableAtoms.insert(replaceSel_.atomIndices().begin(),
                              replaceSel_.atomIndices().end());
        // TODO: It could be nice to check that removableAtoms contains full
        // residues, since we anyways remove whole residues instead of
        // individual atoms.
    }

    if (bBox_)
    {
        ePBC_ = epbcXYZ;
        clear_mat(box_);
        box_[XX][XX] = newBox_[XX];
        box_[YY][YY] = newBox_[YY];
        box_[ZZ][ZZ] = newBox_[ZZ];
    }
    if (det(box_) == 0)
    {
        gmx_fatal(FARGS, "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }

    t_topology        *top_insrt;
    std::vector<RVec>  x_insrt;
    snew(top_insrt, 1);
    {
        int         ePBC_dummy;
        matrix      box_dummy;
        readConformation(insertConfFile_.c_str(), top_insrt, &x_insrt,
                         NULL, &ePBC_dummy, box_dummy, "molecule");
        if (top_insrt->atoms.nr == 0)
        {
            gmx_fatal(FARGS, "No molecule in %s, please check your input",
                      insertConfFile_.c_str());
        }
        if (top_->name == NULL)
        {
            top_->name = top_insrt->name;
        }
        if (positionFile_.empty())
        {
            center_molecule(&x_insrt);
        }
    }

    /* add nmol_ins molecules of atoms_ins
       in random orientation at random place */
    insert_mols(nmolIns_, nmolTry_, seed_, defaultDistance_, scaleFactor_,
                top_, &x_, removableAtoms, top_insrt->atoms, x_insrt,
                ePBC_, box_, positionFile_, deltaR_, enumRot_);

    /* write new configuration to file confout */
    fprintf(stderr, "Writing generated configuration to %s\n",
            outputConfFile_.c_str());
    write_sto_conf(outputConfFile_.c_str(), *top_->name, &top_->atoms,
                   as_rvec_array(x_.data()), NULL, ePBC_, box_);

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n",
            top_->atoms.nr, top_->atoms.nres);

    done_top(top_insrt);
    sfree(top_insrt);

    return 0;
}

}   // namespace

const char InsertMoleculesInfo::name[]             = "insert-molecules";
const char InsertMoleculesInfo::shortDescription[] =
    "Insert molecules into existing vacancies";
ICommandLineOptionsModulePointer InsertMoleculesInfo::create()
{
    return ICommandLineOptionsModulePointer(new InsertMolecules());
}

} // namespace gmx
