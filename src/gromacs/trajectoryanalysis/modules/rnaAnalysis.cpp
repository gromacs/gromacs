/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::RnaAnalysis.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders G�rden�s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "rnaAnalysis.h"

#include "gromacs/options.h"

namespace gmx
{

namespace analysismodules
{


// Constructor
RnaAnalysis::RnaAnalysis() : adata_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(0, 1);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "rnaAnalysis");

    offsetAtom         = 10; // probably be C1'.
    hydrogenRmsd_      = true;
    sugarRmsd_         = true;
    phosphateRmsd_     = true;
    printAllBasePairs_ = false;
    statistic          = false;
    oneBPList          = true;
    Xpm                = false;
    extraOffset        = 0.0;
    addRMSD_           = 0.0;
    datailedInfo_      = false;
    iAtoms             = NULL;
    bondDist_          = 0.0;
    atomprop_          = gmx_atomprop_init();
    rnaLog             = NULL;
    massSugar_         = 0.0;
}

// Destructor
RnaAnalysis::~RnaAnalysis()
{
    // Destroy C structures where there is no automatic memory release
    // C++ takes care of memory in classes (hopefully)
    gmx_atomprop_destroy(atomprop_);
    for (auto &ri : nucleotideInfo_)
    {
        delete ri;
    }
    for (auto &temp : templateDB_)
    {
        delete temp;
    }
}

void
RnaAnalysis::initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Please note that we are still working on improving this tool!\n",
        "Reads and analyzes an RNA structure (-s [<.gro/..>]) or an MD simulation trajectory (-f [<.xtc/.trr/..>] -s [<.gro/..>])",
        "to detect all bases that are in contact (e.g., form base pairs).\n\n",
        "This is can be done in two ways:\n",
        "1. Check if hydrogen bonds between two RNA bases in the structure match a certain distance criteria (-templateDB hbonds).\n",
        "2. compare RNA structure bases to a certain base pair template database (-templateDB bps).\n\n",
        "The data is displayed on the terminal as raw data. More information about",
        "hydrogen bond distances and RMSD to template structures can be displayed(-moreInfo true).\n",
        "The data can also printed to an xpm file can be coverted to an eps file displaying a matrix",
        "describing all base pairs (-matrix true).\n",
        "If you want, you can also store all detected base pairs to pdb files (-pdbFile yes),",
        "you can also specify the directory in which the pdb files are stored (-outPath <string>)."

    };

    // Add the descriptive text (program help text) to the options
    settings->setHelpText(ConstArrayRef<const char *>::fromArray(desc, 12));

    // Add option for optional output file
    options->addOption(FileNameOption("g").filetype(eftGenericData).outputFile()
                           .store(&rnaLogFile).defaultBasename("result")
                           .description("output file").required());

    // Add option for optional plot file
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&rnaGraphFile).defaultBasename("result")
                           .description("2D representation file").required());

    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select").valueCount(1)
                           .store(&sel_).defaultSelectionText("all").onlyAtoms()
                           .description("Use default (all); select which atoms should be used for analysis"));

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);

    // Template databases in which base pairs are stored
    options->addOption(StringOption("templateDB").store(&DB).defaultValue("bps")
                           .description("define the template DB you want to use (currently just bps available)"));

    // Print all results about RNA base pairs during analysis
    options->addOption(BooleanOption("printAllBasePairs").store(&printAllBasePairs_).defaultValue(false)
                           .description("print all possible base pairs during analysis"));

    // Include hydrogen atoms in RMSD calculation
    options->addOption(BooleanOption("hydro").store(&hydrogenRmsd_).defaultValue(true)
                           .description("include the hydrogen atoms of the nucleotide in RMSD calculations"));

    // Include  phosphate group atoms in RMSD calculation
    options->addOption(BooleanOption("phos").store(&phosphateRmsd_).defaultValue(true)
                           .description("include the phosphate atoms of the nucleotide in RMSD calculation"));

    // Include ONLY base atoms of nucleotides in RMSD calculations
    options->addOption(BooleanOption("sugar").store(&sugarRmsd_).defaultValue(true)
                           .description("include the sugar atoms of the nucleotide in RMSD calculation"));

    // Distance of hydrogen bonds
    options->addOption(DoubleOption("bondDist").store(&bondDist_).defaultValue(0.33)
                           .description("distance between hydrogen bond donor and acceptor atom"));

    // Scale for mass of sugar.
    options->addOption(DoubleOption("massSuagr").store(&massSugar_).defaultValue(0.05)
                           .description("scale for mass of sugar when calculating RMSD"));

    // Scale for mass of phosphate.
    options->addOption(DoubleOption("massPhosphate").store(&massPhosphate_).defaultValue(0.05)
                           .description("scale for mass of phosphate when calculating RMSD"));


    // Maximum RMSD for between identified base pair and template base pair
    options->addOption(DoubleOption("addRMSD").store(&addRMSD_).defaultValue(0.0)
                           .description(
                               "increase the maximum RMSD (cut-off value) by offset for identifying base pairs (default = 0.0)"));

    // More information about output data
    options->addOption(BooleanOption("detailedInfo").store(&datailedInfo_).defaultValue(false)
                           .description(
                               "details of H-bond distances and RMSD values between identified and template base pairs"));


    // Extra offset used when looking for base pairs in close proximity
    options->addOption(DoubleOption("addOffset").store(&extraOffset).defaultValue(0.5)
                           .description(
                               "increase distance by offset when searching for atoms within a certain range"));

    // Matrix output file
    options->addOption(BooleanOption("matrix").store(&Xpm).defaultValue(false)
                           .description(
                               "matrix xpm output file is printed, one file for each frame of the trajectory (so far only for templateDB bps)"));

    // Saves coordinates of base pairs in PDB files
    options->addOption(BooleanOption("pdbFiles").store(&statistic).defaultValue(false)
                           .description("save coordinates of all identified RNA base pairs in PDB files"));

    // Output path to directory in which PDB files are stored
    options->addOption(StringOption("outPath").store(&outPath).defaultValue("")
                           .description(
                               "output path to directory which already exists (default current directory)"));

}

//! Check if a ResInfo is valid (by the type and the number of atoms)
static bool validateResInfo(const ResInfo &ri)
{
    GMX_RELEASE_ASSERT((ri.residueType() != Thymine), "Thymine not supported!");
    if (NULL != debug)
    {
        fprintf(debug, "residueType %c, nAtoms %zu\n", ri.residueType(), ri.nAtoms());
    }
    int dn = (int) numberOfNucleotideAtoms[ri.residueType()] - (int) ri.nAtoms();
    return (dn == 0 ||
            (ri.endDirectionality() == FivePrimeEnd && dn == 2) ||
            (ri.endDirectionality() == ThreePrimeEnd && dn == -1));
}

void RnaAnalysis::readTemplates()
{
    // Store all file names of RNA base pair templates within input directory
    // in vector pdbTemplates
    std::string dbName = "RNA-" + DB + ".pdb";
    std::string str1   = "amber99sb-ildn.ff/" + dbName;
    rvec        x;

    // Set path and load file pointer
    real maxDist = 0;

    // Initialize RNA template data base to compare given RNA structure to
    BasePair  *bp = NULL;
    char      *fn = gmxlibfn(str1.c_str());
    TextReader fp {
        std::string(fn)
    };
    sfree(fn);
    std::string line;
    while (fp.readLine(&line))
    {
        if (NULL == bp)
        {
            bp = new BasePair(hydrogenRmsd_, sugarRmsd_, phosphateRmsd_);
        }
        std::vector<std::string> elements = splitString(line);
        if (elements.size() == 0)
        {
            continue;
        }
        int         nElem   = 0;
        std::string pdbType = elements[nElem++];

        if (pdbType.compare("MODEL") == 0)
        {
            continue;
        }
        else if ((pdbType.compare("REMARK") == 0) && (elements.size() >= 3))
        {
            bp->setNucleotides(elements[nElem]);
            bp->setBondTypes(elements[nElem++]);
            bp->setBondSubtypes(elements[nElem++]);
            bp->setTemplateScore(std::strtod(elements[nElem++].c_str(), NULL));
            bp->setModeled(elements[nElem++].compare("1") == 0);
        }
        else if ((pdbType.compare("CRYST1") == 0) && (elements.size() >= 7))
        {
            matrix box;
            int    ePBC;
            read_cryst1(line.c_str(), &ePBC, box);
            bp->setPBC(ePBC, box);
        }
        else if ((pdbType.compare("ATOM") == 0) &&
                 (elements.size() >= 8))
        {
            size_t      anum  = std::atoi(elements[nElem++].c_str());
            std::string aname = elements[nElem++];
            std::string rname = elements[nElem++];
            size_t      rnum  = std::atoi(elements[nElem++].c_str());
            // Convert coordinates from �ngstr�m to nanometer
            x[XX] = 0.1 * std::strtod(elements[nElem++].c_str(), NULL);
            x[YY] = 0.1 * std::strtod(elements[nElem++].c_str(), NULL);
            x[ZZ] = 0.1 * std::strtod(elements[nElem++].c_str(), NULL);
            real m;
            if (gmx_atomprop_query(atomprop_, epropMass, "???",
                                   aname.c_str(), &m))
            {
                m = getMass(m, aname);
                bp->addBasePairAtom(x, anum, aname, rname, rnum, m);
            }
            else
            {
                GMX_THROW(InvalidInputError(line));
            }

        }
        else if ((pdbType.compare("CONECT") == 0) &&
                 (elements.size() >= 3))
        {
            int ai = std::atoi(elements[nElem++].c_str());
            int aj = std::atoi(elements[nElem++].c_str());
            bp->addBondAtoms(ai, aj);
        }
        else if (pdbType.compare("ENDMDL") == 0)
        {
            // Assign the center of mass to origin of all atoms
            // in three dimensions.
            bp->resetToOrigin();

            // Set the distance between atoms
            bp->setAtomDist(offsetAtom);
            if (bp->maximumDistance() > maxDist)
            {
                maxDist = bp->maximumDistance();
            }
            templateDB_.push_back(bp);
            bp = NULL;
        }
        else
        {
            char buf[256];
            snprintf(buf, sizeof(buf), "Unknown pdbType %s", pdbType.c_str());
            GMX_THROW(InvalidInputError(buf));
        }
    }
    fp.close();
    if (NULL != bp)
    {
        delete bp;
    }
    fprintf(rnaLog, "There are %zu RNA base pair types in template database %s.\n",
            templateDB_.size(),
            dbName.c_str());
    fprintf(rnaLog, "Will use %g nm for neighborsearching\n",
            maxDist + extraOffset);

    // Initiate the neighborsearching code
    nb_.setCutoff(maxDist + extraOffset);
    // TODO: Did not work before, but we don't know why....
    nb_.setMode(AnalysisNeighborhood::eSearchMode_Grid);
}

void RnaAnalysis::readTopology(const TopologyInformation &top)
{
    // Initialize variables
    iAtoms = &(top.topology()->atoms);
    std::vector<int> rnaRes;
    const char      *grpnames;

    // TODO: residuetypes.cpp should be refactored to implement a function which
    // provides the residuetypes data, thus we don't have to do such init every time.
    // Set up the residue type
    gmx_residuetype_t *rt;
    gmx_residuetype_init(&rt);
    GMX_RELEASE_ASSERT((rt != NULL), "Problems setting up the residuetypes");

    // Get all the RNA residues
    for (int i = 0; i < iAtoms->nres; i++)
    {
        gmx_residuetype_get_type(rt, *iAtoms->resinfo[i].name, &grpnames);
        if (strcmp("RNA", grpnames) == 0)
        {
            rnaRes.push_back(i);
        }
    }

    // Set all the values of the nucleotide
    int  cAtom = 0;
    bool taken = false;

    for (size_t i = 0; i < rnaRes.size(); i++)
    {
        ResInfo *ri = NULL;

        for (; cAtom < iAtoms->nr; )
        {
            if (taken && (rnaRes[i] != iAtoms->atom[cAtom].resind))
            {
                // New nucleotide
                taken = false;
                break;
            }
            // Test if the current atom is within the right nucleotide
            if (rnaRes[i] == iAtoms->atom[cAtom].resind)
            {
                if (!taken)
                {
                    taken = true;
                    ri    = new ResInfo(*iAtoms->resinfo[rnaRes[i]].name[0],
                                        static_cast<size_t>(iAtoms->resinfo[rnaRes[i]].nr),
                                        static_cast<size_t>(cAtom),
                                        static_cast<size_t>(iAtoms->resinfo[rnaRes[i]].chainnum));
                }
                // Add an atom.
                // TODO: make sure the order of atoms in the input is the same as the template
                ri->addResInfoAtom(*iAtoms->atomname[cAtom],
                                   !isValid(*iAtoms->atomname[cAtom]));
                cAtom++;
            }
        }
        if (NULL != ri)
        {
            if (!validateResInfo(*ri))
            {
                delete ri;
                continue;
            }
            nucleotideInfo_.push_back(std::move(ri));
            if (NULL != debug)
            {
                fprintf(debug, "Just pushed a residue %c%zu with %zu atoms.\n",
                        ri->residueType(), ri->residueNumber(), ri->nAtoms());
            }
        }
    }
    fprintf(rnaLog, "There are %zu RNA bases in input file.\n", rnaRes.size());

    gmx_residuetype_destroy(rt);
}

void RnaAnalysis::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation        &top)
{
    // Open log file
    rnaLog = gmx_fio_fopen(rnaLogFile.c_str(), "w");

    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);
    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data etc.
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(rnaGraphFile);
    plotm_->setTitle("RNA Analysis");
    plotm_->setXAxisIsTime();
    plotm_->setYLabel("RNA Analysis (%)");
    plotm_->appendLegend("RNA Analysis");
    data_.addModule(plotm_);

    readTopology(top);

    readTemplates();
}

void RnaAnalysis::printTriplets(const std::list<PairInfo> &tempPairs)
{
    fprintf(rnaLog, "\n");
    bool swap1, swap2;
    for (auto p1 = tempPairs.begin(); p1 != tempPairs.end(); p1++)
    {
        for (auto p2 = p1; p2 != tempPairs.end(); p2++)
        {
            if (p1 == p2)
            {
                continue;
            }
            swap1 = p1->base1() == p2->base1() || p1->base1() == p2->base2();
            swap2 = p1->base2() == p2->base2() || p1->base1() == p2->base2();
            if (swap1 || swap2 || (p1->base2() == p2->base1()))
            {
                printTriplet(*p1, *p2, swap1, swap2);
            }
        }
    }
}

void RnaAnalysis::printTriplet(const PairInfo &pair1, const PairInfo &pair2, const bool swap1, const bool swap2)
{
    ResInfo &base1 = *nucleotideInfo_[swap1 ? pair1.base2() : pair1.base1()];
    ResInfo &base2 = *nucleotideInfo_[swap1 ? pair1.base1() : pair1.base2()];
    ResInfo &base3 = *nucleotideInfo_[swap2 ? pair2.base1() : pair2.base2()];
    fprintf(rnaLog, "triplet %c%zu - %c%zu - %c%zu",
            base1.residueTypeChar(), base1.residueNumber(),
            base2.residueTypeChar(), base2.residueNumber(),
            base3.residueTypeChar(), base3.residueNumber());
    if (DB.compare("bps") == 0)
    {
        fprintf(rnaLog, "\t\t%s - %s",
                (*pair1.bpi())->pairName(swap1 ^ pair1.swap()).c_str(),
                (*pair2.bpi())->pairName(swap2 ^ pair2.swap()).c_str());
    }
    fprintf(rnaLog, "\n");
}

void RnaAnalysis::printBasePairs(std::list<PairInfo> &tempPairs)
{
    tempPairs.remove_if([&](const PairInfo &pair)
                        {
                            bool a = false, b = false;
                            if (pair.RMSD() < 0.071)
                            {
                                return false;
                            }
                            for (const auto &pair2: tempPairs)
                            {
                                if (&pair == &pair2 || pair.score() < pair2.score())
                                {
                                    continue;
                                }
                                if (pair2.base1() == pair.base1() || pair2.base2() == pair.base1())
                                {
                                    a = true;
                                }
                                if (pair2.base1() == pair.base2() || pair2.base2() == pair.base2())
                                {
                                    b = true;
                                }
                            }
                            return a && b;
                        });
    for (const auto &pair : tempPairs)
    {
        printBasePair(pair);
    }
}

void RnaAnalysis::printBasePair(const PairInfo &tempPair)
{
    // TODO: to set the chain id correctly.
    fprintf(rnaLog, "base pair %zu%c%zu - %zu%c%zu",
            nucleotideInfo_[tempPair.base1()]->chainID(),
            nucleotideInfo_[tempPair.base1()]->residueTypeChar(),
            nucleotideInfo_[tempPair.base1()]->residueNumber(),
            nucleotideInfo_[tempPair.base2()]->chainID(),
            nucleotideInfo_[tempPair.base2()]->residueTypeChar(),
            nucleotideInfo_[tempPair.base2()]->residueNumber());

    BasePairPtrIterator bpi = tempPair.bpi();
    if (DB.compare("bps") == 0)
    {
        fprintf(rnaLog, "\t\t%s %s-%s",
                isomerismName((*bpi)->isomerismType()),
                (*bpi)->bondName(tempPair.swap()).c_str(),
                (*bpi)->subType(tempPair.swap()).c_str());
        // For debugging.
        // fprintf(rnaLog, "\t%f\t%f", tempPair.RMSD(), tempPair.score());
    }
    if (datailedInfo_)
    {
        fprintf(rnaLog, "\nRMSD to template: %f\n", tempPair.RMSD());

        const ResInfo *base1, *base2;
        if (tempPair.swap())
        {
            base2 = nucleotideInfo_[tempPair.base1()];
            base1 = nucleotideInfo_[tempPair.base2()];
        }
        else
        {
            base1 = nucleotideInfo_[tempPair.base1()];
            base2 = nucleotideInfo_[tempPair.base2()];

        }
        for (size_t con = 0; con < (*bpi)->nrBondAtoms(); con++)
        {
            size_t atomId1 = findSameAtom(*base1, (*bpi)->bondIndex(con, 0));
            size_t atomId2 = findSameAtom(*base2, (*bpi)->bondIndex(con, 1));

            atomId1 = base1->atomStart() + atomId1;
            atomId2 = base2->atomStart() + atomId2;

            // Get atom coordinates when Hbond exists and the calculate the distance
            ConstArrayRef<rvec> x           = sel_.coordinates();
            real                currentDist = std::sqrt(distance2(x[atomId1], x[atomId2]));

            if (currentDist <= bondDist_)
            {
                fprintf(rnaLog, "bond distance %s-%s: %f\n",
                        (*bpi)->bondIndex(con, 0),
                        (*bpi)->bondIndex(con, 1), currentDist);
            }
        }
    }
    fprintf(rnaLog, "\n");
}

void RnaAnalysis::print2DRepresentation()
{
    // TODO
}

void RnaAnalysis::getPlane(const Selection &sel, const ResInfo &base, rvec *out, rvec *center)
{
    size_t            n             = 3;
    const std::string ring_atoms[6] = {"N1", "N3", "C5"};
    (*center)[0] = 0;
    (*center)[1] = 0;
    (*center)[2] = 0;
    rvec v1, v2;
    for (size_t i = 0; i < n; i++)
    {
        size_t atomi = base.searchAtom(ring_atoms[i]);
        auto  &x     = sel.coordinates()[base.atomStart() + atomi];
        rvec_inc(*center, x);
        if (i == 0)
        {
            copy_rvec(x, v1);
            copy_rvec(x, v2);
        }
        else if (i == 1)
        {
            rvec_dec(v1, x);
        }
        else if (i == 2)
        {
            rvec_dec(v2, x);
        }
    }
    (*center)[0] /= n;
    (*center)[1] /= n;
    (*center)[2] /= n;
    cprod(v1, v2, *out);

    // return the normal vector.
    real nlen = norm(*out);
    (*out)[0] /= nlen;
    (*out)[1] /= nlen;
    (*out)[2] /= nlen;
}

void RnaAnalysis::getPlaneAngleDistance(const Selection &sel, const ResInfo &base1, const ResInfo &base2,
                                        real *angle, real *distance)
{
    rvec plane1, plane2;
    rvec center1, center2;
    getPlane(sel, base1, &plane1, &center1);
    getPlane(sel, base2, &plane2, &center2);
    *angle = gmx_angle(plane1, plane2);
    if (isnan(*angle))
    {
        *angle = 0;
    }
    if (*angle > M_PI_2)
    {
        *angle = static_cast<real>(M_PI) - *angle;
        for (size_t i = 0; i < DIM; i++)
        {
            plane2[i] = -plane2[i];
        }
    }
    rvec t = {0};
    rvec_add(plane1, plane2, t);
    rvec x = {0};
    x[0] = (t[0] * center2[0] + t[1] * (t[1] / t[0] * center1[0] - center1[1] + center2[1]) +
            t[2] * (t[2] / t[0] * center1[0] - center1[2] + center2[2])) /
        (t[0] + t[1] * t[1] / t[0] + t[2] * t[2] / t[0]) - center1[0];
    x[1]      = t[1] / t[0] * x[0];
    x[2]      = t[2] / t[0] * x[0];
    *distance = norm(x);
    if (isnan(*distance))
    {
        *distance = 0;
    }
}

real RnaAnalysis::analyzeBasePair(BasePairPtrIterator bpi,
                                  const ResInfo      &base1,
                                  const ResInfo      &base2,
                                  t_pbc              *pbc,
                                  const Selection    &sel,
                                  real               *RMSD,
                                  rvec               *vec,
                                  size_t              vecSize)
{
    auto bond = (*bpi)->isHBond(base1, base2,
                                pbc, bondDist_, sel);
    *RMSD = (*bpi)->computeRootMeanSquareDeviation(base1, base2, vec, vecSize);
    auto mirrored_vec = new rvec[vecSize];
    for (size_t i = 0; i < vecSize; i++)
    {
        mirrored_vec[i][XX] = vec[i][XX];
        mirrored_vec[i][YY] = vec[i][YY];
        mirrored_vec[i][ZZ] = -vec[i][ZZ];
    }
    *RMSD = std::min(*RMSD, (*bpi)->computeRootMeanSquareDeviation(base1, base2, mirrored_vec, vecSize));
    delete[] mirrored_vec;
    if (NULL != debug)
    {
        fprintf(debug, "template: %s; RMSD = %g\n",
                (*bpi)->templateName().c_str(), *RMSD);
    }
    auto score = *RMSD;
    if (score > 0.05)
    {
        score += (bond.second - bond.first) * 0.005;
        real plane_angle, plane_distance;
        getPlaneAngleDistance(sel, base1, base2, &plane_angle, &plane_distance);
        score += exp(plane_angle / M_PI_2 - 1) / 8;
        score += (exp(plane_distance) - 1) / 4;
    }
    if ((bond.first < 2 && bond.second > 0) || score > (*bpi)->getTemplateScore())
    {
        score = 1001.0;
    }
    return score;
}

//! Fill a vector with valid coordinates of a base pair
static size_t fillVector(const ResInfo      &ri1,
                         const ResInfo      &ri2,
                         ConstArrayRef<rvec> src,
                         rvec               *dest)
{
    size_t index = 0;
    for (size_t i = 0; i < ri1.nAtoms(); i++)
    {
        if (!ri1.invalid(i))
        {
            copy_rvec(src[i + ri1.atomStart()], dest[index++]);
        }
    }
    for (size_t i = 0; i < ri2.nAtoms(); i++)
    {
        if (!ri2.invalid(i))
        {
            copy_rvec(src[i + ri2.atomStart()], dest[index++]);
        }
    }
    return index;
}

void RnaAnalysis::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    std::list<PairInfo> pairs;
    rvec                vec[200];
    rvec                vec_swap[200];
    AnalysisDataHandle  dh  = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_);

    // Saves all the coordinates (only done when making statistic, stats needs to be true!)
    if (statistic)
    {
        size_t stop;
        frame  frTemp;
        frTemp.vec = new rvec[iAtoms->nr]();
        for (size_t r = 0; r < nucleotideInfo_.size(); r++)
        {
            stop = nucleotideInfo_[r]->atomStart() + nucleotideInfo_[r]->nAtoms();
            for (size_t i = nucleotideInfo_[r]->atomStart(); i < stop; i++)
            {
                copy_rvec(sel.coordinates()[i], frTemp.vec[i]);
            }
        }
        coords.push_back(frTemp);
    }

    GMX_RELEASE_ASSERT(NULL != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Set the default coordinates to search from
    rvec *repVec = new rvec[nucleotideInfo_.size()]();
    int   i      = 0;
    for (std::vector<ResInfo *>::iterator ri = nucleotideInfo_.begin();
         (ri < nucleotideInfo_.end()); ++ri, ++i)
    {
        copy_rvec(sel.coordinates()[(*ri)->atomStart() + offsetAtom], repVec[i]);
    }

    AnalysisNeighborhoodPositions pos(repVec, nucleotideInfo_.size());
    // Use neighborsearching tools
    AnalysisNeighborhoodSearch    nbsearch = nb_.initSearch(pbc, pos);
    // Find the first reference position within the cutoff.
    AnalysisNeighborhoodPair      pair;

    for (size_t base1i = 0; base1i < nucleotideInfo_.size(); base1i++)
    {
        // Set the mode to be grid mode
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(repVec[base1i]);

        // Start the search
        while (pairSearch.findNextPair(&pair))
        {
            size_t base2i = pair.refIndex();
            if (base2i <= base1i)
            {
                continue;
            }
            ResInfo  base1 = *nucleotideInfo_[base1i];
            ResInfo  base2 = *nucleotideInfo_[base2i];

#ifdef DEBUG
            if (base1.residueNumber() == 77 && base2.residueNumber() == 99)
            {
                auto t = 0;
                // For debugging.
                // This if statement is just for setting a breakpoint when debugging.
            }
#endif
            PairInfo tempPair(base1i, base2i, 1000.0, 1000.0, templateDB_.end());

            // Concatenate the two arrays and set values
            // Exclude atoms with atom types which should not be used during analysis
            // (e.g. hydrogen atoms)
            size_t dimerSize = fillVector(base1, base2,
                                          sel.coordinates(), vec);
            size_t dimerSize_Swap = fillVector(base2, base1,
                                               sel.coordinates(), vec_swap);
            GMX_RELEASE_ASSERT((dimerSize == dimerSize_Swap),
                               "Inconsistency filling vectors in two different orders");

            // Loops through the base pair template database
            for (BasePairPtrIterator bpi = templateDB_.begin();
                 (bpi < templateDB_.end()); ++bpi)
            {
                // Now check whether we have nucleotides matching the template
                BasePairComparison bpc =
                    (*bpi)->checkTemplate(base1.residueType(),
                                          base2.residueType());
                if (BasePairComparison_MisMatch == bpc)
                {
                    continue;
                }

#ifdef DEBUG
                if (    // (*bpi)->nucleotideType(0) == Uracil && (*bpi)->nucleotideType(1) == Guanine &&
                    (*bpi)->isomerismType() == Trans && (*bpi)->bondName().compare("HW") == 0)
                {
                    auto p = 2;
                    // For debugging.
                    // also for breakpoints.
                }
#endif
                real match_score = 1000.0;
                real match_rmsd  = 1000.0;
                if ((bpc == BasePairComparison_Match) ||
                    (bpc == BasePairComparison_Both))
                {
                    // Test normal base pair
                    match_score = analyzeBasePair(bpi, base1, base2, pbc,
                                                  sel, &match_rmsd,
                                                  vec, dimerSize);
                }
                real swap_score = 1000.0;
                real swap_rmsd  = 1000.0;
                if ((bpc == BasePairComparison_Swap) ||
                    (bpc == BasePairComparison_Both))
                {
                    // Test swapped base pair
                    swap_score = analyzeBasePair(bpi, base2, base1, pbc,
                                                 sel, &swap_rmsd,
                                                 vec_swap, dimerSize);
                }
                bool swap  = swap_score < match_score;
                real score = std::min(match_score, swap_score);
                if (score < tempPair.score())
                {
                    tempPair.setBPI(bpi);
                    tempPair.setScore(score);
                    tempPair.setRMSD(swap ? swap_rmsd : match_rmsd);
                    tempPair.setSwap(swap);

                    // Prints out list of base pairs every time the criteria is fullfilled
                    // Could be more than once for one single base pair
                    if (printAllBasePairs_)
                    {
                        printBasePair(tempPair);
                    }
                }
            }       // Loop over templates

            // Only true if a templated base pair is found that
            // matches RMSD and hydrogen bond distance criterion
            // for one of the templates.
            if (tempPair.bpi() < templateDB_.end())
            {
                // Save the current pair
                pairs.push_back(tempPair);

            } // end if

        }     // while search for pairs

    }         // end for each residue in PDB structure

    // Prints out base pair list once (with the template that matches the base pair the best)
    if (oneBPList)
    {
        printBasePairs(pairs);
    }
    printTriplets(pairs);

    fprintf(rnaLog, "\nFound %zu RNA base pairs that match the given criterion.\n"
            "Double, triple, etc. count might be possible for one single base pair.\n",
            pairs.size());
    found.push_back(pairs);
    nbsearch.reset();
    // Magic
    dh.finishFrame();

    delete[] repVec;
}

void
RnaAnalysis::finishAnalysis(int /* nframes*/)
{
    if (NULL != rnaLog)
    {
        gmx_fio_fclose(rnaLog);
    }
}

void
RnaAnalysis::writeOutput()
{
    if ((Xpm && (DB == "bps")) || statistic)
    {
        // Set up the default colour
        int   adjustment;
        int   atomId[200];
        t_rgb rlo;
        t_rgb rhi;
        rhi.r = 1.0, rhi.g = 0.0, rhi.b = 1.0;
        rlo.r = 0.0, rlo.g = 0.0, rlo.b = 0.0;
        std::vector<real> axis;

        // Set up axes for colour map, both axis are the same length
        for (size_t i = 0; i < nucleotideInfo_.size(); i++)
        {
            axis.push_back(nucleotideInfo_[i]->residueNumber());
        }

        // Strings
        char                titel[]  = "RNA Analysis";
        char                legend[] = "legend";
        char                label[]  = "residues";
        int                 nlevels  = 15;
        std::vector<real *> matrix;

        for (size_t x = 0; x < nucleotideInfo_.size(); x++)
        {
            matrix.push_back(new real[nucleotideInfo_.size()]());
        }

        /* Open files for storing all detected RNA base pairs */
        std::vector<FILE *> pdbStatistic;
        std::vector<size_t> model_nr;

        if (statistic)
        {
            for (std::vector<BasePair *>::iterator bpi = templateDB_.begin();
                 (bpi < templateDB_.end()); ++bpi)
            {
                char buf[256];

                // Final outpath depends on template
                if (outPath.size() > 0)
                {
                    char temp[256];
                    strcpy(temp, outPath.c_str());
                    strcat(temp, "/");

                    if (DB == "bps")
                    {
                        snprintf(buf, sizeof(buf), "%s%c%c_%s_%s.pdb",
                                 temp,
                                 (*bpi)->nucleotideTypeChar(0),
                                 (*bpi)->nucleotideTypeChar(1),
                                 (*bpi)->bondName().c_str(),
                                 isomerismName((*bpi)->isomerismType()));
                    }
                    else
                    {
                        snprintf(buf, sizeof(buf), "%s%c%c_temp.pdb",
                                 temp,
                                 (*bpi)->nucleotideTypeChar(0),
                                 (*bpi)->nucleotideTypeChar(1));
                    }
                }
                else
                {
                    if (DB == "bps")
                    {
                        snprintf(buf, sizeof(buf), "%c%c_%s_%s.pdb",
                                 (*bpi)->nucleotideTypeChar(0),
                                 (*bpi)->nucleotideTypeChar(1),
                                 (*bpi)->bondName().c_str(),
                                 isomerismName((*bpi)->isomerismType()));
                    }
                    else
                    {
                        snprintf(buf, sizeof(buf), "%c%c_temp.pdb",
                                 (*bpi)->nucleotideTypeChar(0),
                                 (*bpi)->nucleotideTypeChar(1));
                    }
                }
                pdbStatistic.push_back(gmx_fio_fopen(buf, "a+"));
                model_nr.push_back(1);
            }
        }


        // Statistic list
        std::list<frame>::const_iterator frVec = coords.begin();

        // Set up list pointers
        std::list<PairInfo>::const_iterator posI;
        std::list<PairInfo>::const_iterator posEnd;
        std::list < std::list < PairInfo>>::const_iterator frameI;
        std::list < std::list < PairInfo>>::const_iterator frameEnd = found.end();

        int ii = 0;
        // Start running through the whole trajectory
        for (frameI = found.begin(); frameI != frameEnd; ++frameI)
        {
            // Calculate the current marix
            posEnd = frameI->end();
            for (posI = frameI->begin(); posI != posEnd; ++posI)
            {
                BasePairPtrIterator bpi       = posI->bpi();
                size_t              pairIndex = bpi - templateDB_.begin();
                // Saves the positions
                if (statistic)
                {
                    adjustment = 0;
                    size_t dimerSize = (nucleotideInfo_[posI->base1()]->nAtoms() +
                                        nucleotideInfo_[posI->base2()]->nAtoms());
                    // Concatenate the two arrays and set values to them
                    for (size_t i = 0; i < dimerSize; i++)
                    {
                        // Exclude invalid atom types ( hydrogen etc)
                        if (nucleotideInfo_[posI->base1()]->nAtoms() > i)
                        {
                            if (nucleotideInfo_[posI->base1()]->invalid(i))
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the first part of the array
                            atomId[(i - adjustment)] = i + nucleotideInfo_[posI->base1()]->atomStart();
                        }
                        else
                        {
                            if (nucleotideInfo_[posI->base2()]->invalid(i - nucleotideInfo_[posI->base1()]->nAtoms()))
                            {
                                adjustment++;
                                continue;
                            }
                            // Sets the second part of the array
                            atomId[(i - adjustment)] = i + nucleotideInfo_[posI->base2()]->atomStart() -
                                nucleotideInfo_[posI->base1()]->nAtoms();
                        }
                    }
                    dimerSize = dimerSize - adjustment;

                    // Write a pdb file
                    write_pdbfile_indexed(pdbStatistic[pairIndex], " ", iAtoms,
                                          frVec->vec, (*bpi)->getEpbc(),
                                          *(*bpi)->box(),
                                          'A',
                                          model_nr[pairIndex],
                                          dimerSize, atomId, NULL, false);
                    model_nr[pairIndex]++;
                }
                if (Xpm && (DB == "bps"))
                {
                    // Get the index of the bond
                    matrix[posI->base2()][posI->base1()] = (*bpi)->bondPairIndex() + 1;
                    // Only got 16 colours
                    if (matrix[posI->base2()][posI->base1()] > 16)
                    {
                        matrix[posI->base2()][posI->base1()] = 17;
                    }
                    // Isomeri colour either 0 or 1
                    matrix[posI->base1()][posI->base2()] = static_cast<int>((*bpi)->isomerismType()) * 10;
                }

            }
            // Print the current matrix
            if (Xpm && (DB == "bps"))
            {
                FILE              *file;
                char               temp[256] = "";

                std::ostringstream ss;
                ss << ii;
                std::string        s = ss.str();
                char const        *c = s.c_str();
                // Open a file to write to
                if (outPath != "")
                {
                    strcpy(temp, outPath.c_str());
                    strcat(temp, "/");
                    strcat(temp, c);
                    strcat(temp, "_Matrix.xpm");
                }
                else
                {
                    strcat(temp, c);
                    strcat(temp, "_Matrix.xpm");
                }
                file = gmx_fio_fopen(temp, "w");
                GMX_RELEASE_ASSERT((file != NULL), "Could not create an .xpm file");

                write_xpm_split(file, 0, titel, legend,
                                label, label, nucleotideInfo_.size(), nucleotideInfo_.size(),
                                axis.data(), axis.data(), matrix.data(), 0.0, 16.0, &nlevels, rlo, rhi, 0.0, 16.0,
                                &nlevels, true, rlo, rhi);

                // Reset the matix
                for (posI = frameI->begin(); posI != posEnd; ++posI)
                {
                    matrix[posI->base2()][posI->base1()] = 0.0;
                    matrix[posI->base1()][posI->base2()] = 0.0;
                }

                ii++;
                gmx_fio_fclose(file);

            }
            if (statistic)
            {
                // We dont need the current frame anymore
                delete[] frVec->vec;
                ++frVec;
            }


        }

        // Close files and release memory
        if (statistic)
        {
            for (size_t i = 0; i < templateDB_.size(); i++)
            {
                gmx_fio_fclose(pdbStatistic[i]);
            }
        }

        for (auto &each : matrix)
        {
            delete[] each;
        }
    }
    if (Xpm && (DB != "bps"))
    {
        printf("\nCalculation of the matrix is only available if the templateDB flag is set to \"bps\"");
    }
}

// find atom with the same name
size_t RnaAnalysis::findSameAtom(const ResInfo &base, const char *name)
{
    for (size_t x = 0; (x < base.nAtoms()); x++)
    {
        if (base.atomName(x).compare(name) == 0)
        {
            return x;
        }
    }
    char buf[256];
    snprintf(buf, sizeof(buf), "Can not find atom %s", name);
    GMX_THROW(APIError(buf));
}

// Return true if its a valid atom, needs name as input
bool RnaAnalysis::isValid(const char *name)
{
    // Check for atom types we should ignore
    if ((!phosphateRmsd_ && strchr(name, 'P') != NULL) ||
        (!hydrogenRmsd_ && strchr(name, 'H') != NULL) ||
        (!sugarRmsd_ && strchr(name, '\'') != NULL))
    {
        return false;
    }

    bool set = false;
    for (size_t i = 0; i < strlen(name); i++)
    {
        set = set || isalpha(name[i]);
    }

    if (!set)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s does not seem to be an atom\n", name);
        GMX_THROW(APIError(buf));
    }
    return set;
}

real RnaAnalysis::getMass(real m, const std::string &atomname)
{
    if (atomname.find('\'') != std::string::npos)
    {
        return m * static_cast<real>(massSugar_);
    }
    else if (atomname.find('P') != std::string::npos)
    {
        return m * static_cast<real>(massPhosphate_);
    }
    return m;
}

const char RnaAnalysisInfo::name[]             = "rnaAnalysis";
const char RnaAnalysisInfo::shortDescription[] =
    "Analyze RNA PDB structures and trajectories";

TrajectoryAnalysisModulePointer RnaAnalysisInfo::create()
{
    return TrajectoryAnalysisModulePointer(new RnaAnalysis());
}

} // namespace analysismodules

} // namespace gmx
