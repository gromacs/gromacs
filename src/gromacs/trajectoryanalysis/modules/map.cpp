/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led
 * by
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
 * Implements gmx::analysismodules::map.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_trajectoryanalysis
 */

// #include "gmxpre.h"
//
// #include "map.h"
//
// #include "gromacs/analysisdata/analysisdata.h"
// #include "gromacs/externalpotential/modules/densityfitting/densityspreader.h"
// #include "gromacs/externalpotential/modules/densityfitting/emscatteringfactors.h"
// #include "gromacs/fileio/pdbio.h"
// #include "gromacs/fileio/griddataio.h"
// #include "gromacs/math/do_fit.h"
// #include "gromacs/math/vec.h"
// #include "gromacs/math/quaternion.h"
// #include "gromacs/math/griddata/rotatedgrid.h"
// #include "gromacs/math/griddata/encompassinggrid.h"
// #include "gromacs/math/griddata/operations/modifygriddata.h"
// #include "gromacs/math/griddata/operations/gridinterpolator.h"
// #include "gromacs/math/griddata/operations/realfieldmeasure.h"
// #include "gromacs/options/basicoptions.h"
// #include "gromacs/options/filenameoption.h"
// #include "gromacs/options/ioptionscontainer.h"
// #include "gromacs/topology/atomprop.h"
// #include "gromacs/topology/atoms.h"
// #include "gromacs/topology/topology.h"
// #include "gromacs/trajectory/trajectoryframe.h"
// #include "gromacs/trajectoryanalysis/analysissettings.h"
// #include "gromacs/utility/exceptions.h"
//
//
// namespace gmx
// {
//
// namespace analysismodules
// {
//
// namespace
// {
//
// class Map : public TrajectoryAnalysisModule
// {
//     public:
//         Map()  = default;
//         ~Map() = default;
//
//         virtual void initOptions(IOptionsContainer          *options,
//                                  TrajectoryAnalysisSettings *settings);
//         virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
//                                   const TopologyInformation        &top);
//         virtual void optionsFinished(TrajectoryAnalysisSettings *settings);
//
//         virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
//                                   TrajectoryAnalysisModuleData *pdata);
//
//         virtual void finishAnalysis(int nframes);
//         virtual void writeOutput();
//
//     private:
//         void set_finitegrid_from_box(matrix box, rvec translation);
//         void set_box_from_frame(const t_trxframe &fr, matrix box, rvec translation);
//         void frameToDensity_(const t_trxframe &fr, int nFr);
//
//         std::string                            fnmapinput_;
//         std::string                            fnmapoutput_;
//
//         float                                  sigma_   = 0.4;
//         float                                  n_sigma_ = 5;
//         GridDataReal3D                         inputdensity_;
//         GridDataReal3D                         outputDensityBuffer_;
//         real                                   spacing_       = 0.2;
//         bool                                   bRigidBodyFit_ = true;
//         std::vector<float>                     weight_;
//         int                                    every_                    = 1;
//         real                                   expAverage_               = 1.;
//         int                                    nFr_                      = 0;
//         bool                                   bFitFramesToTopStructure_ = false;
//         int                                    nAtomsReference_;
//         std::vector<RVec>                      referenceX_;
//         matrix                                 referenceBox;
//         std::vector<real>                      fitWeights;
//         bool                                   bUseBox_ = false;
//         std::unique_ptr<DensitySpreader>       spreader_;
//         Quaternion::QVec                       orientation_ = {{1, 0, 0, 1}};
// };
//
// void Map::initOptions(IOptionsContainer          *options,
//                       TrajectoryAnalysisSettings *settings)
// {
//
//     static const char *const desc[] = {
//         "[THISMODULE] is a tool to read in and write out (electron) density "
//         "maps.",
//         "With this tool you can produce a density map from an input ",
//         "coordinate file to be embedded in a .tpr file with grompp for",
//         "example.[PAR]", "Possible uses are:[PAR]",
//         "* Provide an input structure with [TT]-f[tt] and output a density map "
//         "based on",
//         "the grid settings. Note that you can also input a whole trajectory, in "
//         "that",
//         "case, a series of map files will be output. Use and index file to "
//         "spread just",
//         "a part of the atoms[BR]", "* Provide a map with [TT]-mi[tt] and and "
//         "output some map characteristics[BR]",
//         "* Provide a [TT].tpr[tt] file with [TT]-s[tt] and output the embedded "
//         "map with [TT]-mo[tt]"
//     };
//
//     settings->setHelpText(desc);
//     settings->setFlag(TrajectoryAnalysisSettings::efUseTopX, true);
//
//     options->addOption(FileNameOption("mi")
//                            .filetype(eftgriddata)
//                            .inputFile()
//                            .store(&fnmapinput_)
//                            .defaultBasename("ccp4in")
//                            .description("CCP4 density map input file"));
//     options->addOption(FileNameOption("mo")
//                            .filetype(eftgriddata)
//                            .outputFile()
//                            .store(&fnmapoutput_)
//                            .defaultBasename("ccp4out")
//                            .description("CCP4 density map output file"));
//     options->addOption(FloatOption("sigma").store(&sigma_).description(
//                                "Create a simulated density by replacing the atoms by Gaussian functions "
//                                "of width sigma (nm)"));
//     options->addOption(FloatOption("N_sigma").store(&n_sigma_).description(
//                                "How many Gaussian width shall be used for spreading?"));
//     options->addOption(FloatOption("spacing").store(&spacing_).description(
//                                "Spacing of the density grid (nm)"));
//     options->addOption(IntegerOption("every").store(&every_).description(
//                                "Analyse only -every frame."));
//     options->addOption(BooleanOption("rigidBodyFit")
//                            .store(&bRigidBodyFit_)
//                            .description("Use rigid body fitting in all steps."));
//     options->addOption(FloatOption("orientation").store(orientation_.data()).vector().valueCount(4));
//     options->addOption(BooleanOption("useBox").store(&bUseBox_).description(
//                                "Use the box information in the structure file to setup the grid."));
//     options->addOption(
//             FloatOption("expaverage")
//                 .store(&expAverage_)
//                 .description(
//                     "Factor for exponential averaging new_map = f*curr+(1-f)*old."));
// }
//
// void Map::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/,
//                        const TopologyInformation       &top)
// {
//     if (top.topology())
//     {
//         auto atomprop = gmx_atomprop_init();
//         get_pdb_atomnumber(&(top.topology()->atoms), atomprop);
//         for (int i_atom = 0; i_atom < top.topology()->atoms.nr; i_atom++)
//         {
//             weight_.push_back(atomicNumber2EmScatteringFactor(
//                                       top.topology()->atoms.atom[i_atom].atomnumber));
//         }
//         gmx_atomprop_destroy(atomprop);
//     }
//
//     if (bFitFramesToTopStructure_)
//     {
//         if (top.topology())
//         {
//
//             nAtomsReference_ = top.topology()->atoms.nr;
//             fitWeights.resize(nAtomsReference_, 1);
//
//             rvec *x = as_rvec_array(referenceX_.data());
//             top.getTopologyConf(&x, referenceBox);
//
//             reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, x,
//                     fitWeights.data());
//             referenceX_.assign(x, x + nAtomsReference_);
//         }
//         else
//         {
//             GMX_THROW(InvalidInputError("Please provide a topology to fit to."));
//         }
//     }
// }
//
// void Map::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
// {
//
//     if (!fnmapinput_.empty())
//     {
//         MrcFile ccp4inputfile;
//         inputdensity_ = ccp4inputfile.read(fnmapinput_);
//     }
// }
//
// void Map::frameToDensity_(const t_trxframe &fr, int nFr)
// {
//     if (nFr == 1)
//     {
//         std::unique_ptr < IGrid < DIM>> outputGrid;
//         if (fnmapinput_.empty())
//         {
//             // Guess the extend of the map from the structure
//             outputGrid = encompassingGridFromCoordinates({fr.x, fr.x+fr.natoms}, spacing_, n_sigma_*sigma_).duplicate();
//         }
//         else
//         {
//             // copy the grid properties from reference grid
//             outputGrid  = inputdensity_.getGrid().duplicate();
//         }
//
//         spreader_ = std::unique_ptr<DensitySpreader>(
//                     new DensitySpreader(*outputGrid, 1, n_sigma_, sigma_));
//
//         outputDensityBuffer_ = GridDataReal3D(*outputGrid);
//         ModifyGridData(outputDensityBuffer_).zero();
//     }
//
//     std::vector<RVec> coordinates(fr.x, fr.x + fr.natoms);
//
//     const auto       &spreaddensity = spreader_->spreadLocalAtoms(coordinates, weight_);
//
//     if (nFr == 1)
//     {
//         std::copy(std::begin(spreaddensity), std::end(spreaddensity), std::begin(outputDensityBuffer_));
//     }
//     else
//     {
//
//         // auto linearAveraging = [this](const real & current, const real &
//         // average){return (average * (nFr_ - 1) + current) / nFr_;};
//         auto exponentialAveraging = [this](const real &current,
//                                            const real &average) {
//                 return expAverage_ * current + (1 - expAverage_) * average;
//             };
//         std::transform(std::begin(spreaddensity), std::end(spreaddensity), std::begin(spreaddensity), std::begin(outputDensityBuffer_), exponentialAveraging);
//     }
// }
//
// void Map::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/,
//                        TrajectoryAnalysisModuleData * /*pdata*/)
// {
//     if (frnr % every_ == 0)
//     {
//         ++nFr_;
//         if (bFitFramesToTopStructure_)
//         {
//             reset_x(nAtomsReference_, nullptr, nAtomsReference_, nullptr, fr.x,
//                     fitWeights.data());
//             do_fit(nAtomsReference_, fitWeights.data(),
//                    as_rvec_array(referenceX_.data()), fr.x);
//         }
//
//         if (nFr_ == 1)
//         {
//             if (weight_.size() < std::size_t(fr.natoms))
//             {
//                 weight_.assign(fr.natoms, 1);
//                 fprintf(stderr, "\nSetting atom weights to unity.\n");
//             }
//         }
//
//         frameToDensity_(fr, nFr_);
//
//         if (!fnmapoutput_.empty())
//         {
//             if (frnr == 0)
//             {
//                 MrcFile().write( fnmapoutput_, outputDensityBuffer_);
//             }
//             MrcFile().write( fnmapoutput_.substr(0, fnmapoutput_.size() - 5) + std::to_string(frnr) + ".ccp4", outputDensityBuffer_);
//         }
//
//     }
// }
//
// void Map::finishAnalysis(int /*nframes*/) {}
//
// void Map::writeOutput() {}
//
// }       // namespace
//
// const char MapInfo::name[]             = "map";
// const char MapInfo::shortDescription[] = "Spread atoms on grid";
//
// TrajectoryAnalysisModulePointer MapInfo::create()
// {
//     return TrajectoryAnalysisModulePointer(new Map);
// }
//
// } // namespace analysismodules
//
// } // namespace gmx
