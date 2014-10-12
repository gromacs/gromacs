/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

// These classes belong to gromacs.Selection, but will remain here until
// t_trxframe and t_pbc are split from TrajectoryAnalysis
class SelectionCollection /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/selection/selectioncollection.h>
%End
public:
    SelectionCollection();
    void evaluate(t_trxframe *fr, t_pbc *pbc);
    void evaluateFinal(int nframes);
};

class ITopologyProvider;
class SelectionOptionBehavior : public IOptionsBehavior /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/selection/selectionoptionbehavior.h>
using gmx::SelectionOptionBehavior;
using gmx::ITopologyProvider;
%End
public:
    SelectionOptionBehavior(SelectionCollection *selections,
                            ITopologyProvider   *topologyProvider);
    void initOptions(IOptionsContainer *options);
    virtual void initBehavior(Options *options);
    virtual void optionsFinishing(Options *options);
    virtual void optionsFinished();
};

class TrajectoryAnalysisRunnerCommon /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/trajectoryanalysis/runnercommon.h>
%End
public:
    explicit TrajectoryAnalysisRunnerCommon(TrajectoryAnalysisSettings *settings);
    ~TrajectoryAnalysisRunnerCommon();

    ITopologyProvider *topologyProvider();

    void initOptions(IOptionsContainer *options, TimeUnitBehavior *timeUnitBehavior);
    void optionsFinished();
    void initTopology();
    void initFirstFrame();
    void initFrameIndexGroup();
    bool readNextFrame();
    void initFrame();

    bool hasTrajectory() const;
    const TopologyInformation &topologyInformation() const;
    t_trxframe                &frame() const;
    t_pbc                     *pbc() const;
};
