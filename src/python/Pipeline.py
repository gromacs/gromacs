#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

import sys, shlex
from gromacs.Commandline import ICommandLineOptionsModule, CommandLineParser, DummyOptionsModuleSettings
from gromacs.Options import Options, TimeUnitBehavior, FileNameOptionManager
from gromacs.TrajectoryAnalysis import TrajectoryAnalysisModule, TrajectoryAnalysisSettings, TrajectoryAnalysisRunnerCommon, SelectionOptionBehavior, SelectionCollection, AnalysisDataParallelOptions, createTrajectoryAnalysisModuleByName

class Pipeline(ICommandLineOptionsModule):

    def __init__(self):
        super(Pipeline, self).__init__()

        self.settings = TrajectoryAnalysisSettings()
        self.common = TrajectoryAnalysisRunnerCommon(self.settings)
        self.selections = SelectionCollection()

        self.modules = []

    def addModule(self, module, argv):
        if isinstance(module, str):
            module = createTrajectoryAnalysisModuleByName(module)
        elif not isinstance(module, TrajectoryAnalysisModule):
            raise TypeError("Invalid module object")

        if isinstance(argv, str):
            argv = ["pipeline_module"] + shlex.split(argv)
        else:
            argv = ["pipeline_module"] + argv

        self.modules.append({
            "module": module,
            "selections": SelectionCollection(),
            "argv": argv
        })

    def init(self, settings):
        if not self.modules:
            raise RuntimeError("No modules added - nothing to do!")

    def initOptions(self, options, settings):
        time_unit_behavior = TimeUnitBehavior()
        selection_option_behavior = SelectionOptionBehavior(self.selections, self.common.topologyProvider())

        settings.addOptionsBehavior(time_unit_behavior)
        settings.addOptionsBehavior(selection_option_behavior)

        common_options = options.addGroup()
        self.common.initOptions(common_options, time_unit_behavior)
        selection_option_behavior.initOptions(common_options)

    def initModuleOptions(self, module):
        options = Options("", "")
        manager = FileNameOptionManager()
        options.addManager(manager)

        selection_option_behavior = SelectionOptionBehavior(module["selections"], self.common.topologyProvider())
        selection_option_behavior.initBehavior(options)

        settings = TrajectoryAnalysisSettings()
        dummy_settings = DummyOptionsModuleSettings()
        settings.setOptionsModuleSettings(dummy_settings)

        selection_option_behavior.initOptions(options)

        module["module"].initOptions(options, settings)
        # TODO check if module's settings conflict with ours

        parser = CommandLineParser(options)
        parser.parse(module["argv"])

        selection_option_behavior.optionsFinishing(options)
        options.finish()

        module["module"].optionsFinished(settings)
        selection_option_behavior.optionsFinished()

    def optionsFinished(self):
        self.common.optionsFinished()

    def run(self):
        self.common.initTopology()
        topology = self.common.topologyInformation();

        for module in self.modules:
            self.initModuleOptions(module)
            module["module"].initAnalysis(self.settings, topology)

        self.common.initFirstFrame()
        self.common.initFrameIndexGroup()

        nframes = 0

        for module in self.modules:
            module["module"].initAfterFirstFrame(self.settings, self.common.frame())
            dataOptions = AnalysisDataParallelOptions()
            module["data"] = module["module"].startFrames(dataOptions, module["selections"])

        while True:
            self.common.initFrame()
            frame = self.common.frame()
            pbc = self.common.pbc()
            self.selections.evaluate(frame, pbc)

            for module in self.modules:
                module["selections"].evaluate(frame, pbc)
                module["module"].analyzeFrame(nframes, frame, pbc, module["data"])
                module["module"].finishFrameSerial(nframes);

            nframes += 1

            if not self.common.readNextFrame():
                break

        if self.common.hasTrajectory():
            print('Analyzed {} frames, last time {:.3f}.'.format(nframes, self.common.frame().time), file=sys.stderr)
        else:
            print('Analyzed topology coordinates', file=sys.stderr)

        self.selections.evaluateFinal(nframes)
        for module in self.modules:
            module["selections"].evaluateFinal(nframes)

            module["module"].finishAnalysis(nframes)
            module["module"].writeOutput()

        return 0
