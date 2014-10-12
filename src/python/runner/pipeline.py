#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

import sys
import shlex
from gromacs import Options
from gromacs import TrajectoryAnalysis

class GromacsPipeline(TrajectoryAnalysis.TrajectoryAnalysisModule):

    def __init__(self,
        name="PythonPipeline",
        description="Pipeline of modules created with Python",
        modules=[],
        keep_datasets=False,
    ):
        super(GromacsPipeline, self).__init__(name, description)

        self.optionsHolder = Options.PyOptionsHolder()

        self.modules = []
        self.options = []

        for module, options in modules:
            if not isinstance(module, TrajectoryAnalysis.TrajectoryAnalysisModule):
                info_name = module + "Info"
                if not hasattr(TrajectoryAnalysis, info_name):
                    raise ValueError("There is no module named {}".format(name))

                module = getattr(TrajectoryAnalysis, info_name).create()

            options_list = [name] + shlex.split(options)

            if keep_datasets:
                for i in range(module.datasetCount()):
                    module.datasetFromIndex(i).requestStorage(-1)

            self.modules.append(module)
            self.options.append(options_list)

    def initOptions(self, options, settings):
        settings.setHelpText(self.description())

    def getBatch(self):
        return self.modules

    def getArgv(self, i):
        return self.options[i]

    def initAnalysis(self, settings, top):
        pass

    def analyzeFrame(self, frnr, frame, pbc, data):
        pass

    def finishAnalysis(self, nframes):
        pass

    def writeOutput(self):
        pass

    def run(self, argv):
        if argv is None:
            argv = sys.argv

        TrajectoryAnalysis.runAsMain(self, argv)

def runPipeline(argv=None, *args, **kwargs):
    pipeline = GromacsPipeline(*args, **kwargs)
    pipeline.run(argv)

    return pipeline
