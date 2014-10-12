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

from gromacs import TrajectoryAnalysis
from runner.pipeline import GromacsPipeline, runPipeline

class Test(TrajectoryAnalysis.TrajectoryAnalysisModule):

    def __init__(self):
        super(Test, self).__init__("test", "test")

    def initOptions(self, options, settings):
        settings.setHelpText(self.description())

    def getBatch(self):
        return self.modules

    def getArgv(self, i):
        return self.options[i]

    def initAnalysis(self, settings, top):
        pass

    def analyzeFrame(self, frnr, frame, pbc, data):
        print("Analyzing frame in Test module")

    def finishAnalysis(self, nframes):
        pass

    def writeOutput(self):
        pass

modules = [
    ("Angle", "-group1 System -oav angles.xvg"),
    (Test(), ""),
    (TrajectoryAnalysis.SasaInfo.create(), "-surface DNA"),
]

pipeline = runPipeline(name="Pipeline", modules=modules, keep_datasets=True)
dataset = pipeline.modules[0].datasetFromIndex(1)
for i in range(dataset.frameCount()):
    print('frame =', i, ', columnCount =', dataset.columnCount(), ', y =', dataset.getDataFrame(i).y(0))
