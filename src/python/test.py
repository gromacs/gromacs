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
from gromacs import Options, TrajectoryAnalysis

class M(TrajectoryAnalysis.TrajectoryAnalysisModule):
    def __init__(self):
        super(M, self).__init__("a", "a")

    def initOptions(self, options, settings):
        print('python: initOptions')

        settings.setHelpText('A stupid test module')

        self.optionsHolder = Options.PyOptionsHolder()

        options.addOption(self.optionsHolder.selectionOption('sel').required())
        options.addOption(self.optionsHolder.fileNameOption('file').defaultBasename('test').description('filename from python to rule them all').outputFile().required().filetype(Options.eftGenericData))
        settings.setFlag(TrajectoryAnalysis.TrajectoryAnalysisSettings.efRequireTop)

        self.angle = TrajectoryAnalysis.AngleInfo.create()

        print('python: inited')

    def getBatch(self):
        print('python: getBatch')
        return [self.angle]

    def getArgv(self, i):
        print('python: getArgv')
        if i == 0:
            #First element of list should be module name -- gets discarded by parser anyway
            return ["gangle", "-group1", "Backbone", "-oav", "angles.xvg"]

    def initAnalysis(self, settings, top):
        print('python: initAnalysis')
        print('There are {} atoms'.format(top.topology().atoms.nr))
        print('Topology name: {}'.format(top.topology().name))

        #Tell GROMACS to keep last frame in storage, required in analyzeFrame()
        self.angle.datasetFromIndex(1).requestStorage(1)

    def analyzeFrame(self, frnr, frame, pbc, data):
        sel = self.optionsHolder['sel']

        dataset = self.angle.datasetFromIndex(1)

        print('frame =', frnr, ', columnCount =', dataset.columnCount(), ', y =', dataset.getDataFrame(frnr).y(0))

    def finishAnalysis(self, nframes):
        print('python: Analyzed {} frames'.format(nframes))

    def writeOutput(self):
        print('python: writeOutput')
        print('file = {}'.format(self.optionsHolder['file']))

        fo = open(self.optionsHolder['file'], 'w')
        fo.write('Test output\n')

TrajectoryAnalysis.runAsMain(M(), sys.argv)
