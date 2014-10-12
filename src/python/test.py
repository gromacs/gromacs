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

        self.optionsHolder = Options.PyOptionsHolder()

        options.setDescription('A stupid test module')
        options.addOption(self.optionsHolder.doubleOption('py').description('option added from python just to check').required())
        options.addOption(self.optionsHolder.selectionOption('sel').description('selection option from python').required())
        options.addOption(self.optionsHolder.booleanOption('fail').description('fail :)'))
        settings.setFlag(TrajectoryAnalysis.TrajectoryAnalysisSettings.efRequireTop)
        print('python: inited')

    def initAnalysis(self, settings, top):
        print('python: initAnalysis')
        print('There are {} atoms'.format(top.topology().atoms.nr))
        print('Topology name: {}'.format(top.topology().name))

    def analyzeFrame(self, frnr, frame, pbc, data):
        sel = self.optionsHolder['sel']
        print('selected atoms {}, {}'.format(sel.atomCount(), sel.coordinates()[0]))
        print('ids: {}'.format(sel.mappedIds()[:5]))
        return
        print('python: Analyzing frame {}, {} atoms'.format(frnr, frame.natoms))
        print(frame.box[0], frame.box[1], frame.box[2])
        print(frame.x[0], frame.v[0], frame.f[0])

    def finishAnalysis(self, nframes):
        print('python: Analyzed {} frames'.format(nframes))

    def writeOutput(self):
        print('python: writeOutput')
        print('py = {}'.format(self.optionsHolder['py']))
        print('sel = {}'.format(self.optionsHolder['sel']))

TrajectoryAnalysis.runAsMain(M(), sys.argv)
