import sys
from gromacs import TrajectoryAnalysis

class M(TrajectoryAnalysis.TrajectoryAnalysisModule):
    def __init__(self):
        super(M, self).__init__(b"a", b"a")

    def initOptions(self, options, settings):
        print('python: initOptions')
        options.setDescription(b'A stupid test module')
        settings.setFlag(TrajectoryAnalysis.TrajectoryAnalysisSettings.efRequireTop)
        print('python: inited')

    def initAnalysis(self, settings, top):
        print('python: initAnalysis')
        print('There are {} atoms'.format(top.topology().atoms.nr))
        print('Topology name: {}'.format(top.topology().name))

    def analyzeFrame(self, frnr, frame, pbc, data):
        print('python: Analyzing frame {}, {} atoms'.format(frnr, frame.natoms))
        print(frame.box[0], frame.box[1], frame.box[2])
        print(frame.x[0], frame.v[0], frame.f[0])

    def finishAnalysis(self, nframes):
        print('python: Analyzed {} frames'.format(nframes))

    def writeOutput(self):
        print('python: writeOutput')

m = M()

runner = TrajectoryAnalysis.TrajectoryAnalysisCommandLineRunner(m)
print(runner.run(sys.argv))
