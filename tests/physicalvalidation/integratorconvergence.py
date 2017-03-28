#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017 by the GROMACS development team, led by
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

# py2 compatibility
from __future__ import absolute_import, division, print_function
# other imports
import numpy as np


def calculate_rmsd(data, time=None, slope=False):
    assert isinstance(data, np.ndarray) and data.ndim == 1
    assert time is None or isinstance(time, np.ndarray) and time.ndim == 1

    avg = data.mean()

    if time is None:
        time = np.range(data.size)

    fit = np.polyfit(time, data, 1)
    f = lambda x: fit[0]*x + fit[1]

    if slope:
        rmsd = 0
        for t, d in zip(time, data):
            rmsd += (d - f(t))**2
        rmsd = np.sqrt(rmsd / data.size)
    else:
        rmsd = data.std()

    return avg, rmsd, fit[0]


def simple_convergence_test(dts, rmsds, tol):
    dt_ratio_2 = (dts[:-1] / dts[1:])**2
    rmsds = rmsds[:-1] / rmsds[1:]
    return np.allclose(rmsds, dt_ratio_2, rtol=tol, atol=0)


def check_convergence(const_traj, verbose=True, slope=False, tol=0.1):
    assert isinstance(const_traj, dict)
    assert len(const_traj) >= 2

    if verbose:
        print('{:65s}'.format('-'*65))
        print('{:>10s} {:>10s} {:>10s} {:>10s} {:^21s}'.format('dt', 'avg', 'rmsd', 'slope', 'ratio'))
        print('{:43s} {:>10s} {:>10s}'.format('', 'dt^2', 'rmsd'))
        print('{:65s}'.format('-'*65))

    prev = None

    results = {}
    for dt, traj in sorted(const_traj.items(), key=lambda x: float(x[0]), reverse=True):
        assert isinstance(traj, np.ndarray)
        assert traj.ndim == 1 or traj.ndim == 2
        dt = float(dt)

        if traj.ndim == 1:
            data = traj
            time = None
        else:
            data = traj[1]
            time = traj[0]

        results[dt] = calculate_rmsd(data, time, slope)

        if verbose:
            if prev is None:
                print('{:10.4g} {:10.2f} {:10.2e} {:10.2e} {:>10s} {:>10s}'.format(dt, results[dt][0], results[dt][1],
                                                                                   results[dt][2], '--', '--'))
                prev = [dt, results[dt][1]]
            else:
                print('{:10.4g} {:10.2f} {:10.2e} {:10.2e} {:10.2f} {:10.2f}'.format(dt, results[dt][0], results[dt][1],
                                                                                     results[dt][2],
                                                                                     prev[0]**2/dt**2,
                                                                                     prev[1]/results[dt][1]))
                prev = [dt, results[dt][1]]

    if verbose:
        print('{:65s}\n'.format('-'*65))

    dts = np.sort(np.array([float(dt) for dt in results.keys()]))[::-1]
    rmsds = np.array([float(results[dt][1]) for dt in dts])

    return simple_convergence_test(dts, rmsds, tol)
