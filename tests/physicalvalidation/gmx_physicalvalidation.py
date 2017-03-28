#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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
import sys
import os
import argparse
import pygmx
import numpy as np
import integratorconvergence as intconv
import shutil
import errno


def test_integratorconvergence(cwd, gmx, it, source_path=None, verb=0, tol=0.1):

    if source_path is not None:
        try:
            shutil.rmtree(os.path.join(cwd, 'integratorconvergence'))
        except OSError as e:
            if e.errno == errno.ENOENT:
                pass
            else:
                raise
        shutil.copytree(os.path.join(source_path, 'integratorconvergence'),
                        os.path.join(cwd, 'integratorconvergence'))

    test_systems = list()

    # system 1
    # 900 H2O molecules
    test_systems.append(pygmx.System.from_deffnm('water', os.path.join(cwd, 'integratorconvergence/system_1')))

    ret_val = 0
    for system in test_systems:

        # read mdp. Different parameter settings to be tested are separated by '%'
        parameters, options = pygmx.read_mdp(system.mdp, sep='%')
        system.parameters = parameters

        # We need a unique starting time step
        if 'dt' in options or 'dt' not in parameters:
            print('ERROR: Cannot determine dt.')
            return 1
        dt = float(parameters['dt'])

        # Check if we have multiple parameter combinations
        if not options:
            systems = [system]
        else:
            systems = pygmx.split_systems(system, parameters, options)

        # loop over parameter combinations
        for s in systems:
            # set conserved quantity
            if 'tcoupl' in s.parameters and s.parameters['tcoupl'] != 'no':
                conserved = 'Conserved-En.'
            else:
                conserved = 'Total-Energy'

            s.quantities = [conserved]

            if dt > 0.001 and ('constraints' not in s.parameters or s.parameters['constraints'] == 'none'):
                print('WARNING: Initial timestep > 1 fs, but no constraints.')

            dts = np.array([dt/2**n for n in range(0, it+1)])
            str_dts = dts.astype(str).tolist()

            s.parameters['dt'] = str_dts
            subsystems = pygmx.split_systems(s, s.parameters, ['dt'], dt=dt)

            results = {}
            for n, ss in enumerate(subsystems):
                ss.run_gromacs(gmx, verb=verb)
                results[ss.parameters['dt']] = np.array([ss.values['time'][2**n-1::2**n],
                                                         ss.values[conserved][2**n-1::2**n]])

            if verb > 0:
                verbose = True

            if verbose:
                if options:
                    print('System \'{:s}\' with parameters'.format(s.deffnm))
                    for key in options:
                        print('  {:25s} = {:s}\n'.format(key, s.parameters[key]))
                else:
                    print('System \'{:s}\''.format(s.deffnm))

            ret_val += not intconv.check_convergence(results, verbose=verbose, tol=tol)
            # intconv.check_convergence(results, verbose=verbose, slope=True)

    return ret_val


def test_pass():
    print('SUCCESS')
    return 0


def test_fail():
    print('FAIL')
    return 1


def main(args):
    parser = argparse.ArgumentParser(
        description='Physical validation suite for GROMACS.',
        prog='gmx_physicalvalidation.py')
    parser.add_argument('testname', type=str,
                        metavar='testname',
                        help='Name of test to be ran.')
    parser.add_argument('--gmx', type=str,
                        metavar='exe', default=None,
                        help='Gromacs executable (default: gmx).')
    parser.add_argument('--double', action='store_true',
                        default=False,
                        help='Use double precision (default: no).')
    parser.add_argument('--wd', type=str,
                        metavar='dir', default=None,
                        help='Working directory (default: current directory)')
    parser.add_argument('--suffix', type=str,
                        metavar='_xxx', default=None,
                        help='Suffix to be used with gmx executable (default: _d if --double, none otherwise)')

    args = parser.parse_args(args)

    source_path = os.path.dirname(os.path.realpath(__file__))
    
    if args.wd is None:
        args.wd = os.getcwd()
    else:
        try:
            os.makedirs(args.wd)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(args.wd):
                pass
            else:
                raise

    if args.double and args.suffix is None:
        args.suffix = '_d'

    if args.gmx is None:
        args.gmx = 'gmx'
        if args.suffix is not None:
            args.gmx += args.suffix

    gmx = pygmx.PyGMX(exe=args.gmx, dp=args.double)

    if args.testname == 'pass':
        exit(test_pass())
    if args.testname == 'fail':
        exit(test_fail())
    if args.testname == 'integratorconvergence':
        if args.double:
            exit(test_integratorconvergence(cwd=args.wd, gmx=gmx, it=4, source_path=source_path, verb=3, tol=0.1))
        else:
            exit(test_integratorconvergence(cwd=args.wd, gmx=gmx, it=1, source_path=source_path, verb=3, tol=0.25))

if __name__ == "__main__":
    main(sys.argv[1:])
