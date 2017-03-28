from __future__ import print_function, division, absolute_import

import sys
import os
import shutil
import json
import argparse
import re
from collections import OrderedDict

from physical_validation import integrator, ensemble, kinetic_energy
from physical_validation.util.gromacs_interface import GromacsInterface
from physical_validation.data.gromacs_parser import GromacsParser


def mkdir_bk(dirname, verbose=False, nobackup=False):
    if os.path.exists(dirname) and nobackup:
        shutil.rmtree(dirname)
    elif os.path.exists(dirname):
        if verbose:
            print('Directory ' + dirname + ' exists. Backing it up.')
        basename = os.path.basename(dirname)
        n = 1
        bk_dir = '#' + basename + '_' + str(n) + '#'
        while os.path.exists(dirname.replace(basename, bk_dir)):
            n += 1
            bk_dir = '#' + basename + '_' + str(n) + '#'
        os.rename(dirname, dirname.replace(basename, bk_dir))
    os.makedirs(dirname)


def file_bk(filename, verbose=False):
    if os.path.exists(filename):
        if verbose:
            print('File ' + filename + ' exists. Backing it up.')
        basename = os.path.basename(filename)
        n = 1
        bk_file = '#' + basename + '_' + str(n) + '#'
        while os.path.exists(filename.replace(basename, bk_file)):
            n += 1
            bk_file = '#' + basename + '_' + str(n) + '#'
        os.rename(filename, filename.replace(basename, bk_file))


def basic_run_cmds(directory, grompp_args=None, mdrun_args=None):
    grompp = '$GROMPPCMD -f system.mdp -p system.top -c system.gro -o system.tpr'
    if grompp_args:
        for arg in grompp_args:
            grompp += ' ' + arg
    mdrun = '$MDRUNCMD -s system.tpr -deffnm system'
    if mdrun_args:
        for arg in mdrun_args:
            mdrun += ' ' + arg
    return [
        'oldpath=$PWD',
        'cd ' + directory,
        grompp,
        mdrun,
        'cd $oldpath'
    ]


class Test(object):
    @classmethod
    def parser(cls):
        raise NotImplementedError

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        raise NotImplementedError

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        raise NotImplementedError

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name):
        raise NotImplementedError

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity):
        raise NotImplementedError


class IntegratorTest(Test):
    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(
            description='Tests the integrator convergence.',
            formatter_class=argparse.RawTextHelpFormatter,
            prog=cls.__name__
        )
        parser.add_argument('-n', '--n_iterations', type=int, default=2,
                            help='The number of different time steps tested.')
        parser.add_argument('-t', '--tolerance', type=float, default=0.1,
                            help=('The relative tolerance accepted. Default: 0.1.\n'
                                  'Example: By default, the test passes if\n'
                                  '         3.6 <= dE(2*dt)/dE(dt) <= 4.4.')
                            )

        return parser

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        args = cls.parser().parse_args(args)
        return cls.prepare(input_dir, target_dir, system_name, nobackup,
                           n_iterations=args.n_iterations)

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        args = cls.parser().parse_args(args)
        return cls.analyze(gmx_parser, system_dir, system_name, base_data, verbosity,
                           tolerance=args.tolerance, n_iterations=args.n_iterations)

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name, nobackup, n_iterations=None):
        # Standard value
        if n_iterations is None:
            n_iterations = cls.parser().get_default('n_iterations')
        # Read base options
        options = GromacsInterface.read_mdp(os.path.join(input_dir, 'system.mdp'))

        # Check if options relevant to integrator tests are set
        # Set to standard if not - I am not happy about hardcoding this here...
        # For future: Find a way to obtain standards from GROMACS build
        if 'nsteps' not in options:
            raise ValueError('System ' + system_name + ' has no \'nsteps\' defined in mdp file. ' +
                             'Running integrator test does not make sense.')
        if 'dt' not in options:
            options['dt'] = str(0.001)
        if 'nstcalcenergy' not in options:
            options['nstcalcenergy'] = str(100)
        if 'nstenergy' not in options:
            options['nstenergy'] = str(1000)
        # if 'nstlist' not in options:
        #     options['nstlist'] = str(10)

        # Prepare folders for iterations
        directories = []
        for n in range(1, n_iterations+1):
            current_dir = os.path.join(target_dir, 'integrator_' + str(n))
            mkdir_bk(current_dir, nobackup=nobackup)
            # update timesteps, length and intervals
            options['dt'] = str(float(options['dt'])*0.5)
            for key in ['nsteps', 'nstcalcenergy', 'nstenergy']:  # , 'nstlist']:
                options[key] = str(int(options[key])*2)
            # write mdp
            GromacsInterface.write_mdp(options, os.path.join(current_dir, 'system.mdp'))
            # copy topology and starting structure
            for suffix in ['gro', 'top']:
                shutil.copy2(os.path.join(input_dir, 'system.' + suffix),
                             current_dir)
            directories.append(current_dir)

        return directories

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity,
                tolerance=None, n_iterations=None):
        # Standard value
        if tolerance is None:
            tolerance = cls.parser().get_default('tolerance')
        if n_iterations is None:
            n_iterations = cls.parser().get_default('n_iterations')

        # list of results from simulations at different time steps
        results = []

        # base data
        if base_data['reduced'] is None:
            current_dir = os.path.join(system_dir, 'base')
            base_data['reduced'] = gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.gro')
            )
        results.append(base_data['reduced'])

        # append data at reduced dt
        for n in range(1, n_iterations+1):
            current_dir = os.path.join(system_dir, 'integrator_' + str(n))
            results.append(gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.gro')
            ))

        # run test
        deviation = integrator.convergence(simulations=results,
                                           verbose=(verbosity > 0),
                                           convergence_test='max_deviation')

        result = deviation <= tolerance
        if result:
            message = 'IntegratorTest PASSED (max deviation = {:.2g}, tolerance = {:.2g})'.format(
                deviation, tolerance
            )
        else:
            message = 'IntegratorTest FAILED (max deviation = {:.2g}, tolerance = {:.2g})'.format(
                deviation, tolerance
            )

        return {'test': result,
                'result': deviation,
                'tolerance': tolerance,
                'message': message}


class EnsembleTest(Test):
    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(
            description='Tests the validity of the potential energy and / or volume ensemble.',
            formatter_class=argparse.RawTextHelpFormatter,
            prog=cls.__name__
        )
        parser.add_argument('--dtemp', nargs='*', type=float, default=None,
                            help='Ensemble validations are made between two simulations at\n'
                                 'different state points.\n'
                                 'dtemp determines the temperature difference between base\n'
                                 'simulation and the additional point. If more than one\n'
                                 'value is given, several tests will be performed.\n'
                                 'By also giving dpress, both temperature and pressure can\n'
                                 'be displaced simultaneously.')
        parser.add_argument('--dpress', nargs='*', type=float, default=None,
                            help='Ensemble validations are made between two simulations at\n'
                                 'different state points.\n'
                                 'dpress determines the pressure difference between base\n'
                                 'simulation and the additional point. If more than one\n'
                                 'value is given, several tests will be performed.\n'
                                 'By also giving dtemp, both temperature and pressure can\n'
                                 'be displaced simultaneously.')
        parser.add_argument('-t', '--tolerance', type=float, default=3,
                            help=('The number of standard deviations a result can be off\n'
                                  'to be still accepted. Default: 3.'))

        return parser

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        args = cls.parser().parse_args(args)
        return cls.prepare(input_dir, target_dir, system_name, nobackup,
                           dtemp=args.dtemp, dpress=args.dpress)

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        args = cls.parser().parse_args(args)
        return cls.analyze(gmx_parser, system_dir, system_name, base_data, verbosity,
                           tolerance=args.tolerance, dtemp=args.dtemp, dpress=args.dpress)

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name, nobackup, dtemp=None, dpress=None):
        # No standard values (system-dependent!)
        if not dtemp and not dpress:
            raise ValueError('Ensemble test for system ' + system_name +
                             ' has no defined temperature or pressure difference.')
        # Pad arrays - assume 0 difference if other difference is set
        if not dtemp:
            dtemp = [0] * len(dpress)
        if not dpress:
            dpress = [0] * len(dtemp)
        if len(dtemp) < len(dpress):
            dtemp.extend([0] * (len(dpress) - len(dtemp)))
        if len(dpress) < len(dtemp):
            dpress.extend([0] * (len(dtemp) - len(dpress)))
        # Check if we do any pressure differences
        no_press = True
        for dp in dpress:
            if dp*dp >= 1e-12:
                no_press = False
                break

        # Read base options
        options = GromacsInterface.read_mdp(os.path.join(input_dir, 'system.mdp'))

        # Check for target temperature
        if 'ref-t' in options:
            ref_t = float(options['ref-t'])
            ref_t_key = 'ref-t'
        elif 'ref_t' in options:
            ref_t = float(options['ref_t'])
            ref_t_key = 'ref_t'
        else:
            raise ValueError('Ensemble test for system ' + system_name + ' selected, ' +
                             'but system has no defined temperature.')

        # Check for target pressure, if relevant
        if no_press:
            ref_p = None
            ref_p_key = None
        elif 'ref-p' in options:
            ref_p = float(options['ref-p'])
            ref_p_key = 'ref-p'
        elif 'ref_p' in options:
            ref_p = float(options['ref_p'])
            ref_p_key = 'ref_p'
        else:
            raise ValueError('Ensemble test with pressure difference for system ' + system_name +
                             ' selected, but system has no defined pressure.')

        # Loop over tests
        directories = []
        for n, (dt, dp) in enumerate(zip(dtemp, dpress)):
            current_dir = os.path.join(target_dir, 'ensemble_' + str(n+1))
            mkdir_bk(current_dir, nobackup=nobackup)
            # change temperature & pressure
            options[ref_t_key] = str(ref_t + dt)
            if not no_press:
                options[ref_p_key] = str(ref_p + dp)
            # write mdp
            GromacsInterface.write_mdp(options, os.path.join(current_dir, 'system.mdp'))
            # copy topology and starting structure
            for suffix in ['gro', 'top']:
                shutil.copy2(os.path.join(input_dir, 'system.' + suffix),
                             current_dir)
            directories.append(current_dir)

        return directories

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity,
                tolerance=None, dtemp=None, dpress=None):
        # No standard values (system-dependent!)
        if not dtemp and not dpress:
            raise ValueError('Ensemble test for system ' + system_name +
                             ' has no defined temperature or pressure difference.')
        # Pad arrays - assume 0 difference if other difference is set
        if not dtemp:
            dtemp = [0] * len(dpress)
        if not dpress:
            dpress = [0] * len(dtemp)
        if len(dtemp) < len(dpress):
            dtemp.extend([0] * (len(dpress) - len(dtemp)))
        if len(dpress) < len(dtemp):
            dpress.extend([0] * (len(dtemp) - len(dpress)))

        nsystems = len(dtemp)

        if tolerance is None:
            tolerance = cls.parser().get_default('tolerance')

        # base data
        if base_data['reduced'] is None:
            current_dir = os.path.join(system_dir, 'base')
            base_data['reduced'] = gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.gro')
            )
        base_result = base_data['reduced']

        # list of results from simulations at different state points
        results = []
        for n in range(nsystems):
            current_dir = os.path.join(system_dir, 'ensemble_' + str(n+1))
            results.append(gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.gro')
            ))

        # run tests
        passed = True
        message = ''
        max_quantiles = -1
        for result, dt, dp in zip(results, dtemp, dpress):
            quantiles = ensemble.check(base_result, result, quiet=(verbosity == 0))
            # filename=os.path.join(system_dir, system_name + '_ens'))
            if any(q > tolerance for q in quantiles):
                passed = False
                if len(quantiles) == 1:
                    message += '\n    --dtemp={:.1f} --dpress={:.1f} : FAILED ({:.1f} quantiles off)'.format(
                        dt, dp, quantiles[0]
                    )
                else:
                    message += '\n    --dtemp={:.1f} --dpress={:.1f} : FAILED ([{:.1f}, {:.1f}] quantiles off)'.format(
                        dt, dp, quantiles[0], quantiles[1]
                    )
            else:
                if len(quantiles) == 1:
                    message += '\n    --dtemp={:.1f} --dpress={:.1f} : PASSED ({:.1f} quantiles off)'.format(
                        dt, dp, quantiles[0]
                    )
                else:
                    message += '\n    --dtemp={:.1f} --dpress={:.1f} : PASSED ([{:.1f}, {:.1f}] quantiles off)'.format(
                        dt, dp, quantiles[0], quantiles[1]
                    )
            max_quantiles = max(max_quantiles, max(quantiles))

        if passed:
            message = ('EnsembleTest PASSED (tolerance: {:.1f} quantiles)'.format(tolerance) +
                       message)
        else:
            message = ('EnsembleTest FAILED (tolerance: {:.1f} quantiles)'.format(tolerance) +
                       message)

        return {'test': passed,
                'result': max_quantiles,
                'tolerance': tolerance,
                'message': message}


class MaxwellBoltzmannTest(Test):
    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(
            description='Tests the validity of the kinetic energy ensemble.',
            formatter_class=argparse.RawTextHelpFormatter,
            prog=cls.__name__
        )
        parser.add_argument('-t', '--tolerance', type=float, default=0.05,
                            help='The alpha value (confidence interval) used in\n'
                                 'the statistical test. Default: 0.05.')

        return parser

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        return cls.prepare(input_dir, target_dir, system_name)

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        args = cls.parser().parse_args(args)
        return cls.analyze(gmx_parser, system_dir, system_name, base_data, verbosity,
                           alpha=args.tolerance)

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name):
        # no additional sims needed, base is enough
        # could check energy writing settings
        return []

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity, alpha=None):
        # Standard value
        if alpha is None:
            alpha = cls.parser().get_default('tolerance')

        # base data
        if base_data['reduced'] is None:
            current_dir = os.path.join(system_dir, 'base')
            base_data['reduced'] = gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.gro')
            )
        base_result = base_data['reduced']

        p = kinetic_energy.mb_ensemble(base_result, verbose=(verbosity > 0))
        # filename=os.path.join(system_dir, system_name + '_mb'))

        if p >= alpha:
            message = 'MaxwellBoltzmannTest PASSED (p = {:g}, alpha = {:f})'.format(p, alpha)
        else:
            message = 'MaxwellBoltzmannTest FAILED (p = {:g}, alpha = {:f})'.format(p, alpha)

        return {'test': p >= alpha,
                'result': p,
                'tolerance': alpha,
                'message': message}


class EquipartitionTest(Test):
    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(
            description='Tests the equipartition of the kinetic energy.',
            formatter_class=argparse.RawTextHelpFormatter,
            prog=cls.__name__
        )
        parser.add_argument('--distribution', default=False, action='store_true',
                            help=('If set, the groups of degrees of freedom will\n'
                                  'be separately tested as Maxwell-Boltzmann\n'
                                  'distributions. Otherwise, only the difference in\n'
                                  'average kinetic energy is checked.'))
        parser.add_argument('-t', '--tolerance', type=float, default=0.1,
                            help=('If --distribution is set:\n'
                                  '    The alpha value (confidence interval) used in the\n'
                                  '    statistical test. Default: 0.1.\n'
                                  'Else:\n'
                                  '    The maximum deviation allowed between groups.\n'
                                  '    Default: 0.1 (10%%)'))
        parser.add_argument('--random_groups', type=int, default=0,
                            help='Number of random division tests attempted.\n'
                                 'Default: 0 (random division tests off).')
        parser.add_argument('--random_divisions', type=int, default=2,
                            help='Number of groups the system is randomly divided in.\n'
                                 'Default: 2.')
        parser.add_argument('--molec_group', type=int, nargs='*', default=None, action='append',
                            help=('List of molecule indeces defining a group. Useful to\n'
                                  'pre-define groups of molecules\n'
                                  '(e.g. solute / solvent, liquid mixture species, ...).\n'
                                  'If not set, no pre-defined molecule groups will be\ntested.\n'
                                  'Note: If the last --molec-group flag is given empty,\n'
                                  'the remaining molecules are collected in this group.\n'
                                  'This allows, for example, to only specify the solute,\n'
                                  'and indicate the solvent by giving an empty flag.'))

        return parser

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name):
        # no additional sims needed, base is enough
        # could check position, velocity & energy writing settings
        return []

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        return cls.prepare(input_dir, target_dir, system_name)

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        args = cls.parser().parse_args(args)
        return cls.analyze(gmx_parser, system_dir, system_name, base_data, verbosity,
                           tolerance=args.tolerance, distribution=args.distribution,
                           random_groups=args.random_groups, random_divisions=args.random_divisions,
                           molec_group=args.molec_group)

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity,
                tolerance=None, distribution=None,
                random_groups=None, random_divisions=None,
                molec_group=None):

        # Standard values
        if tolerance is None:
            tolerance = cls.parser().get_default('tolerance')
        if distribution is None:
            distribution = cls.parser().get_default('distribution')
        if random_groups is None:
            random_groups = cls.parser().get_default('random_groups')
        if random_divisions is None:
            random_divisions = cls.parser().get_default('random_divisions')
        if molec_group is None:
            molec_group = cls.parser().get_default('molec_group')

        # base data
        if base_data['full'] is None:
            current_dir = os.path.join(system_dir, 'base')
            base_data['full'] = gmx_parser.get_simulation_data(
                mdp=os.path.join(current_dir, 'mdout.mdp'),
                top=os.path.join(current_dir, 'system.top'),
                edr=os.path.join(current_dir, 'system.edr'),
                gro=os.path.join(current_dir, 'system.trr')
            )
        base_result = base_data['full']

        result = kinetic_energy.equipartition(base_result,
                                              dtemp=tolerance, alpha=tolerance,
                                              distribution=distribution, molec_groups=molec_group,
                                              random_divisions=random_divisions, random_groups=random_groups,
                                              verbosity=0)

        if distribution:
            res = min(result)
            test = res >= tolerance
            if test:
                message = 'EquipartitionTest PASSED (min p = {:g}, alpha = {:f})'.format(res, tolerance)
            else:
                message = 'EquipartitionTest FAILED (min p = {:g}, alpha = {:f})'.format(res, tolerance)
        else:
            dev_min = min(result)
            dev_max = max(result)
            if (1-dev_min) > (dev_max-1):
                res = dev_min
            else:
                res = dev_max
            test = res <= tolerance
            if test:
                message = 'EquipartitionTest PASSED (max dev = {:g}, tolerance = {:f})'.format(res, tolerance)
            else:
                message = 'EquipartitionTest FAILED (max dev = {:g}, tolerance = {:f})'.format(res, tolerance)

        return {'test': test,
                'result': res,
                'tolerance': tolerance,
                'message': message}


class KinConstraintsTest(Test):
    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(
            description='Tests whether there is kinetic energy in constrained degrees of\nfreedom.',
            formatter_class=argparse.RawTextHelpFormatter,
            prog=cls.__name__
        )
        parser.add_argument('-t', '--tolerance', type=float, default=0.1,
                            help='Tolerance - TODO. Default: 0.1.')

        return parser

    @classmethod
    def prepare(cls, input_dir, target_dir, system_name):
        # no additional sims needed, base is enough
        # could check if there are any constraints in the system
        return []

    @classmethod
    def prepare_parser(cls, input_dir, target_dir, system_name, nobackup, args):
        return cls.prepare(input_dir, target_dir, system_name)

    @classmethod
    def analyze_parser(cls, gmx_parser, system_dir, system_name, base_data, verbosity, args):
        raise NotImplementedError

    @classmethod
    def analyze(cls, gmx_parser, system_dir, system_name, base_data, verbosity):
        raise NotImplementedError


all_tests = OrderedDict([
    ('integrator', IntegratorTest),
    ('ensemble', EnsembleTest),
    ('kin_mb', MaxwellBoltzmannTest),
    ('kin_equipartition', EquipartitionTest),
    ('kin_constraints', KinConstraintsTest)
])


def parse_systems(systems_json, systems_user, source_path):
    # Parse json
    # As the order of the systems and the tests
    # might be meaningful, we need ordered dicts!
    system_list = json.load(systems_json)
    system_dict = OrderedDict()
    for system in system_list:
        system_name = system['dir']
        # do the input files exist?
        input_dir = os.path.join(source_path, system_name, 'input')
        if not (os.path.isdir(input_dir) and
                os.path.exists(os.path.join(input_dir, 'system.mdp')) and
                os.path.exists(os.path.join(input_dir, 'system.top')) and
                os.path.exists(os.path.join(input_dir, 'system.gro'))):
            raise ValueError('System ' + system_name + ' in ' +
                             systems_json.name + ': Input files not found')
        # no need to run systems that we don't test
        if 'tests' not in system:
            raise ValueError('System ' + system_name + ' in ' +
                             systems_json.name + ' has no tests defined')
        test_list = system['tests']
        system['tests'] = OrderedDict()
        for test in test_list:
            test_name = test['test']
            if test_name not in all_tests:
                raise ValueError('Test ' + test_name + ' in ' +
                                 systems_json.name + ' is not a valid test')
            if test_name not in system['tests']:
                if 'args' in test:
                    test['args'] = [test['args'].split()]
                else:
                    test['args'] = [[]]
                system['tests'][test_name] = test
            else:
                if 'args' in test:
                    system['tests'][test_name]['args'].append(test['args'].split())
                else:
                    system['tests'][test_name]['args'].append([])

        system_dict[system_name] = system
        # add standard arguments
        if 'grompp_args' not in system:
            system['grompp_args'] = []
        else:
            system['grompp_args'] = system['grompp_args'].split()
        if 'mdrun_args' not in system:
            system['mdrun_args'] = []
        else:
            system['mdrun_args'] = system['mdrun_args'].split()

    # if user has not chosen otherwise, return full dict of systems
    if not systems_user:
        return system_dict

    # make sure user_dicts matches at least something
    for user_system in systems_user:
        for system in system_dict:
            if re.match(user_system + '$', system):
                break
        else:
            raise ValueError('System ' + user_system +
                             ' used in command line argument is not defined in ' +
                             systems_json.name)

    for system in list(system_dict.keys()):
        # delete systems not selected by user
        for user_system in systems_user:
            if re.match(user_system + '$', system):
                user_key = user_system
                break
        else:
            system_dict.pop(system)

    # return reduced dict of systems
    return system_dict


def main(args):
    parser = argparse.ArgumentParser(
        description='Physical validation suite for GROMACS.',
        prog='gmx_physicalvalidation.py',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Use --tests for details about the available tests and their arguments.'
    )
    parser.add_argument('json', type=open,
                        metavar='systems.json',
                        help='Json file containing systems and tests to be ran.')
    parser.add_argument('--tests', default=False, action='store_true',
                        help='Print details about the available tests and their arguments and exit.')
    parser.add_argument('-v', '--verbosity', type=int,
                        metavar='v', default=0,
                        help='Verbosity level. Default: 0 (quiet).')
    parser.add_argument('--gmx', type=str, metavar='exe', default=None,
                        help=('GROMACS executable. Default: Trying to use \'gmx\'.\n' +
                              'Note: If -p is used, the GROMACS executable is not needed.'))
    parser.add_argument('--bindir', type=str, metavar='dir', default=None,
                        help=('GROMACS binary directory.\n' +
                              'If set, trying to use \'bindir/gmx\' instead of plain \'gmx\'\n' +
                              'Note: If --gmx is set, --bindir and --suffix are ignored.'))
    parser.add_argument('--suffix', type=str, metavar='_s', default=None,
                        help=('Suffix of the GROMACS executable.\n' +
                              'If set, trying to use \'gmx_s\' instead of plain \'gmx\'\n' +
                              'Note: If --gmx is set, --bindir and --suffix are ignored.'))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', '--prepare', action='store_true',
                       default=False,
                       help=('Only prepare simulations and output a \'run.sh\' file\n' +
                             'containing the necessary commands to run the systems.\n' +
                             'This allows to separate running simulations from analyzing them,\n' +
                             'useful e.g. to analyze the results on a different machine.\n' +
                             'Default: If none of \'-p\', \'-r\' or \'-a\' is given,\n' +
                             '         the systems are prepared, ran and analyzed in one call.'))
    group.add_argument('-r', '--run', action='store_true',
                       default=False,
                       help=('Only prepare and run simulations.\n' +
                             'Default: If none of \'-p\', \'-r\' or \'-a\' is given,\n' +
                             '         the systems are prepared, ran and analyzed in one call.'))
    group.add_argument('-a', '--analyze', action='store_true',
                       default=False,
                       help=('Only analyze previously ran simulations.\n' +
                             'This requires that the systems have been prepared using this program.\n' +
                             'Default: If none of \'-p\', \'-r\' or \'-a\' is given,\n' +
                             '         the systems are prepared, ran and analyzed in one call.'))
    parser.add_argument('-s', '--system', action='append', dest='systems',
                        metavar='system',
                        help=('Specify which system to run.\n' +
                              '\'system\' needs to match a system defined in the json file.\n' +
                              'Several systems can be specified by chaining several \'-s\' arguments.\n' +
                              'Defaults: If no \'-s\' argument is given, all systems and tests\n' +
                              '          defined in the json file are ran.\n' +
                              'Note: \'system\' can be a regular expression matching more than one system.'))
    parser.add_argument('--mpicmd', type=str, metavar='cmd', default=None,
                        help='MPI command used to invoke run command')
    parser.add_argument('--wd', '--working_dir', type=str,
                        metavar='dir', default=None,
                        help='Working directory (default: current directory)')
    parser.add_argument('--nobackup', default=False, action='store_true',
                        help='Do not create backups of files or folders.')
    
    if '--tests' in args:
        message = ('Physical validation suite for GROMACS\n'
                   'Available tests and their arguments\n'
                   '=====================================\n\n'
                   'The following tests can be specified in the json file:\n\n')
        for test in all_tests:
            message += '  * ' + test + '\n'

        message += '\nAll tests accept additional arguments:\n'
    
        for test, test_cls in all_tests.items():
            message += '\n'
            test_help = test_cls.parser().format_help()
            test_help = test_help.replace(test_cls.__name__, test).replace('usage: ', '')
            test_help = test_help.replace('  -h, --help            show this help message and exit\n', '')
            test_help = test_help.replace(' [-h]', '')
            test_help = test_help.replace(' ' * (len(test_cls.__name__) + 7), ' ' * len(test))
    
            first_line = test_help.split('\n')[0]
            separator = '-' * 70
            message += test_help.replace(first_line, separator + '\n' + first_line)

        sys.stderr.write(message)
        return

    args = parser.parse_args(args)

    # the input files are expected to be located where this script is
    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'systems')
    # target directory can be current or user chosen
    if args.wd is None:
        target_path = os.getcwd()
    else:
        if not os.path.exists(args.wd):
            os.makedirs(args.wd)
        target_path = args.wd

    # get ordered dict of systems from combination of json file and user choices
    systems = parse_systems(args.json, args.systems, source_path)

    # parse simulation stage to perform
    do_all = not (args.prepare or args.run or args.analyze)
    do_prepare = do_all or args.prepare or args.run
    write_script = args.prepare
    do_run = do_all or args.run
    do_analysis = do_all or args.analyze

    # prepare GROMACS interface
    if args.gmx:
        gmx = args.gmx
    else:
        gmx = 'gmx'
        if args.suffix:
            gmx += args.suffix
        if args.bindir:
            gmx = os.path.join(args.bindir, gmx)
    gmx_interface = None
    gmx_parser = None
    if do_run or do_analysis:
        gmx_interface = GromacsInterface(exe=gmx)
        gmx_parser = GromacsParser(exe=gmx)

    if do_prepare:
        runs = []  # this will contain all information needed to run the system
        for system_name, system in systems.items():
            system_dirs = []  # list of directories with subsystems
            # prepare the base system
            input_dir = os.path.join(source_path, system_name, 'input')
            target_dir = os.path.join(target_path, system_name)
            mkdir_bk(target_dir, nobackup=args.nobackup)
            basedir = os.path.join(target_dir, 'base')
            mkdir_bk(basedir, nobackup=args.nobackup)
            for suffix in ['mdp', 'gro', 'top']:
                shutil.copy2(os.path.join(input_dir, 'system.' + suffix),
                             basedir)
            system_dirs.append(basedir)
            # call prepare method of chosen tests
            for test_name, test in system['tests'].items():
                for test_args in test['args']:
                    system_dirs.extend(
                        all_tests[test_name].prepare_parser(input_dir, target_dir, system_name,
                                                            args.nobackup, test_args)
                    )

            # save run information
            for d in system_dirs:
                runs.append({
                    'dir': d,
                    'grompp_args': system['grompp_args'],
                    'mdrun_args': system['mdrun_args']
                })
        # end of loop over systems

        if write_script:
            script_file = os.path.join(target_path, 'run_simulations.sh')
            if not args.nobackup:
                file_bk(script_file)
            with open(script_file, 'w') as f:
                f.write('# This file was created by the physical validation suite for GROMACS.\n')
                f.write('\n# Define run variables\n')
                f.write('WORKDIR=' + target_path + '\n')
                f.write('GROMPPCMD="' + gmx + ' grompp"\n')
                f.write('MDRUNCMD="' + gmx + ' mdrun"\n')
                f.write('\n# Run systems\n')
                f.write('startpath=$PWD\n')
                f.write('cd $WORKDIR\n')
                for run in runs:
                    for cmd in basic_run_cmds(directory=run['dir'],
                                              grompp_args=run['grompp_args'],
                                              mdrun_args=run['mdrun_args']):
                        f.write(cmd + '\n')
                    f.write('\n')
                f.write('cd $startpath\n')
        # end if write_script

        if do_run:
            # send messages from GROMACS to log
            gmx_log = open(os.path.join(target_path, 'physicalvalidation_gmx.log'), 'w')
            for run in runs:
                gmx_interface.grompp(mdp='system.mdp',
                                     top='system.top',
                                     gro='system.gro',
                                     tpr='system.tpr',
                                     cwd=run['dir'],
                                     args=run['grompp_args'],
                                     stdout=gmx_log,
                                     stderr=gmx_log)
                gmx_interface.mdrun(tpr='system.tpr',
                                    deffnm='system',
                                    cwd=run['dir'],
                                    args=run['mdrun_args'],
                                    stdout=gmx_log,
                                    stderr=gmx_log,
                                    mpicmd=args.mpicmd)
            gmx_log.close()
        # end if do_run
    # end if do_prepare

    if do_analysis:
        title = 'GROMACS PHYSICAL VALIDATION RESULTS'
        width = 70
        indent = int((width - len(title))/2)
        print()
        print(' ' * indent + '=' * len(title))
        print(' ' * indent + title)
        print(' ' * indent + '=' * len(title))
        print()
        passed = True
        for system_name, system in systems.items():
            # save system data if re-used for different test
            # massively reduces run time of multiple tests
            system_data = {
                'reduced': None,
                'full': None
            }
            # system directory
            target_dir = os.path.join(target_path, system_name)

            print('Analyzing system ' + system_name)

            # call analyze method of chosen tests
            for test_name, test in system['tests'].items():
                for test_args in test['args']:
                    try:
                        result = all_tests[test_name].analyze_parser(gmx_parser, target_dir,
                                                                     system_name, system_data,
                                                                     args.verbosity, test_args)
                    except Exception as err:
                        print('    ' + all_tests[test_name].__name__ + ' FAILED (Exception in evaluation)')
                        print('    '*2 + type(err).__name__ + ': ' + str(err))
                    else:
                        for line in result['message'].split('\n'):
                            print('    ' + line)

                        passed = passed and result['test']
            # end loop over tests
            print()
        # end loop over systems
        return int(not passed)
    # end if do_analysis

    # assuming everything is ok if we ended up here
    return 0


if __name__ == "__main__":
    return_value = main(sys.argv[1:])
    sys.exit(return_value)
