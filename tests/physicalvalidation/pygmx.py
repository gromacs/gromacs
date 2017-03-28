# py2 compatibility
from __future__ import absolute_import, division, print_function
# other imports
import subprocess
import os
import sys
import re
import numpy as np
import itertools
import copy
import shutil
import errno


class PyGMX(object):
    def __init__(self, exe=None, dp=None):

        self._exe = None
        self._dp = True

        if exe is None:
            if self.check_exe(quiet=True, exe='gmx'):
                self.exe = 'gmx'
            else:
                print('WARNING: gmx executable not found. Set before attempting to run!')
        else:
            self.exe = exe

        if dp is not None:
            self.dp = dp

    @property
    def exe(self):
        """exe is a string pointing to the gmx executable."""
        return self._exe

    @exe.setter
    def exe(self, exe):
        if self.check_exe(exe=exe):
            if os.path.dirname(exe):
                exe = os.path.abspath(exe)
            self._exe = exe

    @property
    def double(self):
        return self._dp

    @double.setter
    def double(self, dp):
        assert isinstance(dp, bool)
        self._dp = dp

    def check_exe(self, quiet=False, exe=None):
        # could we also check that it is actually Gromacs, not just any existing executable?
        if exe is None:
            exe = self._exe
        try:
            devnull = open(os.devnull)
            subprocess.call([exe], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                # file not found error.
                if not quiet:
                    print('ERROR: gmx executable not found')
                    print(exe)
                return False
        return True

    def run(self, cmd, args, cwd=None, stdin=None, stdout=None, stderr=None):
        if self.exe is None:
            print('ERROR: No gmx executable defined. Set before attempting to run!')
        command = [self.exe, cmd]
        command.extend(args)
        return subprocess.Popen(command, cwd=cwd,
                                stdin=stdin, stdout=stdout, stderr=stderr)

    def grompp(self, mdp, top, gro, tpr=None,
               cwd='.', maxwarn=0, args=None,
               stdin=None, stdout=None, stderr=None):
        cwd = os.path.abspath(cwd)
        assert os.path.exists(os.path.join(cwd, mdp))
        assert os.path.exists(os.path.join(cwd, top))
        assert os.path.exists(os.path.join(cwd, gro))

        if args is None:
            args = []

        if tpr is None:
            tpr = os.path.basename(mdp).replace('.mdp', '') + '.tpr'
        else:
            assert os.path.exists(os.path.join(cwd, os.path.dirname(tpr)))

        args = ['-f', mdp, '-p', top, '-c', gro, '-o', tpr, '-maxwarn', str(maxwarn)] + args
        proc = self.run('grompp', args, cwd=cwd,
                        stdin=stdin, stdout=stdout, stderr=stderr)
        proc.wait()
        return proc.returncode

    def mdrun(self, tpr, edr=None, deffnm=None, cwd='.', args=None,
              stdin=None, stdout=None, stderr=None):
        cwd = os.path.abspath(cwd)
        tpr = os.path.abspath(tpr)
        assert os.path.exists(tpr)
        assert os.path.exists(cwd)

        if args is None:
            args = []

        if deffnm is None:
            deffnm = os.path.basename(tpr).replace('.tpr', '')

        args = ['-s', tpr, '-deffnm', deffnm] + args
        if edr is not None:
            args += ['-e', edr]
        proc = self.run('mdrun', args, cwd=cwd,
                        stdin=stdin, stdout=stdout, stderr=stderr)
        proc.wait()
        return proc.returncode

    def get_quantities(self, edr, quantities, cwd=None,
                       begin=None, end=None, args=None):

        if args is None:
            args = []

        tmp_xvg = 'gmxpy_' + os.path.basename(edr).replace('.edr', '') + '.xvg'
        if cwd is not None:
            tmp_xvg = os.path.join(cwd, tmp_xvg)

        q_dict = {}

        for q in quantities:
            not_found = self.create_xvg(edr, tmp_xvg, [q], cwd=cwd,
                                        begin=begin, end=end, args=args)[1]
            if q in not_found:
                q_dict[q] = None
                continue

            skip_line = re.compile("^[#,@]")
            values = []
            times = []
            with open(tmp_xvg, 'r') as xvg:
                for line in xvg:
                    if skip_line.match(line):
                        continue
                    times.append(float(line.split()[0]))
                    values.append(float(line.split()[1]))

            if 'time' in q_dict:
                if not np.array_equal(np.array(times), q_dict['time']):
                    print('WARNING: Time discrepancy in ' + edr)
            else:
                q_dict['time'] = np.array(times)
            q_dict[q] = np.array(values)

        return q_dict

    def create_xvg(self, edr, xvg, quantities, cwd=None,
                   begin=None, end=None, args=None):
        assert os.path.exists(edr)
        assert os.path.exists(os.path.abspath(os.path.dirname(xvg)))

        if args is None:
            args = []
        
        if self._dp:
            args.append('-dp')
        if begin is not None:
            args.extend(['-b', str(begin)])
        if end is not None:
            args.extend(['-e', str(end)])

        quants = ''
        for q in quantities:
            quants += str(q) + '\n'

        args = ['-f', edr, '-o', xvg] + args
        proc = self.run('energy', args, cwd=cwd,
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        err = proc.communicate(quants.encode(sys.stdin.encoding))[1]

        encoding = sys.stderr.encoding
        if encoding is None:
            encoding = 'UTF-8'

        not_found = []
        if 'does not match anything' in err.decode(encoding):
            for q in quantities:
                if "String '" + q + "' does not match anything" in err.decode(encoding):
                    not_found.append(q)

        return proc.wait(), not_found


class System(object):
    """ Helper class holding all variables of a system needed by the
    physical validation suite."""
    def __init__(self):
        self._cwd = '.'
        self._mdp = None
        self._top = None
        self._gro = None
        self._tpr = None
        self._edr = None
        self._deffnm = None
        self._quantities = None
        self._values = None
        self._parameters = None

    @classmethod
    def from_deffnm(cls, deffnm, cwd):
        assert isinstance(deffnm, str)
        cwd = os.path.abspath(cwd)
        assert os.path.exists(cwd)

        obj = cls()

        obj.cwd = cwd
        obj.deffnm = deffnm

        mdp = deffnm + '.mdp'
        top = deffnm + '.top'
        gro = deffnm + '.gro'
        tpr = deffnm + '.tpr'
        edr = deffnm + '.edr'
        
        obj.mdp = mdp
        obj.top = top
        obj.gro = gro
        obj.tpr = tpr
        obj.edr = edr
        
        return obj
        
    @property
    def cwd(self):
        return self._cwd
    
    @cwd.setter
    def cwd(self, cwd):
        cwd = os.path.abspath(cwd)
        assert os.path.exists(cwd)
        self._cwd = cwd
        
    @property
    def mdp(self):
        return os.path.join(self._cwd, self._mdp)
    
    @mdp.setter
    def mdp(self, mdp):
        assert os.path.exists(os.path.join(self.cwd, mdp))
        self._mdp = mdp
        
    @property
    def top(self):
        return os.path.join(self._cwd, self._top)
    
    @top.setter
    def top(self, top):
        assert os.path.exists(os.path.join(self.cwd, top))
        self._top = top
        
    @property
    def gro(self):
        return os.path.join(self._cwd, self._gro)
    
    @gro.setter
    def gro(self, gro):
        assert os.path.exists(os.path.join(self.cwd, gro))
        self._gro = gro
        
    @property
    def tpr(self):
        return os.path.join(self._cwd, self._tpr)
    
    @tpr.setter
    def tpr(self, tpr):
        assert os.path.exists(os.path.join(self.cwd, os.path.dirname(tpr)))
        self._tpr = tpr
        
    @property
    def edr(self):
        return os.path.join(self._cwd, self._edr)
    
    @edr.setter
    def edr(self, edr):
        assert os.path.exists(os.path.join(self.cwd, os.path.dirname(edr)))
        self._edr = edr

    @property
    def deffnm(self):
        return self._deffnm

    @deffnm.setter
    def deffnm(self, deffnm):
        assert isinstance(deffnm, str)
        self._deffnm = deffnm

    @property
    def quantities(self):
        return self._quantities

    @quantities.setter
    def quantities(self, quantities):
        assert isinstance(quantities, list) or isinstance(quantities, str)
        if isinstance(quantities, str):
            quantities = [quantities]
        self._quantities = quantities

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        assert isinstance(values, dict)
        self._values = values

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, options):
        assert isinstance(options, dict)
        self._parameters = options

    def run_gromacs(self, gmx, verb=2):

        out_file = self.mdp.replace('.mdp', '.stdout')
        err_file = self.mdp.replace('.mdp', '.stderr')
        devnull = open(os.devnull, 'w')

        if verb <= 0:
            stdout = devnull
            stderr = devnull
        elif verb == 1:
            stdout = devnull
            stderr = open(err_file, 'w')
        elif verb == 2:
            stdout = devnull
            stderr = None
        elif verb == 3:
            stdout = open(out_file, 'w')
            stderr = open(err_file, 'w')
        elif verb == 4:
            stdout = open(out_file, 'w')
            stderr = None
        else:
            stdout = None
            stderr = None

        if self._tpr is None:
            self._tpr = self._mdp.replace('.mdp', '.tpr')

        if self._deffnm is None:
            self._deffnm = os.path.basename(self._mdp).replace('.mdp', '')

        gmx.grompp(self.mdp, self.top, self.gro, tpr=self.tpr, cwd=self.cwd, stdout=stdout, stderr=stderr)
        gmx.mdrun(self.tpr, edr=self.edr, deffnm=self.deffnm, cwd=self.cwd, stdout=stdout, stderr=stderr)

        if self._edr is None:
            self.edr = self.deffnm + '.edr'

        if self.quantities is not None:
            self.values = gmx.get_quantities(self.edr, self.quantities, cwd=self.cwd)


def read_mdp(mdp, sep=None):
    """Reads an mdp, allowing for multiple options per parameter
    separated by sep."""
    assert os.path.exists(mdp)

    parameters = {}
    options = []

    with open(mdp, 'r') as f:
        for line in f:
            line = line.split(';')[0]
            if '=' not in line:
                continue
            key, value = line.split('=', 1)
            key = key.strip()
            if sep is None or sep not in value:
                parameters[key] = value.strip()
            else:
                parameters[key] = []
                options.append(key)
                for v in value.split(sep):
                    parameters[key].append(v.strip())

    return parameters, options


def write_mdp(file, parameters):
    with open(file, 'w') as f:
        for key in parameters:
            f.write('{:25s} = {:s}\n'.format(key, parameters[key]))


def split_systems(system, parameters, options, dt=None):
    systems = []

    options_dict = {key: parameters[key] for key in options}

    # Make a list of permutations of the options
    list_of_options = [dict(zip(options_dict, combination))
                       for combination in itertools.product(*options_dict.values())]

    for counter, new_opt in enumerate(list_of_options):
        new_system = copy.deepcopy(system)
        new_system.options = new_opt
        new_path = os.path.join(system.cwd, 'opt_' + str(counter))
        try:
            os.mkdir(new_path)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
        new_system.cwd = new_path
        shutil.copy2(system.top, new_system.cwd)
        shutil.copy2(system.gro, new_system.cwd)
        new_param = copy.deepcopy(parameters)
        if 'dt' in new_opt and 'nsteps' not in new_opt and dt is not None:
            new_opt['nsteps'] = str(int(int(parameters['nsteps']) * dt / float(new_opt['dt'])))
        for key in new_opt:
            new_param[key] = new_opt[key]
        new_system.parameters = new_param
        write_mdp(os.path.join(new_system.cwd, new_system.mdp), new_param)
        systems.append(new_system)

    return systems
