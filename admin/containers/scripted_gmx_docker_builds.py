#!/usr/bin/env python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
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

"""
Generates a set of docker images used for running GROMACS CI on Gitlab.
The images are prepared according to a selection of build configuration targets
that hope to cover a broad enough scope of different possible systems,
allowing us to check compiler types and versions, as well as libraries used
for accelerators and parallel communication systems. Each combinations is
described as an entry in the build_configs dictionary, with the script
analysing the logic and adding build stages as needed.

Based on the example script provided by the NVidia HPCCM repository.

Authors:
    * Paul Bauer <paul.bauer.q@gmail.com>
    * Eric Irrgang <ericirrgang@gmail.com>
    * Joe Jordan <e.jjordan12@gmail.com>

Usage::

    $ python3 scripted_gmx_docker_builds.py --help
    $ python3 scripted_gmx_docker_builds.py --format docker > Dockerfile && docker build .
    $ python3 scripted_gmx_docker_builds.py | docker build -

"""

import argparse
import collections
import typing
from distutils.version import StrictVersion

import hpccm
import hpccm.config
from hpccm.building_blocks.base import bb_base

try:
    import utility
except ImportError:
    raise RuntimeError(
        'This module assumes availability of supporting modules in the same directory. Add the directory to '
        'PYTHONPATH or invoke Python from within the module directory so module location can be resolved.')

# Basic packages for all final images.
_common_packages = ['build-essential',
                    'ca-certificates',
                    'ccache',
                    'git',
                    'gnupg',
                    'libfftw3-dev',
                    'libhwloc-dev',
                    'liblapack-dev',
                    'libx11-dev',
                    'moreutils',
                    'ninja-build',
                    'rsync',
                    'valgrind',
                    'vim',
                    'wget',
                    'xsltproc']

# Extra packages needed to build Python installations from source.
_python_extra_packages = ['build-essential',
                          'ca-certificates',
                          'ccache',
                          'curl',
                          'git',
                          'libbz2-dev',
                          'libffi-dev',
                          'liblzma-dev',
                          'libncurses5-dev',
                          'libncursesw5-dev',
                          'libreadline-dev',
                          'libsqlite3-dev',
                          'libssl-dev',
                          'llvm',
                          'python-openssl',
                          'vim',
                          'wget',
                          'zlib1g-dev']

# Extra packages needed for images for building documentation.
_docs_extra_packages = ['autoconf',
                        'automake',
                        'autopoint',
                        'autotools-dev',
                        'bison',
                        'flex',
                        'ghostscript',
                        'graphviz',
                        'help2man',
                        'imagemagick',
                        'libtool',
                        'linkchecker',
                        'mscgen',
                        'm4',
                        'texinfo',
                        'texlive-latex-base',
                        'texlive-latex-extra',
                        'texlive-fonts-recommended',
                        'texlive-fonts-extra']

# Supported Python versions for maintained branches.
# TODO: Remove '3.5.9' from defaults in master once script in release-2020 diverges.
_python_versions = ['3.5.9', '3.6.10', '3.7.7', '3.8.2']

# Parse command line arguments
parser = argparse.ArgumentParser(description='GROMACS CI image creation script', parents=[utility.parser])

parser.add_argument('--format', type=str, default='docker',
                    choices=['docker', 'singularity'],
                    help='Container specification format (default: docker)')
parser.add_argument('--venvs', nargs='*', type=str, default=_python_versions,
                    help='List of Python versions ("major.minor.patch") for which to install venvs. '
                         'Default: {}'.format(' '.join(_python_versions)))


def base_image_tag(args) -> str:
    # Check if we use CUDA images or plain linux images
    if args.cuda is not None:
        cuda_version_tag = 'nvidia/cuda:' + args.cuda + '-devel'
        if args.centos is not None:
            cuda_version_tag += '-centos' + args.centos
        elif args.ubuntu is not None:
            cuda_version_tag += '-ubuntu' + args.ubuntu
        else:
            raise RuntimeError('Logic error: no Linux distribution selected.')

        base_image_tag = cuda_version_tag
    else:
        if args.centos is not None:
            base_image_tag = 'centos:centos' + args.centos
        elif args.ubuntu is not None:
            base_image_tag = 'ubuntu:' + args.ubuntu
        else:
            raise RuntimeError('Logic error: no Linux distribution selected.')
    return base_image_tag


def get_llvm_packages(args) -> typing.Iterable[str]:
    # If we use the package version of LLVM, we need to install extra packages for it.
    if (args.llvm is not None) and (args.tsan is None):
        return ['libomp-dev',
                'clang-format-' + str(args.llvm),
                'clang-tidy-' + str(args.llvm)]
    else:
        return []


def get_compiler(args, tsan_stage: hpccm.Stage = None) -> bb_base:
    # Compiler
    if args.icc is not None:
        raise RuntimeError('Intel compiler toolchain recipe not implemented yet')

    if args.llvm is not None:
        # Build our own version instead to get TSAN + OMP
        if args.tsan is not None:
            if tsan_stage is not None:
                compiler = tsan_stage.runtime(_from='tsan')
            else:
                raise RuntimeError('No TSAN stage!')
        # Build the default compiler if we don't need special support
        else:
            compiler = hpccm.building_blocks.llvm(extra_repository=True, version=args.llvm)

    elif (args.gcc is not None):
        compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gcc,
                                             fortran=False)
    else:
        raise RuntimeError('Logic error: no compiler toolchain selected.')
    return compiler


def get_mpi(args, compiler):
    # If needed, add MPI to the image
    if args.mpi is not None:
        if args.mpi == 'openmpi':
            use_cuda = False
            if args.cuda is not None:
                use_cuda = True

            if hasattr(compiler, 'toolchain'):
                return hpccm.building_blocks.openmpi(toolchain=compiler.toolchain, cuda=use_cuda, infiniband=False)
            else:
                raise RuntimeError('compiler is not an HPCCM compiler building block!')

        elif args.mpi == 'impi':
            raise RuntimeError('Intel MPI recipe not implemented yet.')
        else:
            raise RuntimeError('Requested unknown MPI implementation.')
    else:
        return None


def get_opencl(args):
    # Add OpenCL environment if needed
    if (args.opencl is not None):
        if args.opencl == 'nvidia':
            if (args.cuda is None):
                raise RuntimeError('Need Nvidia environment for Nvidia OpenCL image')

            return hpccm.building_blocks.packages(ospackages=['nvidia-opencl-dev'])

        elif args.opencl == 'intel':
            return hpccm.building_blocks.packages(
                    apt_ppas=['ppa:intel-opencl/intel-opencl'],
                    ospackages=['opencl-headers', 'ocl-icd-libopencl1',
                                'ocl-icd-opencl-dev', 'intel-opencl-icd'])

        elif args.opencl == 'amd':
            # libelf1 is a necessary dependency for something in the ROCm stack,
            # which they should set up, but seem to have omitted.
            return hpccm.building_blocks.packages(
                    apt_keys=['http://repo.radeon.com/rocm/apt/debian/rocm.gpg.key'],
                    apt_repositories=['deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ xenial main'],
                    ospackages=['ocl-icd-libopencl1', 'ocl-icd-opencl-dev', 'opencl-headers', 'libelf1', 'rocm-opencl'])
    else:
        return None


def get_clfft(args):
    if (args.clfft is not None):
        return hpccm.building_blocks.generic_cmake(
            repository='https://github.com/clMathLibraries/clFFT.git',
            prefix='/usr/local', recursive=True, branch=args.clfft, directory='clFFT/src')
    else:
        return None


def add_tsan_stage(input_args, output_stages: typing.Mapping[str, hpccm.Stage]):
    """Isolate the expensive TSAN preparation stage.

    This is a very expensive stage, but has few and disjoint dependencies, and
    its output is easily compartmentalized (/usr/local) so we can isolate this
    build stage to maximize build cache hits and reduce rebuild time, bookkeeping,
    and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError('Need output_stages container.')
    tsan_stage = hpccm.Stage()
    tsan_stage += hpccm.primitives.baseimage(image=base_image_tag(input_args), _as='tsan')

    tsan_stage += hpccm.building_blocks.packages(ospackages=['git', 'ca-certificates', 'build-essential', 'cmake'])
    # CMake will get duplicated later, but this is an expensive image, and it isn't worth optimizing
    # out that duplication...
    tsan_stage += hpccm.building_blocks.python(python3=True, python2=False, devel=False)

    compiler_branch = 'release_' + str(input_args.llvm) + '0'
    tsan_stage += hpccm.building_blocks.generic_cmake(
        repository='https://git.llvm.org/git/llvm.git',
        prefix='/usr/local', recursive=True, branch=compiler_branch,
        cmake_opts=['-D CMAKE_BUILD_TYPE=Release', '-D LLVM_ENABLE_PROJECTS="clang;openmp;clang-tools-extra"',
                    '-D LIBOMP_TSAN_SUPPORT=on'],
        preconfigure=['export branch=' + compiler_branch,
                      '(cd projects; git clone --depth=1 --branch $branch https://git.llvm.org/git/libcxx.git)',
                      '(cd projects; git clone --depth=1 --branch $branch https://git.llvm.org/git/libcxxabi.git)',
                      '(cd projects; git clone --depth=1 --branch $branch https://git.llvm.org/git/compiler-rt.git)',
                      '(cd ..; git clone --depth=1 --branch $branch https://git.llvm.org/git/openmp.git)',
                      '(cd ..; git clone --depth=1 --branch $branch https://git.llvm.org/git/clang.git)',
                      '(cd ..; git clone --depth=1 --branch $branch https://git.llvm.org/git/clang-tools-extra.git)'],
        postinstall=['ln -s /usr/local/bin/clang++ /usr/local/bin/clang++-' + str(input_args.llvm),
                     'ln -s /usr/local/bin/clang-format /usr/local/bin/clang-format-' + str(input_args.llvm),
                     'ln -s /usr/local/bin/clang-tidy /usr/local/bin/clang-tidy-' + str(input_args.llvm),
                     'ln -s /usr/local/libexec/c++-analyzer /usr/local/bin/c++-analyzer-' + str(input_args.llvm)])
    output_stages['tsan'] = tsan_stage


def prepare_venv(version: StrictVersion) -> typing.Sequence[str]:
    """Get shell commands to set up the venv for the requested Python version."""
    major = version.version[0]
    minor = version.version[1]

    pyenv = '$HOME/.pyenv/bin/pyenv'

    py_ver = '{}.{}'.format(major, minor)
    venv_path = '$HOME/venv/py{}'.format(py_ver)
    commands = ['$({pyenv} prefix `{pyenv} whence python{py_ver}`)/bin/python -m venv {path}'.format(
        pyenv=pyenv,
        py_ver=py_ver,
        path=venv_path
    )]

    commands.append('{path}/bin/python -m pip install --upgrade pip setuptools'.format(
        path=venv_path
    ))
    # Install dependencies for building and testing gmxapi Python package.
    # WARNING: Please keep this list synchronized with python_packaging/requirements-test.txt
    # TODO: Get requirements.txt from an input argument.
    commands.append("""{path}/bin/python -m pip install --upgrade \
            'cmake>=3.9.6' \
            'flake8>=3.7.7' \
            'mpi4py>=2' \
            'networkx>=2.0' \
            'numpy>=1' \
            'pip>=10.1' \
            'pytest>=3.9' \
            'setuptools>=28.0.0' \
            'scikit-build>=0.7'""".format(path=venv_path))

    return commands


def add_python_stages(building_blocks: typing.Mapping[str, bb_base],
                      input_args,
                      output_stages: typing.MutableMapping[str, hpccm.Stage]):
    """Add the stage(s) necessary for the requested venvs.

    One intermediate build stage is created for each venv (see --venv option).

    Each stage partially populates Python installations and venvs in the home
    directory. The home directory is collected by the 'pyenv' stage for use by
    the main build stage.
    """
    if len(input_args.venvs) < 1:
        raise RuntimeError('No venvs to build...')
    if output_stages is None or not isinstance(output_stages, collections.abc.Mapping):
        raise RuntimeError('Need a container for output stages.')

    # Main Python stage that collects the environments from individual stages.
    # We collect the stages individually, rather than chaining them, because the
    # copy is a bit slow and wastes local Docker image space for each filesystem
    # layer.
    pyenv_stage = hpccm.Stage()
    pyenv_stage += hpccm.primitives.baseimage(image=base_image_tag(input_args), _as='pyenv')
    pyenv_stage += building_blocks['compiler']
    pyenv_stage += building_blocks['mpi']
    pyenv_stage += hpccm.building_blocks.packages(ospackages=_python_extra_packages)

    for version in [StrictVersion(py_ver) for py_ver in sorted(input_args.venvs)]:
        stage_name = 'py' + str(version)
        stage = hpccm.Stage()
        stage += hpccm.primitives.baseimage(image=base_image_tag(input_args), _as=stage_name)
        stage += building_blocks['compiler']
        stage += building_blocks['mpi']
        stage += hpccm.building_blocks.packages(ospackages=_python_extra_packages)

        # TODO: Use a non-root user for testing and Python virtual environments.
        stage += hpccm.primitives.shell(commands=[
            'curl https://pyenv.run | bash',
            """echo 'export PYENV_ROOT="$HOME/.pyenv"' >> $HOME/.bashrc""",
            """echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> $HOME/.bashrc""",
            """echo 'eval "$(pyenv init -)"' >> $HOME/.bashrc""",
            """echo 'eval "$(pyenv virtualenv-init -)"' >> $HOME/.bashrc"""])
        pyenv = '$HOME/.pyenv/bin/pyenv'
        commands = ['PYTHON_CONFIGURE_OPTS="--enable-shared" {pyenv} install -s {version}'.format(
            pyenv=pyenv,
            version=str(version))]
        stage += hpccm.primitives.shell(commands=commands)

        commands = prepare_venv(version)
        stage += hpccm.primitives.shell(commands=commands)

        # TODO: Update user home directory.
        pyenv_stage += hpccm.primitives.copy(_from=stage_name, _mkdir=True, src=['/root/'],
                                             dest='/root')

        # Add the intermediate build stage to the sequence
        output_stages[stage_name] = stage

    # TODO: If we activate pyenv for login shells, the `global` "version" should be full-featured.
    # # `version` should be a system installation or pyenv environment (or pyenv-virtualenv)
    # # with the dependencies for all of the Python aspects of CMake-driven builds.
    # commands = '{pyenv} global {version}'.format(
    #             pyenv=pyenv,
    #             version=...)
    # pyenv_stage += hpccm.primitives.shell(commands=commands)

    # Add the aggregating build stage to the sequence. This allows the main stage to copy
    # the files in a single stage, potentially reducing the overall output image size.
    output_stages['pyenv'] = pyenv_stage


def build_stages(args) -> typing.Iterable[hpccm.Stage]:
    """Define and sequence the stages for the recipe corresponding to *args*."""

    # A Dockerfile or Singularity recipe can have multiple build stages.
    # The main build stage can copy files from previous stages, though only
    # the last stage is included in the tagged output image. This means that
    # large or expensive sets of build instructions can be isolated in
    # local/temporary images, but all of the stages need to be output by this
    # script, and need to occur in the correct order, so we create a sequence
    # object early in this function.
    stages = collections.OrderedDict()

    # If we need the TSAN compilers, the early build is more involved.
    if args.llvm is not None and args.tsan is not None:
        add_tsan_stage(input_args=args, output_stages=stages)

    # Building blocks are chunks of container-builder instructions that can be
    # copied to any build stage with the addition operator.
    building_blocks = collections.OrderedDict()

    # These are the most expensive and most reusable layers, so we put them first.
    building_blocks['compiler'] = get_compiler(args, tsan_stage=stages.get('tsan'))
    building_blocks['mpi'] = get_mpi(args, building_blocks['compiler'])

    # Install additional packages early in the build to optimize Docker build layer cache.
    os_packages = _common_packages + get_llvm_packages(args)
    if args.doxygen is not None:
        os_packages += _docs_extra_packages
    building_blocks['ospackages'] = hpccm.building_blocks.packages(ospackages=os_packages)

    building_blocks['cmake'] = hpccm.building_blocks.cmake(eula=True, version=args.cmake)
    building_blocks['opencl'] = get_opencl(args)
    building_blocks['clfft'] = get_clfft(args)

    # Add Python environments to MPI images, only, so we don't have to worry
    # about whether to install mpi4py.
    if args.mpi is not None and len(args.venvs) > 0:
        add_python_stages(building_blocks=building_blocks, input_args=args, output_stages=stages)

    # Create the stage from which the targeted image will be tagged.
    stages['main'] = hpccm.Stage()

    stages['main'] += hpccm.primitives.baseimage(image=base_image_tag(args))
    for bb in building_blocks.values():
        if bb is not None:
            stages['main'] += bb

    # We always add Python3 and Pip
    stages['main'] += hpccm.building_blocks.python(python3=True, python2=False, devel=True)
    stages['main'] += hpccm.building_blocks.pip(upgrade=True, pip='pip3',
                                                packages=['pytest', 'networkx', 'numpy'])

    # Add documentation requirements (doxygen and sphinx + misc).
    if (args.doxygen is not None):
        if (args.doxygen == '1.8.5'):
            doxygen_commit = 'ed4ed873ab0e7f15116e2052119a6729d4589f7a'
        else:
            doxygen_commit = 'a6d4f4df45febe588c38de37641513fd576b998f'
        stages['main'] += hpccm.building_blocks.generic_autotools(
            repository='https://github.com/westes/flex.git',
            commit='f7788a9a0ecccdc953ed12043ccb59ca25714018',
            prefix='/tmp/install-of-flex',
            configure_opts=['--disable-shared'],
            preconfigure=['./autogen.sh'])
        stages['main'] += hpccm.building_blocks.generic_autotools(
            repository='https://github.com/doxygen/doxygen.git',
            commit=doxygen_commit,
            prefix='',
            configure_opts=[
                '--flex /tmp/install-of-flex/bin/flex',
                '--static'],
            postinstall=[
                'sed -i \'/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d\' /etc/ImageMagick-6/policy.xml'])
        stages['main'] += hpccm.building_blocks.pip(pip='pip3', packages=['sphinx==1.6.1'])

    if 'pyenv' in stages and stages['pyenv'] is not None:
        stages['main'] += hpccm.primitives.copy(_from='pyenv', _mkdir=True, src=['/root/.pyenv/'],
                                                dest='/root/.pyenv')
        stages['main'] += hpccm.primitives.copy(_from='pyenv', _mkdir=True, src=['/root/venv/'],
                                                dest='/root/venv')
        # TODO: Update user home directory.
        # TODO: If we activate pyenv for login shells, the `global` "version" should be full-featured.
        # stages['main'] += hpccm.primitives.copy(_from='pyenv', src=['/root/.bashrc'],
        #                                         dest='/root/')

    # Note that the list of stages should be sorted in dependency order.
    for build_stage in stages.values():
        if build_stage is not None:
            yield build_stage


if __name__ == '__main__':
    args = parser.parse_args()

    # Set container specification output format
    hpccm.config.set_container_format(args.format)

    container_recipe = build_stages(args)

    # Output container specification
    for stage in container_recipe:
        print(stage)
