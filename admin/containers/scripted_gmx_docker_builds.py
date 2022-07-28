#!/usr/bin/env python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

"""Building block based Dockerfile generation for CI testing images.

Generates a set of docker images used for running GROMACS CI on Gitlab.
The images are prepared according to a selection of build configuration targets
that hope to cover a broad enough scope of different possible systems,
allowing us to check compiler types and versions, as well as libraries used
for accelerators and parallel communication systems. Each combinations is
described as an entry in the build_configs dictionary, with the script
analysing the logic and adding build stages as needed.

Based on the example script provided by the NVidia HPCCM repository.

Reference:
    `NVidia HPC Container Maker <https://github.com/NVIDIA/hpc-container-maker>`__

Authors:
    * Paul Bauer <paul.bauer.q@gmail.com>
    * Eric Irrgang <ericirrgang@gmail.com>
    * Joe Jordan <e.jjordan12@gmail.com>
    * Mark Abraham <mark.j.abraham@gmail.com>
    * Gaurav Garg <gaugarg@nvidia.com>

Usage::

    $ python3 scripted_gmx_docker_builds.py --help
    $ python3 scripted_gmx_docker_builds.py --format docker > Dockerfile && docker build .
    $ python3 scripted_gmx_docker_builds.py | docker build -

See Also:
    :file:`buildall.sh`

"""

import argparse
import collections
import collections.abc
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
                    'gpg-agent',
                    'less',
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

_opencl_extra_packages = [
    'nvidia-opencl-dev',
    # The following require apt_ppas=['ppa:intel-opencl/intel-opencl']
    'intel-opencl-icd',
    'ocl-icd-libopencl1',
    'ocl-icd-opencl-dev',
    'opencl-headers',
]

_rocm_extra_packages = [
    # The following require
    #             apt_keys=['http://repo.radeon.com/rocm/rocm.gpg.key'],
    #             apt_repositories=['deb [arch=amd64] http://repo.radeon.com/rocm/apt/X.Y.Z/ ubuntu main']
    'clinfo',
    'hipfft',
    'libelf1',
    'rocfft',
    'rocfft-dev',
    'rocm-opencl',
    'rocm-dev',
]

# Extra packages required to build CP2K
_cp2k_extra_packages = [
                        'autoconf',
                        'autogen',
                        'automake',
                        'autotools-dev',
                        'bzip2',
                        'less',
                        'libtool',
                        'make',
                        'nano',
                        'patch',
                        'pkg-config',
                        'python',
                        'python-numpy',
                        'python3',
                        'unzip',
                        'xxd',
                        'zlib1g-dev',
                        'libopenblas-dev'
]

# Extra packages needed to build Intel Compute Runtime
_intel_compute_runtime_extra_packages = ['intel-opencl-icd',
                                         'intel-level-zero-gpu',
                                         'level-zero',
                                         'libmfx1']

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
                        'mscgen',
                        'm4',
                        'openssh-client',
                        'texinfo',
                        'texlive-latex-base',
                        'texlive-latex-extra',
                        'texlive-fonts-recommended',
                        'texlive-fonts-extra',
                        'tex-gyre']

# Parse command line arguments
parser = argparse.ArgumentParser(description='GROMACS CI image creation script',
                                 parents=[utility.parser])

parser.add_argument('--format', type=str, default='docker',
                    choices=['docker', 'singularity'],
                    help='Container specification format (default: docker)')


def base_image_tag(args) -> str:
    """Generate *image* for hpccm.baseimage()."""
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


def hpccm_distro_name(args) -> str:
    """Generate *_distro* for hpccm.baseimage().

    Convert the linux distribution variables into something that hpccm
    understands.

    The same format is used by the lower level hpccm.config.set_linux_distro().
    """
    if args.centos is not None:
        name_mapping = {'7': 'centos7',
                        '8': 'centos8'}
        if args.centos in name_mapping:
            hpccm_name = name_mapping[args.centos]
        else:
            raise RuntimeError('Logic error: unsupported CentOS distribution selected.')
    elif args.ubuntu is not None:
        name_mapping = {'20.04': 'ubuntu20',
                        '18.04': 'ubuntu18',
                        '16.04': 'ubuntu16'}
        if args.ubuntu in name_mapping:
            hpccm_name = name_mapping[args.ubuntu]
        else:
            raise RuntimeError('Logic error: unsupported Ubuntu distribution selected.')
    else:
        raise RuntimeError('Logic error: no Linux distribution selected.')
    return hpccm_name


def get_llvm_packages(args) -> typing.Iterable[str]:
    # If we use the package version of LLVM, we need to install extra packages for it.
    if (args.llvm is not None) and (args.tsan is None):
        packages = [f'libomp-{args.llvm}-dev',
                    f'libomp5-{args.llvm}',
                    'clang-format-' + str(args.llvm),
                    'clang-tidy-' + str(args.llvm)]
        if args.hipsycl is not None:
            packages += [f'llvm-{args.llvm}-dev',
                         f'libclang-{args.llvm}-dev',
                         f'lld-{args.llvm}']
        return packages
    else:
        return []


def get_opencl_packages(args) -> typing.List[str]:
    if (args.doxygen is None) and (args.oneapi is None):
        return _opencl_extra_packages
    else:
        return []


def get_rocm_packages(args) -> typing.List[str]:
    if (args.rocm is None):
        return []
    else:
        return _rocm_extra_packages

def get_cp2k_packages(args) -> typing.List[str]:
    if args.mpi is not None:
        packages = _cp2k_extra_packages + ['libfftw3-mpi-dev']

    if (args.cp2k is None):
        return []
    else:
        return packages

def get_compiler(args, compiler_build_stage: hpccm.Stage = None) -> bb_base:
    # Compiler
    if args.llvm is not None:
        # Build our own version instead to get TSAN + OMP
        if args.tsan is not None:
            if compiler_build_stage is not None:
                compiler = compiler_build_stage.runtime(_from='tsan')
            else:
                raise RuntimeError('No TSAN compiler build stage!')
        # Build the default compiler if we don't need special support
        else:
            # Always use the "upstream" llvm repositories because the
            # main ubuntu repositories stop adding support for new
            # llvm versions after a few llvm releases.
            compiler = hpccm.building_blocks.llvm(version=args.llvm, upstream=True)

    elif args.oneapi is not None:
        if compiler_build_stage is not None:
            compiler = compiler_build_stage.runtime(_from='oneapi')
            # Prepare the toolchain (needed only for builds done within the Dockerfile, e.g.
            # OpenMPI builds, which don't currently work for other reasons)
            oneapi_toolchain = hpccm.toolchain(CC=f'/opt/intel/oneapi/compiler/{args.oneapi}/linux/bin/intel64/icx',
                                               CXX=f'/opt/intel/oneapi/compiler/{args.oneapi}/linux/bin/intel64/icpx')
            setattr(compiler, 'toolchain', oneapi_toolchain)

        else:
            raise RuntimeError('No oneAPI compiler build stage!')

    elif args.gcc is not None:
        if args.cp2k is not None:
            compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gcc,
                                             fortran=True)
        else:
            compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gcc,
                                             fortran=False)
    else:
        raise RuntimeError('Logic error: no compiler toolchain selected.')
    return compiler


def get_gdrcopy(args, compiler):
    if args.cuda is not None:
        if hasattr(compiler, 'toolchain'):
            # Version last updated June 7, 2021
            return hpccm.building_blocks.gdrcopy(toolchain=compiler.toolchain, version="2.2")
        else:
            raise RuntimeError('compiler is not an HPCCM compiler building block!')
    else:
        return None


def get_ucx(args, compiler, gdrcopy):
    if args.cuda is not None:
        if hasattr(compiler, 'toolchain'):
            use_gdrcopy = (gdrcopy is not None)
            # Version last updated June 7, 2021
            return hpccm.building_blocks.ucx(toolchain=compiler.toolchain, gdrcopy=use_gdrcopy, version="1.10.1",
                                             cuda=True)
        else:
            raise RuntimeError('compiler is not an HPCCM compiler building block!')
    else:
        return None


def get_mpi(args, compiler, ucx):
    # If needed, add MPI to the image
    if args.mpi is not None:
        if args.mpi == 'openmpi':
            if hasattr(compiler, 'toolchain'):
                if args.oneapi is not None:
                    raise RuntimeError('oneAPI building OpenMPI is not supported')
                use_cuda = (args.cuda is not None)
                use_ucx = (ucx is not None)
                # Version last updated June 7, 2021
                return hpccm.building_blocks.openmpi(toolchain=compiler.toolchain, version="4.1.1", cuda=use_cuda,
                                                     ucx=use_ucx, infiniband=False)
            else:
                raise RuntimeError('compiler is not an HPCCM compiler building block!')

        elif args.mpi == 'impi':
            # TODO Intel MPI from the oneAPI repo is not working reliably,
            # reasons are unclear. When solved, add packagages called:
            # 'intel-oneapi-mpi', 'intel-oneapi-mpi-devel'
            # during the compiler stage.
            # TODO also consider hpccm's intel_mpi package if that doesn't need
            # a license to run.
            raise RuntimeError('Intel MPI recipe not implemented yet.')
        else:
            raise RuntimeError('Requested unknown MPI implementation.')
    else:
        return None


def get_clfft(args):
    if (args.clfft is not None):
        return hpccm.building_blocks.generic_cmake(
            repository='https://github.com/clMathLibraries/clFFT.git',
            prefix='/usr/local', recursive=True, branch=args.clfft, directory='clFFT/src')
    else:
        return None


def get_heffte(args):
    if (args.heffte is not None):
        return hpccm.building_blocks.generic_cmake(
            cmake_opts=['-D CMAKE_BUILD_TYPE=Release',
                        '-D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda',
                        '-D Heffte_ENABLE_CUDA=ON',
                        '-D Heffte_ENABLE_FFTW=OFF',
                        '-D BUILD_SHARED_LIBS=ON'],
            repository='https://bitbucket.org/icl/heffte.git',
            prefix='/usr/local', recursive=True, commit=args.heffte, directory='heffte')
    else:
        return None

def get_nvhpcsdk(args):
    if (args.nvhpcsdk is not None):
        return hpccm.building_blocks.nvhpc(eula=True, cuda_multi=False, environment=False, mpi=False, version=args.nvhpcsdk)
    else:
        return None

def get_hipsycl(args):
    if args.hipsycl is None:
        return None
    if args.llvm is None:
        raise RuntimeError('Can not build hipSYCL without llvm')

    if args.rocm is None:
        raise RuntimeError('hipSYCL requires the rocm packages')

    cmake_opts = ['-DLLVM_DIR=/opt/rocm/llvm/lib/cmake/llvm',
                  '-DCMAKE_PREFIX_PATH=/opt/rocm/lib/cmake',
                  '-DWITH_ROCM_BACKEND=ON']
    if args.cuda is not None:
        cmake_opts += ['-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda',
                       '-DWITH_CUDA_BACKEND=ON']

    postinstall = [
        # https://github.com/illuhad/hipSYCL/issues/361#issuecomment-718943645
        'for f in /opt/rocm/amdgcn/bitcode/*.bc; do ln -s "$f" "/opt/rocm/lib/$(basename $f .bc).amdgcn.bc"; done'
    ]
    if args.cuda is not None:
        postinstall += [
            # https://github.com/illuhad/hipSYCL/issues/410#issuecomment-743301929
            f'sed s/_OPENMP/__OPENMP_NVPTX__/ -i /usr/lib/llvm-{args.llvm}/lib/clang/*/include/__clang_cuda_complex_builtins.h',
        ]

    return hpccm.building_blocks.generic_cmake(
        repository='https://github.com/illuhad/hipSYCL.git',
        directory='/var/tmp/hipSYCL',
        prefix='/usr/local', recursive=True, commit=args.hipsycl,
        cmake_opts=['-DCMAKE_BUILD_TYPE=Release', *cmake_opts],
        postinstall=postinstall)

def get_cp2k(args):
    if args.cp2k is None:
        return None

    if args.gcc is None:
        raise RuntimeError('CP2K build requires GNU compilers')

    make_commands = ['make ARCH=local VERSION=ssmp libcp2k']
    if args.mpi is not None:
        make_commands += ['make ARCH=local VERSION=psmp libcp2k']
    make_commands += ['rm -rf ./obj']

    return hpccm.building_blocks.generic_build(
        repository='https://github.com/cp2k/cp2k.git',
        branch=f'support/v{args.cp2k}',
        recursive=True,
        build=['mkdir -p /opt/cp2k',
               'cp -rf ./* /opt/cp2k/',
               'cd /opt/cp2k/tools/toolchain',
               './install_cp2k_toolchain.sh  \
                --with-gcc=system            \
                --with-cmake=system          \
                --with-openmpi=system        \
                --with-fftw=system           \
                --with-openblas=system       \
                --with-scalapack=install     \
                --with-gsl=no                \
                --with-elpa=no               \
                --with-spglib=no             \
                --with-spfft=no              \
                --with-cosma=no              \
                --with-libvori=no            \
                --with-sirius=no             \
                --with-hdf5=no               \
                --with-libxc=install         \
                --with-libxsmm=install       \
                --with-libint=install        \
                --libint-lmax=5',
                'rm -rf ./build',
                'cp ./install/arch/local.* ../../arch/',
                'cd ../../']
                 + make_commands)

def add_tsan_compiler_build_stage(input_args, output_stages: typing.Mapping[str, hpccm.Stage]):
    """Isolate the expensive TSAN preparation stage.

    This is a very expensive stage, but has few and disjoint dependencies, and
    its output is easily compartmentalized (/usr/local) so we can isolate this
    build stage to maximize build cache hits and reduce rebuild time, bookkeeping,
    and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError('Need output_stages container.')
    if 'compiler_build' in output_stages:
        raise RuntimeError('"compiler_build" output stage is already present.')
    tsan_stage = hpccm.Stage()
    tsan_stage += hpccm.primitives.baseimage(image=base_image_tag(input_args),
                                             _distro=hpccm_distro_name(input_args),
                                             _as='tsan')

    tsan_stage += hpccm.building_blocks.packages(ospackages=['git', 'ca-certificates', 'build-essential', 'cmake'])
    # CMake will get duplicated later, but this is an expensive image, and it isn't worth optimizing
    # out that duplication...
    tsan_stage += hpccm.building_blocks.python(python3=True, python2=False, devel=False)

    compiler_branch = 'release/' + str(input_args.llvm) + '.x'
    tsan_stage += hpccm.building_blocks.generic_cmake(
        repository='https://github.com/llvm/llvm-project.git',
        directory='/var/tmp/llvm-project/llvm/',
        prefix='/usr/local', recursive=True, branch=compiler_branch,
        cmake_opts=['-D CMAKE_BUILD_TYPE=Release',
                    '-D LLVM_ENABLE_PROJECTS="clang;openmp;clang-tools-extra;compiler-rt;lld"',
                    '-D LIBOMP_TSAN_SUPPORT=on'],
        postinstall=['ln -s /usr/local/bin/clang++ /usr/local/bin/clang++-' + str(input_args.llvm),
                     'ln -s /usr/local/bin/clang-format /usr/local/bin/clang-format-' + str(input_args.llvm),
                     'ln -s /usr/local/bin/clang-tidy /usr/local/bin/clang-tidy-' + str(input_args.llvm),
                     'ln -s /usr/local/share/clang/run-clang-tidy.py /usr/local/bin/run-clang-tidy-'
                     + str(input_args.llvm) + '.py',
                     'ln -s /usr/local/bin/run-clang-tidy-'
                     + str(input_args.llvm) + '.py /usr/local/bin/run-clang-tidy-' + str(input_args.llvm),
                     'ln -s /usr/local/libexec/c++-analyzer /usr/local/bin/c++-analyzer-' + str(input_args.llvm)])
    output_stages['compiler_build'] = tsan_stage


def oneapi_runtime(_from='0'):
    oneapi_runtime_stage = hpccm.Stage()
    oneapi_runtime_stage += hpccm.primitives.copy(_from='oneapi-build',
                                                  files={"/opt/intel": "/opt/intel",
                                                         "/etc/bash.bashrc": "/etc/bash.bashrc"})
    return oneapi_runtime_stage


def add_oneapi_compiler_build_stage(input_args, output_stages: typing.Mapping[str, hpccm.Stage]):
    """Isolate the oneAPI preparation stage.

    This stage is isolated so that its installed components are minimized in the
    final image (chiefly /opt/intel) and its environment setup script can be
    sourced. This also helps with rebuild time and final image size.
    """
    if not isinstance(output_stages, collections.abc.MutableMapping):
        raise RuntimeError('Need output_stages container.')
    if 'compiler_build' in output_stages:
        raise RuntimeError('"compiler_build" output stage is already present.')
    oneapi_stage = hpccm.Stage()
    oneapi_stage += hpccm.primitives.baseimage(image=base_image_tag(input_args),
                                               _distro=hpccm_distro_name(input_args),
                                               _as='oneapi-build')

    version = str(input_args.oneapi)

    # Add required components for the next stage (both for hpccm and Intel's setvars.sh script)
    oneapi_stage += hpccm.building_blocks.packages(ospackages=['wget', 'gnupg2', 'ca-certificates', 'lsb-release'])
    oneapi_stage += hpccm.building_blocks.packages(
        apt_keys=['https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB'],
        apt_repositories=['deb https://apt.repos.intel.com/oneapi all main'],
        # Add minimal packages (not the whole HPC toolkit!)
        ospackages=[f'intel-oneapi-dpcpp-cpp-{version}',
                    f'intel-oneapi-openmp-{version}',
                    f'intel-oneapi-mkl-{version}',
                    f'intel-oneapi-mkl-devel-{version}']
    )
    # Ensure that all bash shells on the final container will have access to oneAPI
    oneapi_stage += hpccm.primitives.shell(
        commands=['echo "source /opt/intel/oneapi/setvars.sh" >> /etc/bash.bashrc',
                  'unlink /opt/intel/oneapi/compiler/latest',
                  f'ln -sf /opt/intel/oneapi/compiler/{version} /opt/intel/oneapi/compiler/latest']
    )
    setattr(oneapi_stage, 'runtime', oneapi_runtime)

    output_stages['compiler_build'] = oneapi_stage


def prepare_venv(version: StrictVersion) -> typing.Sequence[str]:
    """Get shell commands to set up the venv for the requested Python version."""
    major = version.version[0]
    minor = version.version[1]  # type: int

    pyenv = '$HOME/.pyenv/bin/pyenv'

    py_ver = f'{major}.{minor}'
    venv_path = f'$HOME/venv/py{py_ver}'
    commands = [f'$({pyenv} prefix `{pyenv} whence python{py_ver}`)/bin/python -m venv {venv_path}']

    commands.append(f'{venv_path}/bin/python -m pip install --upgrade pip setuptools')
    # Install dependencies for building and testing gmxapi Python package.
    # WARNING: Please keep this list synchronized with python_packaging/src/requirements.txt
    # TODO: Get requirements.txt from an input argument.
    commands.append(f"""{venv_path}/bin/python -m pip install --upgrade \
            'breathe' \
            'cmake>=3.16.3' \
            'flake8>=3.7.7' \
            'gcovr>=4.2' \
            'mpi4py>=3.0.3' \
            'networkx>=2.0' \
            'numpy>1.7' \
            'packaging' \
            'pip>=10.1' \
            'pybind11>2.6' \
            'Pygments>=2.2.0' \
            'pytest>=4.6' \
            'setuptools>=42' \
            'Sphinx>=1.6.3' \
            'sphinxcontrib-plantuml>=0.14' \
            'wheel'""")
    return commands


def add_python_stages(input_args: argparse.Namespace, *,
                      base: str,
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
    pyenv_stage += hpccm.primitives.baseimage(image=base,
                                              _distro=hpccm_distro_name(input_args),
                                              _as='pyenv')
    pyenv_stage += hpccm.building_blocks.packages(ospackages=_python_extra_packages)

    for version in [StrictVersion(py_ver) for py_ver in sorted(input_args.venvs)]:
        stage_name = 'py' + str(version)
        stage = hpccm.Stage()
        stage += hpccm.primitives.baseimage(image=base,
                                            _distro=hpccm_distro_name(input_args),
                                            _as=stage_name)
        stage += hpccm.building_blocks.packages(ospackages=_python_extra_packages)

        # TODO: Use a non-root user for testing and Python virtual environments.
        stage += hpccm.primitives.shell(commands=[
            'curl https://pyenv.run | bash',
            """echo 'export PYENV_ROOT="$HOME/.pyenv"' >> $HOME/.bashrc""",
            """echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> $HOME/.bashrc""",
            """echo 'eval "$(pyenv init -)"' >> $HOME/.bashrc""",
            """echo 'eval "$(pyenv virtualenv-init -)"' >> $HOME/.bashrc"""])
        pyenv = '$HOME/.pyenv/bin/pyenv'
        commands = [f'PYTHON_CONFIGURE_OPTS="--enable-shared" {pyenv} install -s {version}']
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


def add_documentation_dependencies(input_args,
                                   output_stages: typing.MutableMapping[str, hpccm.Stage]):
    """Add appropriate layers according to doxygen input arguments."""
    if input_args.doxygen is None:
        return
    # Always clone the same version of linkchecker (latest release at June 1, 2021)
    output_stages['main'] += hpccm.building_blocks.pip(pip='pip3', packages=[
        'git+https://github.com/linkchecker/linkchecker.git@v10.0.1'])
    output_stages['main'] += hpccm.primitives.shell(
        commands=[
            'sed -i \'/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d\' /etc/ImageMagick-6/policy.xml'])
    if input_args.doxygen == '1.8.5':
        doxygen_commit = 'ed4ed873ab0e7f15116e2052119a6729d4589f7a'
        output_stages['main'] += hpccm.building_blocks.generic_autotools(
            repository='https://github.com/westes/flex.git',
            commit='f7788a9a0ecccdc953ed12043ccb59ca25714018',
            prefix='/tmp/install-of-flex',
            configure_opts=['--disable-shared'],
            preconfigure=['./autogen.sh'])
        output_stages['main'] += hpccm.building_blocks.generic_autotools(
            repository='https://github.com/doxygen/doxygen.git',
            commit=doxygen_commit,
            prefix='',
            configure_opts=[
                '--flex /tmp/install-of-flex/bin/flex',
                '--static'])
    else:
        version = input_args.doxygen
        archive_name = f'doxygen-{version}.linux.bin.tar.gz'
        archive_url = f'https://sourceforge.net/projects/doxygen/files/rel-{version}/{archive_name}'
        binary_path = f'doxygen-{version}/bin/doxygen'
        commands = [
            'mkdir doxygen && cd doxygen',
            f'wget {archive_url}',
            f'tar xf {archive_name} {binary_path}',
            f'cp {binary_path} /usr/local/bin/',
            'cd .. && rm -rf doxygen'
        ]
        output_stages['main'] += hpccm.primitives.shell(commands=commands)


def add_base_stage(name: str,
                   input_args,
                   output_stages: typing.MutableMapping[str, hpccm.Stage]):
    """Establish dependencies that are shared by multiple parallel stages."""
    # Building blocks are chunks of container-builder instructions that can be
    # copied to any build stage with the addition operator.
    building_blocks = collections.OrderedDict()
    building_blocks['base_packages'] = hpccm.building_blocks.packages(
        ospackages=_common_packages)

    # These are the most expensive and most reusable layers, so we put them first.
    building_blocks['compiler'] = get_compiler(input_args, compiler_build_stage=output_stages.get('compiler_build'))
    building_blocks['gdrcopy'] = get_gdrcopy(input_args, building_blocks['compiler'])
    building_blocks['ucx'] = get_ucx(input_args, building_blocks['compiler'], building_blocks['gdrcopy'])
    building_blocks['mpi'] = get_mpi(input_args, building_blocks['compiler'], building_blocks['ucx'])

    # Create the stage from which the targeted image will be tagged.
    output_stages[name] = hpccm.Stage()

    output_stages[name] += hpccm.primitives.baseimage(image=base_image_tag(input_args),
                                                      _distro=hpccm_distro_name(input_args),
                                                      _as=name)
    for bb in building_blocks.values():
        if bb is not None:
            output_stages[name] += bb


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

    # If we need TSAN or oneAPI support the early build is more complex,
    # so that our compiler images don't have all the cruft needed to get those things
    # installed.
    if args.llvm is not None and args.tsan is not None:
        add_tsan_compiler_build_stage(input_args=args, output_stages=stages)
    if args.oneapi is not None:
        add_oneapi_compiler_build_stage(input_args=args, output_stages=stages)

    add_base_stage(name='build_base', input_args=args, output_stages=stages)

    # Add Python environments to MPI images, only, so we don't have to worry
    # about whether to install mpi4py.
    if args.mpi is not None and len(args.venvs) > 0:
        add_python_stages(base='build_base', input_args=args, output_stages=stages)

    # Building blocks are chunks of container-builder instructions that can be
    # copied to any build stage with the addition operator.
    building_blocks = collections.OrderedDict()

    for i, cmake in enumerate(args.cmake):
        building_blocks['cmake' + str(i)] = hpccm.building_blocks.cmake(
            eula=True,
            prefix=f'/usr/local/cmake-{cmake}',
            version=cmake)

    # Install additional packages early in the build to optimize Docker build layer cache.
    os_packages = list(get_llvm_packages(args)) + get_opencl_packages(args) + get_rocm_packages(args) + get_cp2k_packages(args)
    if args.doxygen is not None:
        os_packages += _docs_extra_packages
    if args.oneapi is not None:
        os_packages += ['lsb-release']
    if args.hipsycl is not None:
        os_packages += ['libboost-fiber-dev']
    building_blocks['extra_packages'] = []
    if args.intel_compute_runtime:
        building_blocks['extra_packages'] += hpccm.building_blocks.packages(
            apt_keys=['https://repositories.intel.com/graphics/intel-graphics.key'],
            apt_repositories=[f'deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main']
        )
        os_packages += _intel_compute_runtime_extra_packages
    if args.rocm is not None:
        building_blocks['extra_packages'] += hpccm.building_blocks.packages(
            apt_keys=['http://repo.radeon.com/rocm/rocm.gpg.key'],
            apt_repositories=[f'deb [arch=amd64] http://repo.radeon.com/rocm/apt/{args.rocm}/ ubuntu main']
        )
    building_blocks['extra_packages'] += hpccm.building_blocks.packages(
        ospackages=os_packages,
        apt_ppas=['ppa:intel-opencl/intel-opencl'])

    building_blocks['CP2K'] = get_cp2k(args)

    if args.cuda is not None and args.llvm is not None:
        # Hack to tell clang what version of CUDA we're using
        # based on https://github.com/llvm/llvm-project/blob/1fdec59bffc11ae37eb51a1b9869f0696bfd5312/clang/lib/Driver/ToolChains/Cuda.cpp#L43
        cuda_version_split = args.cuda.split('.')
        # LLVM requires having the version in x.y.z format, while args.cuda be be either x.y or x.y.z
        cuda_version_str = '{}.{}.{}'.format(
            cuda_version_split[0],
            cuda_version_split[1],
            cuda_version_split[2] if len(cuda_version_split) > 2 else 0
        )
        building_blocks['cuda-clang-workaround'] = hpccm.primitives.shell(commands=[
            f'echo "CUDA Version {cuda_version_str}" > /usr/local/cuda/version.txt'
        ])

    building_blocks['clfft'] = get_clfft(args)

    building_blocks['heffte'] = get_heffte(args)

    building_blocks['nvhpcsdk'] = get_nvhpcsdk(args)
    if building_blocks['nvhpcsdk'] is not None:
          nvshmem_lib_path = '/opt/nvidia/hpc_sdk/Linux_x86_64/' + args.nvhpcsdk + '/comm_libs/nvshmem/lib/:$LD_LIBRARY_PATH'
          building_blocks['nvhpcsdk_path'] = hpccm.primitives.environment(variables={'LD_LIBRARY_PATH': nvshmem_lib_path})

    building_blocks['hipSYCL'] = get_hipsycl(args)

    # Add Python environments to MPI images, only, so we don't have to worry
    # about whether to install mpi4py.
    if args.mpi is not None and len(args.venvs) > 0:
        add_python_stages(base='build_base', input_args=args, output_stages=stages)

    # Create the stage from which the targeted image will be tagged.
    stages['main'] = hpccm.Stage()

    stages['main'] += hpccm.primitives.baseimage(image='build_base',
                                                 _distro=hpccm_distro_name(args),
                                                 _as='main')
    for bb in building_blocks.values():
        if bb is not None:
            stages['main'] += bb

    # We always add Python3 and Pip
    stages['main'] += hpccm.building_blocks.python(python3=True, python2=False)

    # Add documentation requirements (doxygen and sphinx + misc).
    if args.doxygen is not None:
        add_documentation_dependencies(args, stages)

    if 'pyenv' in stages and stages['pyenv'] is not None:
        stages['main'] += hpccm.primitives.copy(_from='pyenv', _mkdir=True, src=['/root/.pyenv/'],
                                                dest='/root/.pyenv')
        stages['main'] += hpccm.primitives.copy(_from='pyenv', _mkdir=True, src=['/root/venv/'],
                                                dest='/root/venv')
        # TODO: Update user home directory.
        # TODO: If we activate pyenv for login shells, the `global` "version" should be full-featured.
        # stages['main'] += hpccm.primitives.copy(_from='pyenv', src=['/root/.bashrc'],
        #                                         dest='/root/')

    # Make sure that `python` resolves to something.
    stages['main'] += hpccm.primitives.shell(commands=['test -x /usr/bin/python || '
                                                       'update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && '
                                                       '/usr/bin/python --version'])

    # Note that the list of stages should be sorted in dependency order.
    for build_stage in stages.values():
        if build_stage is not None:
            yield build_stage


if __name__ == '__main__':
    args = parser.parse_args()

    # Set container specification output format
    hpccm.config.set_container_format(args.format)
    # Normally in hpccm the first call to baseimage sets the context
    # for other packages, e.g. which distro version to use. We want to
    # set that early on, so that hpccm can do the right thing to use
    # the right llvm apt repo when we want to use versions of llvm
    # that were never supported by the main apt repo.
    hpccm.config.set_linux_distro(hpccm_distro_name(args))

    container_recipe = build_stages(args)

    # Output container specification
    for stage in container_recipe:
        print(stage)
