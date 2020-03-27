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
                    'ccache',
                    'git',
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

# Parse command line arguments
parser = argparse.ArgumentParser(description='GROMACS CI image creation script', parents=[utility.parser])

parser.add_argument('--format', type=str, default='docker',
                    choices=['docker', 'singularity'],
                    help='Container specification format (default: docker)')


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


def get_compiler(args) -> bb_base:
    # Compiler
    if args.icc is not None:
        raise RuntimeError('Intel compiler toolchain recipe not implemented yet')

    if args.llvm is not None:
        # Build the default compiler if we don't need special support
        if args.tsan is None:
            compiler = hpccm.building_blocks.llvm(extra_repository=True, version=args.llvm)
        # Build our own version instead to get TSAN + OMP
        else:
            compiler_branch = 'release_' + str(args.llvm) + '0'
            compiler = hpccm.building_blocks.generic_cmake(
                repository='https://git.llvm.org/git/llvm.git',
                prefix='/usr/local', recursive=True, branch=compiler_branch,
                cmake_opts=['-D CMAKE_BUILD_TYPE=Release', '-D LLVM_ENABLE_PROJECTS="clang;openmp;clang-tools-extra"',
                            '-D LIBOMP_TSAN_SUPPORT=on'],
                preconfigure=['export branch=' + compiler_branch,
                              '(cd projects; git clone https://git.llvm.org/git/libcxx.git; cd libcxx; git checkout $branch)',
                              '(cd projects; git clone https://git.llvm.org/git/libcxxabi.git; cd libcxxabi; git checkout $branch)',
                              '(cd projects; git clone https://git.llvm.org/git/compiler-rt.git; cd compiler-rt; git checkout $branch)',
                              '(cd ..; git clone https://git.llvm.org/git/openmp.git; cd openmp; git checkout $branch)',
                              '(cd ..; git clone https://git.llvm.org/git/clang.git; cd clang; git checkout $branch)',
                              '(cd ..; git clone https://git.llvm.org/git/clang-tools-extra.git clang-tools-extra; cd clang-tools-extra; git checkout $branch)'],
                postinstall=['ln -s /usr/local/bin/clang++ /usr/local/bin/clang++-' + str(args.llvm),
                             'ln -s /usr/local/bin/clang-format /usr/local/bin/clang-format-' + str(args.llvm),
                             'ln -s /usr/local/bin/clang-tidy /usr/local/bin/clang-tidy-' + str(args.llvm),
                             'ln -s /usr/local/libexec/c++-analyzer /usr/local/bin/c++-analyzer-' + str(args.llvm)])


    elif (args.gnu is not None):
        compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gnu,
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

            return hpccm.building_blocks.openmpi(toolchain=compiler.toolchain, cuda=use_cuda, infiniband=False)
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
            return hpccm.building_blocks.packages(ospackages=['ocl-icd-opencl-dev', 'opencl-headers',
                                                              'beignet-opencl-icd'])
        elif args.opencl == 'amd':
            # Due to the wisdom of AMD, this needs to be done differently for the OS and version! Hurray!
            # And they don't allow wget, so this branch is not taken for now! AMD, please allow me to use wget.
            raise RuntimeError('AMD recipe can not be generated because they do not allow wget for getting the packages.')
            # if args.ubuntu:
            #     if args.ubuntu is not '16.04':
            #         Stage0 += hpccm.building_blocks.generic_build(url='https://www2.ati.com/drivers/linux/ubuntu/'+args.ubuntu+'/amdgpu-pro-18.30-641594.tar.xz',
            #                                                       install=['./amdgpu-install --opencl=legacy --headless -y'])
            #     elif:
            #         Stage0 += hpccm.building_blocks.generic_build(url='https://www2.ati.com/drivers/linux/ubuntu/amdgpu-pro-18.30-641594.tar.xz',
            #                                                       install=['./amdgpu-install --opencl=legacy --headless -y'])
            # elif args.centos:
            #         Stage0 += hpccm.building_blocks.generic_build(url='https://www2.ati.com/drivers/linux/rhel'+args.centos'/amdgpu-pro-18.30-641594.tar.xz',
            #                                                       install=['./amdgpu-install --opencl=legacy --headless -y'])
    else:
        return None


def get_clfft(args):
    if (args.clfft is not None):
        return hpccm.building_blocks.generic_cmake(
            repository='https://github.com/clMathLibraries/clFFT.git',
            prefix='/usr/local', recursive=True, branch=args.clfft, directory='clFFT/src')
    else:
        return None


def get_llvm_packages(args) -> typing.Iterable[str]:
    # If we use the package version of LLVM, we need to install extra packages for it.
    if (args.llvm is not None) and (args.tsan is None):
        return ['libomp-dev',
                'clang-format-' + str(args.llvm),
                'clang-tidy-' + str(args.llvm)]
    else:
        return []


def build_stages(args) -> typing.Iterable[hpccm.Stage]:
    building_blocks = collections.OrderedDict()

    building_blocks['cmake'] = hpccm.building_blocks.cmake(eula=True, version=args.cmake)
    building_blocks['compiler'] = get_compiler(args)
    building_blocks['mpi'] = get_mpi(args, building_blocks['compiler'])
    building_blocks['opencl'] = get_opencl(args)
    building_blocks['clfft'] = get_clfft(args)

    # Create Stage
    main_stage = hpccm.Stage()

    main_stage += hpccm.primitives.baseimage(image=base_image_tag(args))

    # Install additional packages early in the build to optimize Docker build layer cache.
    os_packages = _common_packages + get_llvm_packages(args)
    if args.doxygen is not None:
        os_packages += _docs_extra_packages
    main_stage += hpccm.building_blocks.packages(ospackages=os_packages)

    # We always add Python3 and Pip
    main_stage += hpccm.building_blocks.python(python3=True, python2=False, devel=True)
    main_stage += hpccm.building_blocks.pip(upgrade=True, pip='pip3',
                                            packages=['pytest', 'networkx', 'numpy'])
    for bb in building_blocks.values():
        if bb is not None:
            main_stage += bb

    # Add documentation requirements (doxygen and sphinx + misc).
    if (args.doxygen is not None):
        if (args.doxygen == '1.8.5'):
            doxygen_commit = 'ed4ed873ab0e7f15116e2052119a6729d4589f7a'
        else:
            doxygen_commit = 'a6d4f4df45febe588c38de37641513fd576b998f'
        main_stage += hpccm.building_blocks.generic_autotools(repository='https://github.com/westes/flex.git',
                                                              commit='f7788a9a0ecccdc953ed12043ccb59ca25714018',
                                                              prefix='/tmp/install-of-flex',
                                                              configure_opts=['--disable-shared'],
                                                              preconfigure=['./autogen.sh'])
        main_stage += hpccm.building_blocks.generic_autotools(repository='https://github.com/doxygen/doxygen.git',
                                                              commit=doxygen_commit,
                                                              prefix='',
                                                              configure_opts=['--flex /tmp/install-of-flex/bin/flex',
                                                                              '--static'],
                                                              postinstall=[
                                                                  'sed -i \'/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d\' /etc/ImageMagick-6/policy.xml'])
        main_stage += hpccm.building_blocks.pip(pip='pip3', packages=['sphinx==1.6.1'])

    for build_stage in [main_stage]:
        if build_stage is not None:
            yield build_stage


if __name__ == '__main__':
    args = parser.parse_args()

    container_recipe = build_stages(args)

    # Set container specification output format
    hpccm.config.set_container_format(args.format)

    # Output container specification
    for stage in container_recipe:
        print(stage)
