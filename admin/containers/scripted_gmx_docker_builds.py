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

import hpccm
import hpccm.config

try:
    import utility
except ImportError:
    raise RuntimeError(
        'This module assumes availability of supporting modules in the same directory. Add the directory to '
        'PYTHONPATH or invoke Python from within the module directory so module location can be resolved.')


# Parse command line arguments
parser = argparse.ArgumentParser(description='GROMACS CI image creation script', parents=[utility.parser])

parser.add_argument('--format', type=str, default='docker',
                    choices=['docker', 'singularity'],
                    help='Container specification format (default: docker)')

def main(args) -> hpccm.Stage:
    # Create Stage
    Stage0 = hpccm.Stage()

    # Create string for base image tag
    base_image_tag = str()

    # Check if we use CUDA images or plain linux images
    if (args.cuda is not None):
        cuda_version_tag = 'nvidia/cuda:' + args.cuda + '-devel'
        if (args.centos is not None):
            cuda_version_tag += '-centos' + args.centos
        elif (args.ubuntu is not None):
            if ((args.cuda == '9.0') and (args.ubuntu == '18.04')):
                raise RuntimeError('Can not combine CUDA 9.0 and Ubuntu 18.04')
            cuda_version_tag += '-ubuntu' + args.ubuntu
        else:
            raise RuntimeError('Logic error: no Linux distribution selected.')

        base_image_tag = cuda_version_tag
    else:
        if (args.centos is not None):
            base_image_tag = 'centos:centos' + args.centos
        elif (args.ubuntu is not None):
            base_image_tag = 'ubuntu:' + args.ubuntu
        else:
            raise RuntimeError('Logic error: no Linux distribution selected.')

    Stage0 += hpccm.primitives.baseimage(image=base_image_tag)

    # Install the GROMACS packages we always will need for our builds.
    Stage0 += hpccm.building_blocks.packages(ospackages=['build-essential',
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
                                                         'wget',
                                                         'xsltproc'])

    # Add CMake to image
    Stage0 += hpccm.building_blocks.cmake(eula=True, version=args.cmake)

    # We always add Python3 and Pip
    Stage0 += hpccm.building_blocks.python(python3=True, python2=False)
    Stage0 += hpccm.building_blocks.pip(upgrade=True, pip='pip3',
                                        packages=['pytest', 'networkx', 'numpy'])

    # Compiler
    if (args.icc is not None):
        raise RuntimeError('Intel compiler toolchain recipe not implemented yet')

    if (args.llvm is not None):
        # Build the default compiler if we don't need special support
        if (args.tsan is None):
            if (args.llvm == 3):
                if ((args.ubuntu is not None) and (args.ubuntu == '18.04')):
                    raise RuntimeError('LLVM 3 and Ubuntu 18.04 can cause issues when used together')
                args.llvm = 3.6
            compiler = hpccm.building_blocks.llvm(extra_repository=True, version=args.llvm)
        # Build our own version instead to get TSAN + OMP
        else:
            compiler_branch = 'release_'+str(args.llvm)+'0'
            compiler = hpccm.building_blocks.generic_cmake(repository='https://git.llvm.org/git/llvm.git',
                    prefix='/usr/local', recursive=True, branch=compiler_branch,
                    cmake_opts=['-D CMAKE_BUILD_TYPE=Release', '-D LLVM_ENABLE_PROJECTS="clang;openmp"', '-D LIBOMP_TSAN_SUPPORT=on'],
                    preconfigure=['export branch='+compiler_branch,
                                  '(cd projects; git clone https://git.llvm.org/git/libcxx.git; cd libcxx; git checkout $branch)',
                                  '(cd projects; git clone https://git.llvm.org/git/libcxxabi.git; cd libcxxabi; git checkout $branch)',
                                  '(cd projects; git clone https://git.llvm.org/git/compiler-rt.git; cd compiler-rt; git checkout $branch)',
                                  '(cd ..; git clone https://git.llvm.org/git/openmp.git; cd openmp; git checkout $branch)',
                                  '(cd ..; git clone https://git.llvm.org/git/clang.git; cd clang; git checkout $branch)',
                                  '(cd ../clang/tools; git clone https://git.llvm.org/git/clang-tools-extra.git extra; cd extra; git checkout $branch)'],
                    postinstall=['ln -s /usr/local/bin/clang++ /usr/local/bin/clang++-'+str(args.llvm),
                                 'ln -s /usr/local/bin/clang-format /usr/local/bin/clang-format-'+str(args.llvm),
                                 'ln -s /usr/local/bin/clang-tidy /usr/local/bin/clang-tidy-'+str(args.llvm)])


    elif (args.gnu is not None):
        compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gnu,
                                             fortran=False)
    else:
        raise RuntimeError('Logic error: no compiler toolchain selected.')


    Stage0 += compiler
    # If we use the package version of LLVM, we need to install extra packages for it.
    if (args.llvm is not None) and (args.tsan is None):
        Stage0 += hpccm.building_blocks.packages(ospackages=['libomp-dev',
                                                             'clang-format-'+str(args.llvm),
                                                             'clang-tidy-'+str(args.llvm)])

    # If needed, add MPI to the image
    if (args.mpi is not None):
        if args.mpi == 'openmpi':
            use_cuda = False
            if (args.cuda is not None):
                use_cuda = True

            Stage0 += hpccm.building_blocks.openmpi(toolchain=compiler.toolchain, cuda=use_cuda, infiniband=False)
        elif args.mpi == 'impi':
            raise RuntimeError('Intel MPI recipe not implemented yet.')

    # Add OpenCL environment if needed
    if (args.opencl is not None):
        if args.opencl == 'nvidia':
            if (args.cuda is None):
                raise RuntimeError('Need Nvidia environment for Nvidia OpenCL image')

            Stage0 += hpccm.building_blocks.packages(ospackages=['nvidia-opencl-dev'])

        elif args.opencl == 'intel':
            Stage0 += hpccm.building_blocks.packages(ospackages=['ocl-icd-opencl-dev', 'opencl-headers', 'beignet-opencl-icd'])
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

        if (args.clfft is not None):
            Stage0 += hpccm.building_blocks.generic_cmake(repository='https://github.com/clMathLibraries/clFFT.git',
                                                          prefix='/usr/local', recursive=True, branch=args.clfft, directory='clFFT/src')


    return Stage0


if __name__ == '__main__':
    args = parser.parse_args()

    container_recipe = main(args)

    # Set container specification output format
    hpccm.config.set_container_format(args.format)

    # Output container specification
    print(container_recipe)
