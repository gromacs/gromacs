==========
Containers
==========

|Gromacs| project infrastructure uses Docker containerization to
isolate automated tasks.
A number of images are maintained to provide a breadth of testing coverage.

Scripts and configuration files for building images are stored in the repository
under :file:`admin/containers/`
Images are (re)built manually by |Gromacs| project staff and pushed to
repositories at https://hub.docker.com/u/gromacs

Utilities
=========

.. automodule:: utility
    :members:

HPC container maker
-------------------

We use the `NVidia HPC Container Maker <https://github.com/NVIDIA/hpc-container-maker>`__
package for scripted Dockerfile generation.
See :file:`admin/containers/scripted_gmx_docker_builds.py`.

.. automodule:: scripted_gmx_docker_builds
