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

Refer to :file:`buildall.sh` in the ``master`` branch for the set of images
currently being built.

Utilities
=========

:file:`utility.py`
------------------

.. automodule:: utility
    :members:

:file:`scripted_gmx_docker_builds.py`
-------------------------------------

.. automodule:: scripted_gmx_docker_builds
