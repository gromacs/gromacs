.. This file is a placeholder that is used in case RELENG_PATH does not
   identify the location of the releng repository.

.. _releng-workflow-release:

Release engineering with Gitlab
===============================

.. toctree::
   :hidden:

We are currently switching our build and testing system to use Gitlab
and the integrated CI system, with information for the general system found
in the official `Gitlab documentation <https://docs.gitlab.com/ee/ci/yaml/>`_.
The new configuration for the builds and tests can be found in the file
``.gitlab-ci.yml``, with the templates for configuring is found in the files in the
``admin/ci-templates/`` directory. This section is going to be extended
with individual build information as it comes available. For now we are
using a combination of building with the previous system on Jenkins
and post-submit verification on Gitlab.

.. _releng-triggering-builds:

Triggering builds on Gitlab
---------------------------

Pipelines can be triggered through the web interface, with different
pipelines available through the use of specified environment variables
in the trigger interface.

This section is going to be extended with information for how to trigger
different builds and their individual behaviour.
