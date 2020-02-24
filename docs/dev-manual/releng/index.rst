.. This file is a placeholder that is used in case RELENG_PATH does not
   identify the location of the releng repository.

.. _releng-workflow-release:

Release engineering with Gitlab
===============================

.. toctree::
   :hidden:

We are currently switching our build and testing system to use Gitlab
CI pipelines run on GitLab Runner. This section is going to be extended
with individual build information as it comes available. For now we are
using a combination of building with the previous system on Jenkins
and post-submit verification on Gitlab.

.. seealso:: :doc:`../infrastructure`

.. _releng-triggering-builds:

Triggering builds on Gitlab
---------------------------

Pipelines can be triggered through the web interface, with different
pipelines available through the use of specified environment variables
in the trigger interface.

This section is going to be extended with information for how to trigger
different builds and their individual behaviour.
