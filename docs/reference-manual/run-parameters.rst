Run parameters and Programs
===========================

Online documentation
--------------------

We install standard UNIX man pages for all the programs. If
you have sourced the ``GMXRC`` script in the |Gromacs| binary directory for
your host they should already be present in your ``MANPATH`` environment
variable, and you should be able to type *e.g.* ``man gmx-grompp``. You can
also use the ``-h`` flag on the command line (e.g. :ref:`gmx grompp <gmx grompp>` ``-h``) to see the
same information, as well as ``gmx help grompp``. The list of all programs
are available from :ref:`gmx help <gmx help>`.

File types
----------

Information about different file types can be found
in :doc:`file-formats`.

|Gromacs| files written in XDR format can be read on any architecture with
|Gromacs| version 1.6 or later if the configuration script found the XDR
libraries on your system. They should always be present on UNIX since
they are necessary for NFS support.

Run Parameters
--------------

The descriptions of :ref:`mdp` parameters can be found at
under the link above both in your local |Gromacs| installation,
or :ref:`here <mdp-general>`.

.. raw:: latex

    \clearpage



