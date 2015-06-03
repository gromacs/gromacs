Command-line reference
======================

.. toctree::
   :hidden:
   :glob:

   /onlinehelp/gmx
   /onlinehelp/gmx-*

|Gromacs| includes many tools for preparing, running and analysing
molecular dynamics simulations. These are all structured as part of a single
:command:`gmx` wrapper binary, and invoked with commands like :command:`gmx grompp`.
:ref:`mdrun <gmx mdrun>` is the only other binary that
:ref:`can be built <building just the mdrun binary>`; in the normal
build it can be run with :command:`gmx mdrun`. Documentation for these can
be found at the respective sections below, as well as on man pages (e.g.,
:manpage:`gmx-grompp(1)`) and with :samp:`gmx help {command}` or
:samp:`gmx {command} -h`.

Command-line interface and conventions
--------------------------------------

All |Gromacs| commands require an option before any arguments (i.e., all
command-line arguments need to be preceded by an argument starting with a
dash, and values not starting with a dash are arguments to the preceding
option).  Most options, except for boolean flags, expect an argument (or
multiple in some cases) after the option name.
The argument must be a separate command-line argument, i.e., separated by
space, as in ``-f traj.xtc``.  If more than one argument needs to be given to
an option, they should be similarly separated from each other.
Some options also have default arguments, i.e., just specifying the option
without any argument uses the default argument.
If an option is not specified at all, a default value is used; in the case of
optional files, the default might be not to use that file (see below).

All |Gromacs| options start with a single dash, whether they are single- or
multiple-letter options.  However, two dashes are also recognized (starting
from 5.1).

In addition to command-specific options, some options are handled by the
:command:`gmx`, and can be specified for any command.  See
:doc:`wrapper binary help </onlinehelp/gmx>` for the list of such options.
These options are recognized both before the command name (e.g.,
:command:`gmx -quiet grompp`) as well as after the command name (e.g.,
:command:`gmx grompp -quiet`).
There is also a ``-hidden`` option that can be specified in combination with
``-h`` to show help for advanced/developer-targeted options.

Most analysis tools can process a trajectory with fewer atoms than the
run input or structure file, but only if the trajectory consists of the
first *n* atoms of the run input or structure file.

Handling specific option types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

boolean options
  Boolean flags can be specified like ``-pbc`` and negated like ``-nopbc``.
  It is also possible to use an explicit value like ``-pbc no`` and
  ``-pbc yes``.
file options
  Options that accept files have several features to support default file
  names:

  * If the option/file is marked optional, it is not used unless the option
    is set (or other conditions make the file required).
  * If an option is set, but no value is provided, a default file name is
    always used.
  * If a required option is not set, a default file name is used.

  All such options will accept file names without a file extension.
  The extension is automatically appended in such a case.
  When multiple input formats are accepted, such as a generic structure format,
  the directory will be searched for files of each type with the supplied or
  default name. When no file with a recognized extension is found, an error is given.
  For output files with multiple formats, a default file type will be used.

  Some file formats can also be read from compressed (:file:`.Z` or
  :file:`.gz`) formats.
enum options
  Many commands accept enumerated options (enum), which should be used with one
  of the arguments listed in the option description. The argument may be
  abbreviated, and the first match to the shortest argument in the list will be
  selected.
vector options
  Many commands also use options that may accept a vector of values.
  Either 1 or 3 parameters can be supplied; when only one parameter is supplied
  the two other values are also set to this value.
selection options
  See :doc:`/onlinehelp/selections`.

Commands by name
----------------

.. include:: /fragments/byname.rst

Commands by topic
-----------------

.. include:: /fragments/bytopic.rst

Special topics
--------------

The information in these topics is also accessible through
:samp:`gmx help {topic}` on the command line.

.. toctree::

   /onlinehelp/selections
