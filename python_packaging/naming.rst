Notes on package names and namespaces
=====================================

CMake
-----

* Project name at root of GROMACS repository is ``Gromacs``. Library target is ``libgromacs``.
* "gmxapi" library target in GROMACS repository is importable as ``Gromacs::gmxapi``.
* Kasson Lab "gmxapi" Python package for 2019 has CMake project name ``gmxpy``.
* The example repository with usable model MD restraint code has a project name matching
  the repository name: ``sample_restraint``. Users are encouraged to change the name for their own
  code, but probably don't.

C++ namespaces
--------------

``gmx``
~~~~~~~

The GROMACS C++ code nominally uses the ``::gmx`` namespace, with nested namespaces used for details
that would not be exposed in the public library object.

``gmxapi``
~~~~~~~~~~

The 2019 Kasson Lab contributions attempt to establish a demarcation such that
* the ``libgmxapi`` library defines only symbols in ``::gmxapi``,
* all symbols and types in the top level of ``::gmxapi`` are documented in the public API documentation,
* symbols and types accessible through the installed headers and library are documented and well-defined
  unless explicitly noted as implementation details,
* non-fully-defined types exposed in signatures are in clear foward-declaration headers that
  indicate where and under what circumstances the types _are_ fully defined,

With the consolidation of the installed libraries and the deprecation of the public / installed API,
it is unclear what the policies will be regarding naming, namespace, documentation, symbol visibility,
API level, and any conventions for the mutual implications of these details.

``gmxpy``
~~~~~~~~~

Used internally in the Python bindings module supporting the ``gmx`` Python package. No installed headers.

``plugin``
~~~~~~~~~~

Used internally in the ``sample_restraint`` code. Users are encouraged to change the name for their own
code, but probably don't.

Python packages
---------------

Note the distinction between ``import``able *import packages* and (installable) *distribution packages*
(which can provide multiple import packages).

A Python distribution package is usually automatically named by a packaging tool from the
``name`` argument to the ``setup`` command in the ``setup.py`` file.
This name conventionally matches the primary package installed with the distribution,
but is not required to and often does not, especially when multiple packages are installed
from the same ``setup.py`` file.

It is generally convenient if the primary package name, distribution name, and repository
name are the same. It is generally inconvenient for a package name to be different
than the name of the subdirectory providing the package source. The installed filesystem
name of a module or package is generally required to match the importable module name.

As of 2019, the Kasson Lab Python package has a base module named ``gmx`` (from https://github.com/kassonlab/gmxapi.git),
anticipating that functionality would resemble that of the ``gmx`` command line tool and that the ``gmx``
Python module might ultimately serve as a replacement for the CLI tool when invoking the module's
``__main__`` method. Alternatively, though, the Python package could provide an "entry point" named "gmx"
through Python package management metadata, regardless of the package or distribution name.

C++ Python extension modules
----------------------------

Python bindings to C or C++ code are provided for runtime "import" with dynamic linking
shared object libraries. When accompanied by a pure Python component, the extension module
is commonly imported as a submodule of the pure Python package, where the extension module
is named the same as the package, but prefixed with an underscore (``\_``).

The ``gmx`` package uses a C++ extension in a submodule named ``core``. The expectation was
that the "gmx" name would imply user-level functionality familiar to CLI users, while the
bindings extension module would map closely to the new public C++ API. But there is also
the case for a compiled module for acceleration or support of the ``gmx`` package that does
not provide a stable API for external use. ``gmx.core`` could plausibly be split into a module
named to match the C++ public interface and another with a subscripted name to match the
Python user interface.

PyPI distributions
------------------

A PyPI distribution archive uses the distribution name defined in the ``setup.py`` file.

My ``uvaerici`` PyPI account owns projects with stub distributions for
``gmx`` (linked to the gmxapi.readthedocs.org page, which is linked to https://github.com/gmxapi/gmxapi),
``gromacs`` (linked to https://github.com/kassonlab/gromacs-gmxapi),
and ``gmxapi`` (linked to https://github.com/kassonlab/gmxapi).
We should think about what to do with those names, what they should mean, and who should own them.
Note, though, that we should control project names corresponding to the
first component of any entry point group names that we adopt (see below).

PyPI projects can have multiple "collaborators" with varying degrees of access.
Right now, my ``eirrgang`` PyPI account is also an owner of these projects,
so I can use ``twine`` to publish updates to the distributions from various credentials
and can extend that authority to other people or automated tools.

Plugins and entry points
------------------------

We have avoided imposing a burden on gmxapi extension authors to
formally package their code, instead leaving it to the user to make sure
that they have prepared the environment for a Python interpreter to be able
to ``import`` additional code. We could provide more formal Python plugin
support as an alternative to extension authors or as a route to implementing
some of the core interoperability in a unified namespace without hard
packaging dependencies.

The Python Packaging Authority specifies that entry point group names
"should use names starting with a PyPI name owned by the consumer project,"
(https://packaging.python.org/specifications/entry-points/#data-model)
so we currently have the flexibility to use ``pkg_resources`` for registering
different providers of functionality, such as through ``setup.py`` files, with kwargs such as
* ``entry_points={'gmx.tool': 'rmsf = somemodule:some_rmsf_function'}``
* ``entry_points={'gmxapi': 'mdrun = gmxapi._builtin_operations.mdrun'}``
* ``entry_points={'gromacs.plugins': 'myspecialmethod = myplugin:myspecialmethod'}``
* ``entry_points={'gmxapi.context': 'mpi = gromacs._gromacs:MpiContext'}``

Note (from https://packaging.python.org/specifications/entry-points/):
"If different distributions provide the same name, the consumer decides how to handle such conflicts."
In other words, there may be multiple equivalent string values ``entry_point.name`` in
``for entry_point in pkg_resources.iter_entry_points('gmxapi.context')``

Python naming proposal (partial)
--------------------------------

Summary:
* GROMACS bundles a Python package installed and imported as "gmxapi".
* For brevity and consistency, we update our Python convention to ``import gmxapi as gmx``.
* We can extend the ``gmxapi`` package bundled with GROMACS through the Python packaging system,
  if we find that experimental functionality develops so much conflict with the package in the
  master branch that separate packages need to be maintained.
* We can preemptively name our experimental package "gmxapi_kl"

Reserve the distribution name "gromacs" for the core library such that, in the future,
``pip install gromacs`` or ``conda install gromacs`` explicitly refers to a
full GROMACS installation with public API enabled.
(or some sort of stub or metapackage referring to the GROMACS installation).
We could also use this as an unadvertised package to help us test experimental
versions of the Python package against commits in or targeted for the GROMACS
``master`` branch.

Use the distribution name "gmxapi" for the Python distribution package.

Avoid "gmx" in Python packaging unless/until
* the distinction is better defined between the ``gmx`` C++ name space,
  the ``gmx`` command line tool name space, and the ``gmx`` Python name space
* Python package compatibility policies are not confusing with respect to those of anything else named ``gmx``

As long as no module or package is actually named "gmx," I think it is okay to
use "gmx" as an abbreviation for "GROMACS" or "gmxapi", and to document conventions
such as ``import gmxapi as gmx``

For gmxapi "operations" or features in the ``gromacs`` or ``gmxapi`` gmxapi operation name spaces,
use ``pkg_resources`` entry points to allow us to provide pre-release or experimental
operation implementations without naming conflicts or confusing naming conventions.

Allow GROMACS 2019 and the gmxapi 0.0.7 release to remain synonymous with the
distribution package and import package named ``gmx``.
I can make a source distribution for PyPI if it makes sense to fill that space.