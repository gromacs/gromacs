# Documentation strategy

## Documentation types

User documentation: Mostly extracted from the most user-facing module docstrings,
but curated in the framework of RST docs to facilitate organization of classes and
functions into sections, potentially mixing auto-extracted documentation with sample
documentation for proposed or incomplete features (to facilitate collaborative design
and code review).

User reference documentation: Function and class docstrings from the most user-facing public entities.

API reference: function and class documentation for public entities not included by
default in the top-level namespace. Should link to user docs when higher-level collaborators exist.

Design documentation and additional developer docs: module documentation for modules
not included in regular documentation.

Documentation for maintainers: static RST for versioning, releases, packaging, testing infrastructure, source code maintenance.

Note that the RST sources are also used by the GROMACS webpage-sphinx CMake target.
