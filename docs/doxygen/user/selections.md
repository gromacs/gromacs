Dynamic selections {#page_selections}
==================

The \ref module_selection module provides a mechanism that allows selections
specified as text, and the engine evaluates them to atoms, or more generally to
a set of positions, for one or more sets of coordinates.  The selected atoms
can depend on the trajectory frame.  This allows writing general-purpose
analysis tools that only operate on positions, and get a lot of flexibility for
free from the selection engine.  For example, such tools can readily operate on
centers of mass of groups in addition to individual atoms as long as they do
not require access to atomic properties.

For people familiar with VMD, the selection syntax is quite familiar, but there
are some differences.  Not all the keywords supported by VMD are there, and
there are some extensions related to the support to evaluate to center-of-mass
positions in addition to individual atoms.
For old-time \Gromacs users, tools that support selections do not generally
need `make_ndx`.

Structural overview
===================

Central concepts useful for understanding the selection engine are explained
below.  A graph represents the relations between the different parts, and a
textual description of the user-visible components and other concepts follows.
The graph also includes an overview of how the selection engine integrates into
the \ref page_analysisframework.  When using selections from the analysis
framework, the parts in gray are managed by the framework.
When using selections outside the framework, it is either possible to use only
the core components (shown in the graph as a box), or to also use the selection
option mechanisms.  In both cases, the caller is responsible of managing all
the objects owned by the framework in the graph.

\dot
    digraph selection_overview {
        subgraph cluster_framework {
            label = "analysis framework"
            analysisframework [label="framework", fillcolor=grey75, style=filled]
            options [label="options collection", fillcolor=grey75, style=filled]
            analysistool [label="analysis tool"]
        }
        subgraph cluster_core {
            label = "core engine"
            labelloc = b
            selectioncollection [label="selection collection", fillcolor=grey75, style=filled]
            selectiondata [label="internal selection data", fillcolor=grey75, style=filled]
            selection [label="selection object"]
        }
        selectionoption [label="selection option"]
        selectionoptionmanager [label="selection option manager", fillcolor=grey75, style=filled]

        selectioncollection -> selection [label="creates"]
        selectioncollection -> selectiondata [label="owns and updates"]
        selectionoption -> selection [label="returns"]
        selection -> selectiondata [label="reads data", constraint=false]
        selectionoptionmanager -> selectionoption [label="provides values to"]
        selectionoptionmanager -> selectioncollection [label="gets selection objects"]
        analysistool -> selectionoption [label="declares"]
        analysistool -> selection [label="reads data from"]
        analysistool -> options [label="declares options"]
        analysisframework -> selectionoptionmanager [label="owns"]
        analysisframework -> selectioncollection [label="owns"]
        analysisframework -> options [label="owns"]
        options -> selectionoption [label="contains"]
    }
\enddot

 - _selection_: Evaluates to a single list of _selection positions_.
   Note in particular that the output is positions, not a list of atoms.
   A tool can accept one or more selections, and expect different semantics for
   different selections.
 - _dynamic selection_: The word _dynamic_ refers to selections for which the
   set of positions (instead of only the positions themselves) depends on the
   input coordinates.
 - _selection position_: A single coordinate as returned by a selection.
   This can correspond to an individual atom, but also to a collective
   coordinate such as a center of mass of a group of atoms.
   In addition to the output coordinates, the position provides information
   about the atoms that constitute it, and metadata that allows one to
   associate positions between different frames if different positions
   are returned at the same time.
 - _selection collection_: Group of selections that are processed as a unit
   against the same topology and sets of coordinates.
   In the analysis framework, there is always a single selection collection
   managed by the framework.
 - _selection variable_: When providing selections through text, it is possible
   to create variables and use them as part of selections.  This makes it
   easier to write repetitive selections by making complex common
   subexpressions into variables.  This also provides optimization
   opportunities for the selection engine: the variable value is not repeatedly
   evaluated.  Variables always exist in the context of a selection collection.
 - _selection object_: When a selection is _parsed_ (see below), the selection
   collection gives a handle to the selection as a _selection object_.  This
   handle is valid for the lifetime of the selection collection, and can be
   used to access information about the selection.  Operations on the selection
   collection (_compilation_ and _evaluation_, see below) alter the values
   returned by the selection objects.
 - _selection option_: A special type of command-line option that directly
   returns selection objects.  This higher-level construct is used by the
   analysis framework to provide a convenient interface for analysis tools:
   they can simply declare one or more selection options, and get a list of
   _selection objects_ as a return value for each of these.  Other parts of the
   selection engine are managed by the framework.

Core selection engine
=====================

The core of the selection engine is the _selection collection_ object.
The graph below shows how it handles selections.  The operations that the
collection object supports and their sequence is shown in the boxes in the
middle.  Inputs are shown at top, and outputs at the bottom.

\dot
    digraph selection_process {
        subgraph cluster_collection {
            label = "selection collection"
            subgraph actions {
                rank = same
                create [shape=box]
                parse [shape=box]
                compile [shape=box]
                evaluate [shape=box]
                evaluatefinal [label="finish evaluation",shape=box]

                create -> parse
                parse:ne -> parse:nw
                parse -> compile
                compile -> evaluate
                evaluate:ne -> evaluate:nw
                evaluate -> evaluatefinal
            }
            selectiondata [label="internal selection data"]
        }

        selectiontext [label="selection text"]
        topology [label="topology/\natom count"]
        frame [label="frame coordinates"]

        selection [label="selection object"]

        selectiontext -> parse
        parse:s -> selection:nw [label="returns"]
        parse -> selectiondata [label="creates"]
        topology -> compile
        compile -> selectiondata [label="initializes\npositions"]
        frame -> evaluate
        evaluate -> selectiondata [label="sets\npositions"]
        evaluatefinal -> selectiondata [label="resets to\npost-compilation\nstate"]
        selectiondata -> selection [label="reads data", dir=back]
    }
\enddot

 - _parsing_: after creating an empty selection collection,
   selections need to be parsed from text.  As a result, the selection
   collection initializes an internal data object to hold some basic
   information about the selections, and returns _selection objects_ as a
   handle to this data.  It is possible to parse more than one set of
   selections into the same collection by calling the parsing methods more than
   once.  The input string to parsing can also contain variable declarations,
   which get added into the collection and can be used in later calls to the
   parser.

 - _compilation_: when all selections are parsed, the whole selection
   collection is compiled.  This analyzes the provided selections, and
   evaluates all parts that do not depend on atom coordinates (e.g.,
   (sub)selections based on atom or residue names).  After compilation,
   the coordinates in the output positions are not initialized, but all other
   information is initialized as if all atoms satisfied any dynamic conditions.
   This means that any subsequent evaluation will return a subset of the
   positions returned at this point.  The caller can use this information to
   check the selections for validity and allocate memory for its own
   processing.
   Compilation also allocates all the memory necessary to do the evaluation.

   In the figure, topology is shown as input to the compilation, but generally
   it can be set at any point before the compilation.  If the selection text
   does not require any information from the topology for evaluation, it is
   sufficient to set only the atom count.

 - _evaluation_: after the selections are compiled, they can be evaluated for
   one or more sets of atom coordinates.  This updates the set of positions
   accessible through the selection objects.  For dynamic selections, the group
   of positions can change; for static selections, only the coordinates of the
   positions are updated.

 - _final evaluation_: This returns the selections to the state they were after
   compilation, i.e., to the maximum possible set of positions.  The
   coordinates of the positions are again uninitialized, but other information
   is available.  The caller can use this information to do post-processing
   and, e.g., produce labels in its output based on the selection positions.

\if internal
Internal implementation
=======================

Implementation details of different parts of the module are discussed on
separate pages.

  - \subpage page_module_selection_custom
  - \subpage page_module_selection_parser
  - \subpage page_module_selection_compiler
  - \subpage page_module_selection_insolidangle
\endif
