Module dependency graph {#page_modulegraph}
=======================

The graph below shows the dependencies between the source code modules,
computed from include statements in the code.
For documented modules (those that do not have a gray background), clicking on
the module takes you to the module documentation.
Legend for the graph can be found below the graph.

\ifnot xml
\dotfile module-deps.dot
\endif

Legend
======

The graph below annotates the colors and line styles used in the module
dependency graph above.  More detailed textual annotation is below the graph.

\dot
digraph legend {
    node [fontname="FreeSans",fontsize=10,height=.2,shape=box]
    edge [fontname="FreeSans",fontsize=10]
    rankdir = "LR"
    subgraph cluster_nodes {
        label = "Nodes"
        legacy    [label="undocumented", fillcolor=grey75, style="filled"]
        analysis  [label="analysis", fillcolor="0 .2 1", style="filled"]
        utility   [label="utility", fillcolor=".08 .2 1", style="filled"]
        mdrun     [label="mdrun", fillcolor=".75 .2 1", style="filled"]
        installed [label="installed", color=".66 .5 1", penwidth=3]
    }
    subgraph cluster_edges {
        label = "Edges"
        node [label="<any>"]
        invalid1 -> invalid2 [label="invalid", color=red]
        legacy1 -> legacy2 [label="legacy", color=grey75]
        legacy2 [label="undoc"]
        public1 -> public2 [label="public", color=black]
        public1 [label="public"]
        public2 [label="public"]
        library1 -> library2 [label="library", color=".66 .8 .8"]
        library1 [label="library"]
        pubimpl1 -> pubimpl2 [color=black, style=dashed]
        pubimpl1 [label="internal"]
        pubimpl2 [label="public"]
        libimpl1 -> libimpl2 [color=".66 .8 .8", style=dashed]
        libimpl1 [label="internal"]
        libimpl2 [label="library"]
        test1 -> test2 [label="test", color=".33 .8 .8", style=dashed]
        test1 [label="test"]
    }
    legacy -> invalid1 [style="invis"]
}
\enddot

Node colors:
<dl>
<dt>gray background</dt>
<dd>undocumented module</dd>
<dt>orange background</dt>
<dd>documented utility modules</dd>
<dt>red background</dt>
<dd>documented analysis modules</dd>
<dt>violet background</dt>
<dd>documented MD execution modules</dd>
<dt>light blue border</dt>
<dd>module contains public API (installed) headers</dd>
</dl>

Edge colors (an edge with a certain color indicates that types above it in the
list are not present):
<dl>
<dt>red</dt>
<dd>invalid dependency</dd>
<dt>gray</dt>
<dd>legacy dependency
(dependency on undocumented file, or to undocumented directories)</dd>
<dt>solid black</dt>
<dd>public header depends on the other module</dd>
<dt>solid blue</dt>
<dd>library header depends on the other module</dd>
<dt>dashed blue</dt>
<dd>source file depends on library header in the other module</dd>
<dt>dashed black</dt>
<dd>source file depends on public header in the other module</dd>
<dt>dashed green</dt>
<dd>test file depends on the other module</dd>
</dl>
