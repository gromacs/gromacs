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
