Run parameters and Programs
===========================

Online documentation
--------------------

:raw-latex:`\index{online documentation}` More documentation is
available online from the GROMACS web site,
http://manual.gromacs.org/documentation.

In addition, we install standard UNIX man pages for all the programs. If
you have sourced the :raw-latex:`\tt`GMXRC script in the GROMACS binary
directory for your host they should already be present in your
:raw-latex:`\tt`MANPATH environment variable, and you should be able to
type *e.g.* :raw-latex:`\tt`man gmx-grompp. You can also use the
:raw-latex:`\tt`-h flag on the command line (e.g. :raw-latex:`\tt`gmx
grompp -h) to see the same information, as well as :raw-latex:`\tt`gmx
help grompp. The list of all programs are available from
:raw-latex:`\tt`gmx help.

File types:raw-latex:`\index{file type}`:raw-latex:`\index{type, file|see{file type}}`
--------------------------------------------------------------------------------------

Table :raw-latex:`\ref{tab:form}` lists the file types used by GROMACS
along with a short description, and you can find a more detail
description for each file in your HTML reference, or in our online
version.

GROMACS files written in XDR:raw-latex:`\index{XDR}` format can be read
on any architecture with GROMACS version 1.6 or later if the
configuration script found the XDR libraries on your system. They should
always be present on UNIX since they are necessary for NFS support.

Run Parameters:raw-latex:`\index{run parameter}`:raw-latex:`\index{parameter, run|see{run parameter}}`
------------------------------------------------------------------------------------------------------

The descriptions of :raw-latex:`\tt`.mdp parameters can be found at
http://manual.gromacs.org/current/mdp-options.html or in your
installation at :raw-latex:`\tt`share/gromacs/html/mdp-options.html

