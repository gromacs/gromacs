:orphan:

GROMACS
=======

Description
-----------

**GROMACS** is a full-featured suite of programs to perform molecular dynamics
simulations - in other words, to simulate the behavior of systems with hundreds
to millions of particles, using Newtonian equations of motion.
It is primarily used for research on proteins, lipids, and polymers, but can be
applied to a wide variety of chemical and biological research questions.

Synopsis
--------

The following commands make up the GROMACS suite.  Please refer to their
individual man pages for further details.

.. include:: bytopic.rst

Additional documentation
------------------------

Consult the manual at <http://www.gromacs.org> for an introduction to molecular
dynamics in general and GROMACS in particular, as well as an overview of the
individual programs.

References
----------

The development of GROMACS is mainly funded by academic research grants.
To help us fund development, the authors humbly ask that you cite the GROMACS papers:

H.J.C. Berendsen, D. van der Spoel and R. van Drunen.  **GROMACS: A message-passing
parallel molecular dynamics implementation**.  Comp. Phys. Comm. *91*, 43-56 (1995)

Erik Lindahl, Berk Hess and David van der Spoel.  **GROMACS 3.0: A package for
molecular simulation and trajectory analysis**.  J. Mol. Mod. *7*, 306-317 (2001)

B. Hess, C. Kutzner, D. van der Spoel, and E. Lindahl.  **GROMACS 4: Algorithms for
Highly Efficient, Load-Balanced, and Scalable Molecular Simulation**.
J. Chem. Theory Comput. *4*, 3, 435-447 (2008),
<http://dx.doi.org/10.1021/ct700301q>

Authors
-------

Current developers:

David van der Spoel <spoel@gromacs.org>  
Berk Hess <hess@gromacs.org>  
Erik Lindahl <lindahl@gromacs.org>

A full list of present and former contributors
is available at <http://www.gromacs.org>

This manual page is largely based on the GROMACS online reference, and the text
was prepared originally by Nicholas Breen <nbreen@ofb.net>.

Bugs
----

GROMACS has no major known bugs, but be warned that it stresses your CPU more
than most software.  Systems with slightly flaky hardware may prove unreliable
while running heavy-duty simulations.  If at all possible, please try to
reproduce bugs on another machine before reporting them.
