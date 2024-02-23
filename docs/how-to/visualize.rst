.. _reference manual: gmx-manual-parent-dir_

.. _gmx-visualize:

Visualization Software
----------------------

Some programs that are useful for visualizing either a trajectory file and/or a coordinate file are:

* `VMD`_ - a molecular visualization program for displaying, animating, and analyzing
  large biomolecular systems using 3-D graphics and built-in scripting. Reads |Gromacs| trajectories.
* `PyMOL`_ - capable molecular viewer with support for animations, high-quality rendering, crystallography,
  and other common molecular graphics activities. Does not read |Gromacs| trajectories in default
  configuration, requiring conversion to PDB or similar format. When compiled with `VMD`_ plugins,
  :ref:`trr` & :ref:`xtc` files can be loaded.
* `Rasmol`_ - the derivative software `Protein Explorer`_ (below) might be a better alternative, but
  the Chime component requires windows. `Rasmol`_ works fine on Unix.
* `Protein Explorer`_ - a `RasMol`_\ -derivative, is the easiest-to-use and most powerful software
  for looking at macromolecular structure and its relation to function. It runs on Windows or Macintosh/PPC computers.
* `Chimera`_ - A full featured, Python-based visualization program with all sorts of features for
  use on any platform. The current version reads |Gromacs| trajectories.
* `Molscript`_ - This is a script-driven program form high-quality display of molecular 3D structures
  in both schematic and detailed representations. You can get an academic license for free from Avatar.

Topology bonds vs Rendered bonds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Remember that each of these visualization tools is only looking at the coordinate file you gave it.
Thus it's not using your topology which is described in either your :ref:`top` file or your 
:ref:`tpr` file. Each of these programs makes their own guesses about where the chemical bonds 
are for rendering purposes, so do not be surprised if the heuristics do not always match your topology.

.. _Rasmol: http://www.umass.edu/microbio/rasmol/index2.htm
.. _Protein Explorer: http://www.umass.edu/microbio/rasmol/
.. _Chimera: http://www.rbvi.ucsf.edu/chimera/
.. _Molscript: https://github.com/pekrau/MolScript


Extracting Trajectory Information
---------------------------------

There are several techniques available for finding information in |Gromacs|
trajectory (:ref:`trr`, :ref:`xtc`, :ref:`tng`) files.

* use the |Gromacs| trajectory analysis utilities
* use :ref:`gmx traj` to write a :ref:`xvg` file and read that in an external program as above
* write your own C code using ``gromacs/share/template/template.cpp`` as a template
* use :ref:`gmx dump` and redirect the shell output to a file and read that in an external
  program like MATLAB, or Mathematica or other spreadsheet software.

External tools to perform trajectory analysis
---------------------------------------------

In recent years several external tools have matured sufficiently to analyse diverse sets
of trajectory data from several simulation packages. Below is a short list of tools (in an alphabetical order)
that are known to be able to analyse |Gromacs| trajectory data.

* `LOOS <http://loos.sourceforge.net/>`__
* `MDAnalysis <https://www.mdanalysis.org/>`__
* `MDTraj <http://mdtraj.org/>`__
* `Pteros <https://github.com/yesint/pteros/>`__


Plotting Data
-------------

The various |Gromacs| analysis utilities can generate :ref:`xvg` files. These are text files
that have been specifically formatted for direct use in Grace. You can, however, in
all |Gromacs| analysis programs turn off the Grace specific codes by running the programs
with the ``-xvg none`` option. This circumvents problems with tools like gnuplot and Excel (see below).

Note that Grace uses some embedded backslash codes to indicate superscripts, normal script, etc. in units. So "Area (nm\S2\N)" is nm squared. 

Software
^^^^^^^^

Some software packages that can be used to graph data in a :ref:`xvg` file:

* Grace - WYSIWYG 2D plotting tool for the X Window System and M\*\ tif. Grace runs on practically
  any version of Unix-like OS, provided that you can satisfy its library dependencies (Lesstif is a
  valid free alternative to Motif). It is also available for the other common operation systems.
* gnuplot - portable command-line driven interactive data and function plotting utility for UNIX,
  IBM OS/2, MS Windows, DOS, Macintosh, VMS, Atari and many other platforms. Remember to use::

    set datafile commentschars "#@&"

  to avoid gnuplot trying to interpret Grace-specific commands in the :ref:`xvg` file or use
  the ``-xvg none`` option when running the analysis program. For simple usage,::

    plot "file.xvg" using 1:2 with lines

  is a hack that will achieve the right result.
* Matplotlib - a popular Python library for visualization. A simple script that will plot the data
  in ``file.xvg`` and show the result on the screen

  .. code-block:: python

      import numpy as np
      import matplotlib.pyplot as plt
      x, y = np.loadtxt("file.xvg", comments=["@", "#", "&"], unpack=True)
      plt.plot(x, y)
      plt.show()

* MS Excel - change the file extension to .csv and open the file (when prompted, choose to ignore the
  first 20 or so rows and select fixed-width columns, if you are using German MS Excel version, you
  have to change decimal delimiter from "," to ".", or use your favourite \*nix tool.
* Sigma Plot A commercial tool for windows with some useful analysis tools in it.
* R - freely available language and environment for statistical computing and graphics which provides
  a wide variety of statistical and graphical techniques: linear and nonlinear modelling, statistical
  tests, time series analysis, classification, clustering, etc.
* SPSS A commercial tool (Statistical Product and Service Solutions), which can also plot and analyse data.


Micelle Clustering
------------------

This is necessary for the :ref:`gmx spatial` tool if you have a fully-formed single aggregate and
want to generate the spatial distribution function for that aggregate or for solvent around that aggregate.

Clustering to ensure that the micelle is not split across a :ref:`periodic boundary condition <gmx-pbc>`
border is an essential step prior to calculating properties such as the radius of gyration and the
radial distribution function. Without this step your results will be incorrect (a sign of this error
is unexplained huge fluctuations in the calculated value when the visualized trajectory looks fine).

Three steps are required:

* use :ref:`trjconv <gmx trjconv>` ``-pbc cluster`` to obtain a single frame that has all of the
  lipids in the unit cell. This must be the first frame of your trajectory. A similar frame
  from some previous timepoint will not work.
* use :ref:`grompp <gmx grompp>` to make a new :ref:`tpr` file based on the frame that was output from the step above.
* use :ref:`trjconv <gmx trjconv>` ``-pbc nojump`` to produce the desired trajectory using the newly produced :ref:`tpr` file.

More explicitly, the same steps are:

::

 gmx trjconv -f a.xtc -o a_cluster.gro -e 0.001 -pbc cluster
 gmx grompp -f a.mdp -c a_cluster.gro -o a_cluster.tpr
 gmx trjconv -f a.xtc -o a_cluster.xtc -s a_cluster.tpr -pbc nojump


