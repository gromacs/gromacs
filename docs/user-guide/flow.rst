Flow Chart
==========

This is a flow chart of a typical |Gromacs| MD run of a protein
in a box of water.
A more detailed example is available in :doc:`getting-started`.
Several steps of energy minimization may be necessary,
these consist of cycles: :ref:`gmx grompp` -> :ref:`gmx mdrun`.

.. digraph:: flowchart

   node [ shape=box, width=1.5 ]

   input_pdb [
     label="eiwit.pdb"
     tooltip="Protein Databank file"
     URL="file-formats.html#pdb"
     shape=none, width=0, height=0, margin=0
     group=input
   ]
   pdb2gmx [
     label="Generate a GROMACS topology\ngmx pdb2gmx"
     tooltip="Convert PDB file to GROMACS coordinate file and topology"
     URL="../onlinehelp/gmx-pdb2gmx.html"
     width=3
     group=main
   ]

   input_pdb -> pdb2gmx [ headport=e ]

   editconf [
     label="Enlarge the box\ngmx editconf"
     tooltip="Adjust box size and placement of molecule"
     URL="../onlinehelp/gmx-editconf.html"
   ]

   pdb2gmx -> editconf [
     label="conf.gro"
     labeltooltip="GROMACS coordinate file containing molecules from PDB file"
     URL="file-formats.html#gro"
   ]

   solvate [
     label="Solvate protein\ngmx solvate"
     tooltip="Fill box with water (solvate molecule)"
     URL="../onlinehelp/gmx-solvate.html"
     width=3
     group=main
   ]

   pdb2gmx -> solvate [
     label="topol.top"
     labeltooltip="GROMACS ascii topology file"
     URL="file-formats.html#top"
   ]
   editconf -> solvate [
     label="conf.gro"
     labeltooltip="GROMACS coordinate file with adjusted box etc."
     URL="file-formats.html#gro"
   ]

   input_mdp [
     label="grompp.mdp"
     tooltip="Parameter file from grompp (controls all MD parameters)"
     URL="file-formats.html#mdp"
     shape=none, width=0, height=0, margin=0
     group=input
   ]
   grompp [
     label="Generate mdrun input file\ngmx grompp"
     tooltip="Process parameters, coordinates and topology and write binary topology"
     URL="../onlinehelp/gmx-grompp.html"
     width=3
     group=main
   ]

   input_pdb -> input_mdp [ style=invis, minlen=3 ]

   input_mdp -> grompp [ headport=e, weight=0 ]
   solvate -> grompp [
     label="conf.gro"
     labeltooltip="GROMACS coordinate file with water molecules added"
     URL="file-formats.html#gro"
   ]
   solvate -> grompp [
     label="topol.top"
     labeltooltip="GROMACS ascii topology file with water molecules added"
     URL="file-formats.html#top"
   ]

   mdrun [
     label="Run the simulation (EM or MD)\ngmx mdrun"
     tooltip="The moment you have all been waiting for! START YOUR MD RUN"
     URL="../onlinehelp/gmx-mdrun.html"
     width=3
     group=main
   ]

   grompp -> mdrun [
     label="topol.tpr"
     labeltooltip="Portable GROMACS binary run input file (contains all information to start MD run)"
     URL="file-formats.html#tpr"
   ]
   mdrun -> mdrun [
     label="Continuation\nstate.cpt"
     labeltooltip="Checkpoint file"
     URL="file-formats.html#cpt"
   ]

   analysis [
     label="Analysis\ngmx ...\ngmx view"
     tooltip="Your favourite GROMACS analysis tool"
     URL="cmdline.html#commands-by-topic"
   ]

   mdrun -> analysis [
     label="traj.xtc / traj.trr"
     labeltooltip="Portable compressed trajectory / full precision portable trajectory"
     URL="file-formats.html#xtc"
   ]

   energy [
     label="Analysis\ngmx energy"
     tooltip="Energy plots, averages and fluctuations"
     URL="../onlinehelp/gmx-energy.html"
   ]

   mdrun -> energy [
     label="ener.edr"
     labeltooltip="Portable energy file"
     URL="file-formats.html#edr"
   ]
