Flow Chart
==========

This is a flow chart of a typical GROMACS MD run of a protein
in a box of water.
A more detailed example is available in :doc:`getting-started`.
Several steps of energy minimization may be necessary,
these consist of cycles: :ref:`gmx grompp` -> :ref:`gmx mdrun`.

.. raw:: html

   <CENTER>
   <TABLE BORDER=0 CELLMARGIN=0 CELLPADDING=0 CELLSPACING=0>

   <TR>
   <TD COLSPAN="2" STYLE="text-align:right"> <A HREF="file-formats.html#pdb" onMouseOver="window.status='Protein Databank file'; return true">eiwit.pdb</A></TD>
   <TD STYLE="text-align:right"><IMG SRC=../images/flow_leftrightdown.gif></TD>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   </TR>

   <TR>
   <TD COLSPAN=2>Generate a GROMACS topology </TD>
   <TD></TD>
   <TD BGCOLOR=#777777 COLSPAN=3 STYLE="text-align:center"><A HREF=../programs/gmx-pdb2gmx.html onMouseOver="window.status='Convert PDB file to GROMAX coordinate file and topology'; return true"><B>gmx pdb2gmx</B></A></TD>
   <TD><IMG SRC=../images/flow_vrule.gif></TD>
   </TR>


   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif BORDER=0></TD>
   <TD WIDTH=20></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><A HREF="file-formats.html#gro" onMouseOver="window.status='GROMACS coordinate file containing molecules from PDB file'; return true">conf.gro</A></TD>
   <TD></TD>
   <TD STYLE="text-align:center"> <A HREF="file-formats.html#top" onMouseOver="window.status='GROMACS ascii topology file'; return true">topol.top</A> </TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif BORDER=0></TD>
   <TD></TD>
   <TD ROWSPAN=5 COLSPAN=1 STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif><BR><IMG SRC=../images/flow_vline.gif><BR><IMG SRC=../images/flow_vline.gif><BR><IMG SRC=../images/flow_vline.gif><BR><IMG SRC=../images/flow_down.gif></TD>
   </TR>

   <TR>
   <TD COLSPAN=2>Enlarge the box</TD>
   <TD></TD>
   <TD STYLE="text-align:center" BGCOLOR=#777777><A HREF=../programs/gmx-editconf.html onMouseOver="window.status='Adjust boxsize and placement of molecule'; return true"><B>gmx editconf</B></A></TD>
   <TD></TD>
   <TD><IMG SRC=../images/flow_vrule.gif></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   <TD></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"> <A HREF="file-formats.html#gro" onMouseOver="window.status='GROMACS coordinate file with adjusted box etc.'; return true">conf.gro</A> </TD>
   <TD></TD>
   </TR>


   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   <TD></TD>
   </TR>

   <TR>
   <TD COLSPAN=2>Solvate protein</TD>
   <TD></TD>
   <TD COLSPAN=3 STYLE="text-align:center" BGCOLOR=#777777>&nbsp;<A HREF=../programs/gmx-solvate.html onMouseOver="window.status='Fill box with water (solvate molecule)'; return true"><B>gmx solvate</B></A>&nbsp;</TD>
   <TD><IMG SRC=../images/flow_vrule.gif></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><A HREF="file-formats.html#gro" onMouseOver="window.status='GROMACS coordinate file with water molecules added'; return true">conf.gro</A></TD>
   <TD></TD>
   <TD STYLE="text-align:center"> <A HREF="file-formats.html#top" onMouseOver="window.status='GROMACS ascii topology file with water molecules added'; return true">topol.top</A> </TD>
   </TR>

   <TR>
   <TD COLSPAN=2 STYLE="text-align:right"><A HREF="file-formats.html#mdp" onMouseOver="window.status='Parameter file for grompp (controls all MD parameters)'; return true">grompp.mdp</A></TD>
   <TD STYLE="text-align:right">&nbsp;<IMG SRC=../images/flow_leftrightdown.gif></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   <TD></TD>
   </TR>


   <TR>
   <TD COLSPAN=2>Generate mdrun input file</TD>
   <TD></TD>
   <TD COLSPAN=3 STYLE="text-align:center" BGCOLOR=#777777><A HREF=../programs/gmx-grompp.html onMouseOver="window.status='Process parameters, coordinates and topology and write binary topology'; return true"><B>gmx grompp</B></A></TD>
   <TD><IMG SRC=../images/flow_vrule.gif></TD>
   <TD></TD>
   <TD></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   <TD ROWSPAN=3 STYLE="text-align:right"><IMG SRC=../images/flow_rightleftdown.gif></TD>
   <TD STYLE="text-align:center;vertical-align:bottom">Continuation</TD>
   </TR>

   <TR>
   <TD COLSPAN=2></TD>
   <TD></TD>
   <TD COLSPAN=3 STYLE="text-align:center"> <A HREF="file-formats.html#tpr" onMouseOver="window.status='Portable GROMACS binary run input file (contains all information to start MD run)'; return true">topol.tpr</A></TD>
   <TD STYLE="text-align:center"><A HREF="file-formats.html#cpt" onMouseOver="window.status='Checkpoint file'; return true">state.cpt</A></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   <TD ROWSPAN=2 STYLE="text-align:center">
   <IMG SRC=../images/flow_vline.gif><BR>
   <IMG SRC=../images/flow_leftup.gif></TD>
   </TR>

   <TR>
   <TD COLSPAN=2>Run the simulation (EM or MD)</TD>
   <TD></TD>
   <TD  COLSPAN=3 STYLE="text-align:center" BGCOLOR=#777777>&nbsp;<A HREF=../programs/gmx-mdrun.html onMouseOver="window.status='The moment you have all been waiting for! START YOUR MD RUN'; return true"><B>gmx mdrun</B></A>&nbsp;</TD>
   <TD></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_vline.gif></TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"> <A HREF="file-formats.html#xtc" onMouseOver="window.status='Portable compressed trajectory'; return true">traj.xtc</A> /
   <A HREF="file-formats.html#trr" onMouseOver="window.status='Full precision portable trajectory'; return true">traj.trr</A> </TD>
   <TD></TD>
   <TD STYLE="text-align:center"> <A HREF="file-formats.html#edr" onMouseOver="window.status='Portable energy file'; return true">ener.edr</A> </TD>
   </TR>

   <TR>
   <TD></TD>
   <TD></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   <TD></TD>
   <TD STYLE="text-align:center"><IMG SRC=../images/flow_down.gif></TD>
   </TR>

   <TR>
   <TD COLSPAN=2>Analysis</TD>
   <TD></TD>
   <TD STYLE="text-align:center" BGCOLOR=#777777><A HREF="../programs/bytopic.html" onMouseOver="window.status='Your favourite GROMACS analysis tool'; return true"><B>g_...</B></A><BR><A HREF=../programs/gmx-view.html onMouseOver="window.status='gmx view, the GROMACS trajectory viewer'; return true"><B>gmx view</B></A></TD>
   <TD></TD>
   <TD STYLE="text-align:center" BGCOLOR=#777777><A HREF=../programs/gmx-energy.html onMouseOver="window.status='Energy plots, averages and  fluctuations'; return true"><B>gmx energy</B></A></TD>
   <TD><IMG SRC=../images/flow_vrule.gif></TD>
   </TR>


   </TABLE>
   </CENTER>

