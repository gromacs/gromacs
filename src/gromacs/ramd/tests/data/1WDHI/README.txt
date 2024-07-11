gmx grompp -f pull.mdp -c 1WDHI_box.gro -p topol.top -o topol.tpr
gmx dump -s topol.tpr | less
gmx mdrun -v
