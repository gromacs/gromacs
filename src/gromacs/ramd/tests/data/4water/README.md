# Generate tpr file

```
cd /workspaces/gromacs-ramd/src/gromacs/ramd/tests/data/4water
/workspaces/gromacs-ramd/build/bin/gmx grompp -f ramd.mdp -c /workspaces/gromacs-ramd/src/testutils/simulationdatabase/4water.gro -p /workspaces/gromacs-ramd/src/testutils/simulationdatabase/4water.top -n /workspaces/gromacs-ramd/src/testutils/simulationdatabase/4water.ndx -o topol.tpr
```