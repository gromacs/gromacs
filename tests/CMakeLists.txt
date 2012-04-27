add_test(NAME TestExec_mdrun-h
         COMMAND mdrun -h)
add_dependencies(tests mdrun)
