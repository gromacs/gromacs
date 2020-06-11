Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Configuration-time trivalue options changed from autodetection to boolean on/off
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
To simplify the CMake configuration and avoid having multiple settings that
change outside of the users direct control we have removed the support for
automatically setting booleans. GMX_BUILD_HELP and GMX_HWLOC are now
disabled by default, while GMX_LOAD_PLUGINS is enabled by default.
