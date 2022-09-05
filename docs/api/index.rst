=======
C++ API
=======

.. toctree::
    :hidden:

    gmxlibs

Public C++ application programming interfaces are available for |Gromacs| installations
depending on the detected environment and user options when the |Gromacs| build
is configured with CMake.

* :doc:`gmxlibs`
    * CMake target ``Gromacs::gmxapi``,
      enabled by :cmake:`GMXAPI`
      (default, when :cmake:`BUILD_SHARED_LIBS` on non-Windows platforms),
      provides :file:`gmxapi/` headers and ``::gmxapi`` C++ namespace.
    * CMake target ``Gromacs::libgromacs``,
      enabled by :cmake:`GMX_INSTALL_LEGACY_API` (default ``OFF``),
      provides :file:`gromacs/` headers and ``::gmx`` C++ namespace.

* :ref:`nblib`: Enabled by :cmake:`GMX_INSTALL_NBLIB_API`.
  (default, when :cmake:`BUILD_SHARED_LIBS` on non-Windows platforms)

.. only:: html

    * `Legacy API <../doxygen/html-user/index.xhtml>`_:
      Enabled by :cmake:`GMX_INSTALL_LEGACY_API`. (default ``OFF``)
