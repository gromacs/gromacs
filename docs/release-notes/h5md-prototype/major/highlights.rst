Highlights
^^^^^^^^^^

|Gromacs| H5MD prototype was released on August 30th, 2024, as a milestone in the
`MDDB <https://mddbr.eu>`_ EU project.

Patch releases of the prototype may have been made since then, please use the updated versions!

The H5MD prototype is based on the |Gromacs| development branch (as per August 29th,
2024), soon to evolve into the 2025 beta release. The H5MD prototype branch will not
describe any of the features or improvements to be delivered in the |Gromacs| 2025 release.

:ref:`h5md` is a new (in |Gromacs|) file specification using the HDF5 file format.
An overview of the H5MD features, supported in |Gromacs|, are:

* Read and write coordinates, simulation box, velocities and forces.
   * Lossy compression of coordinates using the SZ3 HDF5 plugin.
   * Lossless compression.
* Continue from checkpoints.
* Write connectivity records.
* Write a custom |Gromacs| topology record.
* Write a separate topology record compatible with the VMD H5MD plugin (**not feature complete**).

Keep in mind that this is a prototype release. There are no guarantees that the
H5MD files produced by this version will be compatible with future versions of
|Gromacs|, regarding the file contents and the compression formats.

HDF5 libraries (version >= 1.10.1) are required to access the H5MD features.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!
