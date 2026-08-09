# RCSB mmCIF Dictionary Suite

![example workflow](https://github.com/rcsb/cpp-dict-pack/workflows/Testing/badge.svg)

## Introduction

This repository/package contains the source for the mmCIF dictionary suite.  Included are a number of tools for validation of dictionaries and data.

## Installation

### Download the source code

* As a distribution with a filename like mmcif-dict-suite-vX.XXX-prod-src.tar.gz
```bash
   gzip -d -c mmcif-dict-suite-vX.XXX-prod-src.tar.gz | tar xf -
```
   
* From github 

```bash
git clone  --recurse-submodules  https://github.com/rcsb/cpp-dict-pack.git

```

### Building
To build the dictionary suite, you need to have the following in your path or development environment

* CMake
* python (version 2/3)
* C++ compiler
* bash
* csh
* flex
* bison

Typically, one created a build tree, uses cmake to configure, and then build.

```
mkdir build
cd build
cmake .. -D<options>
make
```

This will build the tools, and compile the dictionaries.

#### Configuration options
As CMake is used, the following command line settings may be used:

* MINIMAL\_DICTS:  This indicates that a subset of dictionaries should be downloaded and built. (use ```-DMINIMAL_DICTS=ON```).
* ALWAYS\_DOWNLOAD\_DICTS: When working with a pre-packaged tar distribution, only a subset of the dictionaries are provided. Use this option to force download of dictionaries from GitHub. (use ```-DALWAYS_DOWNLOAD_DICTS=ON```).

### Using the dictionary suite

#### Retrieving the latest dictionary

You may retrieve the latest version of the mmcif-pdbx dictionary by doing the following (from within the build directory)

```bash
mkdir PDBx
cd PDBx
wget http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
../bin/DictToSdb -ddlFile ../dicts/dict-mmcif_ddl/mmcif_ddl.dic \
      -dictFile mmcif_pdbx_v50.dic -dictSdbFile mmcif_pdbx_v50.sdb

```

If errors are found, parsing errors are stored in the file *mmcif_pdbx_v50.dic-parser.log* and validation errors, against the DDL, are stored in the file *mmcif_pdbx_v50.dic-diag.log*
    
#### Validating a mmCIF/PDBx file against the PDBx dictionary

This is an example of PDBx file validation, of entry 3q45, against the PDBx dictionary.

```bash
wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/q4/3q45.cif.gz
gunzip 3q45.cif.gz
./bin/CifCheck -f 3q45.cif -dictSdb PDBx/mmcif_pdbx_v50.sdb
```

If errors are found, parsing errors are stored in the file *3q45.cif-parser.log* and validation errors, against the dictionary, are stored in the file *3q45.cif-diag.log*.


## Notes for developers
Some compromises were made in building the suite and keeping dictionaries consistent

* Normally when you run *cmake*, dictionaries are downloaded from GitHub - unless already present in the *dicts* subdirectory.  If you remove a dictionary and run cmake, it will retrieve if need be.

* If you are doing dictionary development, a symlink to a development dictionary maintained in the dicts directory will be sufficient.  Dependencies will be updated properly.

* To create a distribution, use `make dist`

* *everything* option to make will build the sdb/odb/xml files.  It is protected so the *-j* option to make will work without collision.

