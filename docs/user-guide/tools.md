# Tools in GROMACS  

GROMACS includes many tools for preparing, running and analysing
molecular dynamics simulations. Documentation for these can be found
at the links below.

* Tools documentation [by name](programs/byname.html)  
* Tools documentation [by topic](programs/bytopic.html)  

## Common options and behaviour of GROMACS tools

Optional files are not used unless the option is set, in contrast to
non optional files, where the default file name is used when the
option is not set.

All GROMACS tools will accept file options without a file extension
or filename being specified. In such cases the default filenames will
be used. With multiple input file types, such as generic structure
format, the directory will be searched for files of each type with the
supplied or default name. When no such file is found, or with output
files the first file type will be used.

All GROMACS tools with the exception of [mdrun], [trjcat] and
[eneconv] check if the command line options are valid.  If this is not
the case, the program will be halted.

All GROMACS tools have 4 hidden options:

       option  type  default  description
-------------  ----  -------  ------------------
`-[no]hidden`  bool  [`yes`]  [hidden] Print description for hidden options
`-[no]quiet`   bool  [` no`]  [hidden] Do not print help info
`-[no]debug`   bool  [` no`]  [hidden] Write file with debug information

Many tools accept enumerated options (enum), which should be used with
one of the arguments listed in the option description. The argument
may be abbreviated, and the first match to the shortest argument in
the list will be selected.

Many tools also use options that may accept a vector of values. Either
1 or 3 parameters can be supplied; when only one parameter is supplied
the two other values are also set to this value.

All GROMACS tools can read compressed or g-zipped files. There might
be a problem with reading compressed [.tng] or [.xtc] files, but these
will not compress very well anyway.

Most GROMACS tools can process a trajectory with fewer atoms than the
run input or structure file, but only if the trajectory consists of
the first n atoms of the run input or structure file.


