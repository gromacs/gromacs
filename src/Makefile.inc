## Process this file with automake to produce Makefile.in
#
# Don't edit - this file is generated automatically from Makefile.am
#

#################
# We need a second compile command producing executables 
# that can be run on the local host to make the innerloops.
# In most cases this is identical to the normal one, but not
# for cross-compilation. We also need it on parallel machines 
# where the MPI executables cannot be run outside a batch queue.

BUILD_COMPILE = $(BUILD_CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)

#################
# This might be bad - but I don't know any other way to enable
# us to type make <progname> in subdirs right now
#AM_CPPFLAGS = -DHAVE_CONFIG_H


