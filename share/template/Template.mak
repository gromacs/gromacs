# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 template makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the template one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
LIBS         = @LIBS@
LDFLAGS      = @LDFLAGS@
CFLAGS	     = @CFLAGS@	
CC           = @CC@
LD           = $(CC)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

template:	template.o
		$(LD) $(LDFLAGS) -o $@ template.o $(LIBS)
