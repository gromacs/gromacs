# Sun/SOLARIS
#	@(#) Makefile.sol 1.47 4/16/97
#
#	GROMACS - Groningen Machine for Chemical Simulation
#	Copyright (c) 1990, 1991, 1992, Groningen University
#
#	Makefile for gromacs on sparc processor
#
# 	See README file for info
#
SYSDEFS		= -DHAVE_IDENT -DHAVE_STRDUP -DHAVE_STRCASECMP

# Sun C compilers
CC		= cc
CCC		= CC
F77		= f77

# DO NOT PUT OPTIMIZATION HIGHER FOR SparcStation 5, it will produce
# erroneous results! -xO2 is safe.
CFLAGS		= -v -xO2 -xunroll=3 -Xa -K PIC
#CFLAGS		= -v -g -Xa -K PIC
CCFLAGS		= -O4
FFLAGS		= -O4 -PIC

# Generic linking stuff
LDFLAGS		= -L$(LIBDIR) -L/usr1/local/lib
LD		= $(CC)  $(LDFLAGS) -z nodefs  
FLD		= $(F77) $(LDFLAGS)
CCLD		= $(CCC) $(LDFLAGS) -z nodefs 

XLIBS		= -lsocket -lX11
SYSLIBS		= $(LINKLIB) -lm -lnsl 
ARFLAGS		= ur
SHAREIT		= (cd $(LIBDIR); ar x $(LIB); cc $(LDFLAGS) -o $(SOLIB) -G *.o; $(RM) *.o)
RANLIB		= echo
X11INC		= -I/usr1/local/include

#
#	USER MODIFIABLE SECTION
#
# If you want to use fortran innerloops set this to yes
# For most machines this will improve the performance quite a bit
# because C-compilers are not as good as Fortran compilers
USEF77		= yes
SYSDEFS		+= -DF77UNDERSCORE
#
# If you want to run in *** P A R A L L E L ***
# please select either PVM or MPI, check with your local hacker
# to see what's installed at your site. If you have neither,
# set both to no. If you select both, something will break!!!
#
USE_PVM3	= no
USE_MPI		= no
#
# If you want to use compressed data files set this to yes
# This uses the xdr libraries of your UNIX system, which is virtually
# Allways present, because it is necessary for NFS (network file system)
#
USE_XDR		= yes
#
# Graphics Library
# Set this on if you have Open GL and the Viewkit library
# (standard with SGI, available for money elsewhere)
# This is used in the Open GL trajectory viewer.
HAVE_GL		= no
#
# Note that these variables are also used in Makefile.std
# If something does not work, please check out the link command line
# in that file (e.g. for PVM)
#

