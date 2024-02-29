## Version 3 (do not change this one)
MIFVN = 3
### USER - SPECIFIC PART OF THE MAKEFILE ####################
# list of all source files (without directory part and suffix)
# needed from source except include and module files in
# sequence of dependency; Main program file should be the last file
# (defaults to all files in source)
#
SOURCE_FILES =
#
# as above for source_<v> files:
#
SOURCE_V_FILES = main
#-----------------------------------------------------------------------
# list of compiler options for each operating system
#
#    $(OS) comes from outside, do not set it here
#    $(F90C) is the Fortran 90/95 compiler command
#    $COMP is the compiler vendor 
#    (up to now SUN, DEC, ABS, NEC, NAG, EPC, PGI, CRAY tested)
#    special handling for SUN, ABS, EPC implemented
#    $OPTIONS are the compiler options without library, module, and include pathes
#    $LINK_OPTS are appended to the command line at the link stage
#
# defaults
COMP      = SUN
F90C      = f90
OPTIONS   = -O -xildoff -ansi -xarch=sse3 -m64
LINK_OPTS = 
ifeq ($(OS),SUPER-UX)
         COMP      = NEC 
         F90C      = f90
         OPTIONS   = -C hopt
         LINK_OPTS = -llapack -lblas
endif
ifeq ($(OS),Linux)
         COMP      = INTEL
         F90C      = /opt/intel/bin/ifort
         OPTIONS   = -assume byterecl -O3 -m64 -u -align -qopenmp -mkl=parallel -nbs -warn errors -warn general -warn declarations -warn truncated_source -warn alignments -warn general -warn uninitialized
         LINK_OPTS = 
endif
ifeq ($(OS),SunOS)
         COMP      = SUN
         F90C      = f95
         OPTIONS   =  -w3 -xildoff -ansi -g -xs -O3 -xtarget=nehalem -xarch=sse3 -m64 -xlibmopt -xknown_lib=blas -xvector=simd -xdepend -libmil -ftrap=common,no%underflow -xlic_lib=sunperf
         LINK_OPTS = 
endif
#ifeq ($(OS),?)
#        F90C      = 
#        COMP      = 
#        OPTIONS   = 
#        LINK_OPTS = 
#endif
#----------------------------------------------------------------------
# primary directory- and filename assignments:
#
#    version independent source files in $(SRC)
#       "        "       module   "   in $(SRC)/$(MOD)      (should not be .)
#       "        "       include  "   in $(SRC)/$(INC)
#    version dependent   source   "   in $(SRC)_<v>         (can be empty)
#       "        "       module   "   in $(SRC)_<v>/$(MOD)  (should not be .)
#       "        "       include  "   in $(SRC)_<v>/$(INC)
#    objects files in $(OBJD)/$(OS) resp. $(OBJD)/$(OS)_<v>
#    executable file is $(EXCD)/$(OS)/$(MAIN)_<v>
#    module library files are constructed from $LIB and $LIBV as follows:
#        $LIB=A[X] B[Y]  and 
#        $LIBV=n m'      accesses Libraries A/$(OS)_n/libX.a and B/$(OS)_m/libY.a
#    created module library file (gmake <v> Library) is $(MLIBD)/$(OS)_<v>/lib$(MAIN).a,
#    including all .mod files
#
#    !!! DON'T DELETE EMPTY VARIABLES - BEWARE OF TRAILING BLANKS !!!
#
SRC   = source
MOD   = modules
INC   = 
OBJD  = objects
EXCD  = OS
#
# must be distinct from <v> if $(EXD) is . No % allowed for EPC compiler !
MAIN  = co2_cool
LIB   = 
MLIBD = lib
#
# set to -v if verbose behaviour of the dependency generator is desired
VERB = 
##### END #####
