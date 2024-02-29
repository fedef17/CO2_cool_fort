#======================================================================#
# DarkStar Makefile for gmake >= 3.70 , F90/95                         #
#----------------------------------------------------------------------#
# STATUS: gra 29-Nov-2007      Version 3      CREATED: gra 18.Sep 1996 #
#----------------------------------------------------------------------#
# CALL:   gmake [<v> [Library] [LIBV='<v1> <v2> ...']], where <v>,<v1> #
#         is any string (except "_", source__ is always removed !) of  #
#         existing source_<v> directories; defaults to first one found #
#         or "_", if none exists (then all source files have to be in  #
#         "source").The keyword "Library" generates a library of all   #
#         modules, the LIBV='<v1> <v2>' specifies linking of libs with #
#         different versions.                                          #
#----------------------------------------------------------------------#
# FILES:  Generates and uses .<MAIN>.<v>.dep, .<OS>.Version_<v>        #
#         The user-specific part has to be in Makefile.inc (included)  #
#----------------------------------------------------------------------#
# USES:   Dependency - generator Make_Depp                             #
#----------------------------------------------------------------------#
# NOTES AND RESTRICTIONS:                                              #
#         This Makefile is able to manage different versions of a      #
#         source package, with version independent files of the source #
#         in "source" and version dependent parts in "source_<v>",     #
#         which may have complex dependencies on each other.           #
#         Version dependent files are recognized automatically just by #
#         being in the source_<v>/ directory. A file of the same name  #
#         can be left in source/, thus different versions can have     #
#         different version-dependent files. Defaults:                 #
#         Module files have to reside in source/modules/*.f90 or in    #
#         source_<v>/modules/*.f90 (user module files), include files  #
#         in source/includes or source[_<v>]/includes, library modules #
#         in <dir>/<OS>_<v>/libX.a, other source files in              #
#         source[_v]/*.f90 (see declarations below for modifications). #
#         For module and include files the source_<v> directory is     #
#         searched first, then the source directory. Library include   #
#         files are not supported.                                     #    
#         Library module files are not part of the make process,       #
#         they are considered to be up to date.                        #
#         The executable name must be distinct from the version name,  #
#         or it must be in a subdirectory (see EXD ).                  #
#         Module filenames and USE statements must be case-identical.  #
#         Module filenames and other filenames must be distinct.       #
#         There may only be one use statement per line.                #
#         Preprocessor #include directives are not handled.            #
#         The target 'Library' generates a library of all module object#
#         files in $MLIBD/$OS_<v>/lib$(MAIN).a and places all compiler #
#         dependent module interface files (.mod, .M, .d) therein.     #
#         With >make CLEANUP<, all objects and .* files are deleted.   #
#         This Makefile supports working on multiple architectures     #
#         (e.g., NFS-mounted directories on heterogeneous cluster).    #
#----------------------------------------------------------------------#
# AUTHOR: Dr.U.Grabowski,Institut f"ur Meteorologie und Klimaforschung #
#         Forschungszentrum Karlsruhe, Germany.                        #
#         Any suggestions to improve this code are welcome.            #
#         Send mail to Udo.Grabowski@imk.fzk.de                        #
#----------------------------------------------------------------------#
# LEGAL:  This code can be used and exchanged freely except for        #
#         commercial purposes.                                         #
#         There is no guarantee that this code is bug-free.            #
#----------------------------------------------------------------------#
# ACKNOWLEDGEMENTS:                                                    #
#         This code has been developed as part of the DarkStar Project #
#         at the Astronomical Institute Basel, Switzerland,            #
#         supported by the Swiss National Science Foundation,          #
#         project no. 2100-045568                                      #
#----------------------------------------------------------------------# 
# MAINTENANCE HISTORY:                                                 #
# gra 29-Nov-07:  workaround for SUN missing __xpg6 problem            #
# gra 05-Jul-02:  corrected bug for absolute path names in libraries   #
# gra 17-Jun-02:  corrected module flags for more than one library     #
# gra 13-Jun-02:  corrected deletion of directories in level 0 pass    #
# gra 27-May-02:  remove executable only on build;leave other versions #
# gra 23-May-02:  Version 3:generated lib has version tag in directory;#
#                 module libraries attachable by version number        #
# gra 29-Dec-01:  make CLEANUP safer                                   #
# gra 11-Apr-01:  restricted CLEANUP to managed executable             #
# gra 02-Jul-00:  added version to name of executable                  #
# gra 20-Mar-00:  added $LINK_OPTS to linkage stage                    #
# gra 20-Oct-99:  adapted to Cray T3E f90 compiler                     #
# gra 18-Jun-99:  changes for EPC compiler                             # 
# gra 17-Jun-99:  can now build module library; safer CLEANUP; Usage   #
# gra 23-Apr-99:  forced at least a linking process for consistency    #
# gra 25-Feb-99:  changed <v> default to "_"; temporary creation of    #
#                 "source__" to enable usage with single source dir.   #
# gra 23-Feb-99:  added support for multiple architecture libraries;   #
#                 spawned user dependent part to Makefile.inc          #
# gra 12-Feb-99:  added support for multiple architectures             #
# gra 30-Nov-98:  changed intermediate target name to $@__             #
# gra 11-Nov-98:  corrected bug in version control                     #
# gra 22-Oct-98:  added support for EPC compiler, script epc_work.sh   #
# gra 06-Jun-98:  new target CLEANUP                                   #
# gra 18-Apr-98:  explicit 'cd' to module directories                  #
# gra 28-Mar-98:  does not create unused objects_<v> anymore;          #
#                 default compilation of first directory found         #
# gra 25-Mar-98:  corrected for missing source_<v> directories         #
# gra 24-Mar-98:  optional source file lists; any version names <v>    #
# gra  4-Mar-98:  changed handling of include files                    #
# gra 18-Feb-98:  corrected bug in objects list                        #
# gra 12-Feb-98:  corrected include directory problem for SUN          #
# gra  6-Feb-98:  automatic recognition of version dependent files     #
# gra 28-Jan-98:  automatic library name, removed hardwired directory  #
#                 and >EMPTY< target; dummy .old no longer required    #
# gra 22-Jan-98:  >phony< for 1-9,test targets                         #
# gra 17-Jan-98:  added >test< target; modified dependency for main    #
# gra 11-Jan-98:  unified targets, some simplifications                #
# gra 08-Jan-98:  completely restructured, gmake conformant            #
# gra 06-Jan-98:  corrected include dependency handling                #
# gra 17-Aug-97:  corrected version dependency handling                #
# gra 11-Aug-97:  can now deal with module files in source_<v>/modules #
# gra 25-Apr-97:  changed for updated MAKE_DEPP                        #
# gra 13-Apr-97:  now handles up to 9 different source versions        #
# gra 21-Oct-96:  added library dependence on rlib_F90 routines        #
#======================================================================#
# Makefile version ;operating system indicator
MV := 3
OS := $(shell uname -s)
export OS
OSS = $(strip $(OS))
#    $(MAK) should be set to <Makefile-name> (this file)
#    $(MAKI) should be set to  <include-makefile-name>
#    $(MAKD) is the name of the dependency generator
#
MAK  = Makefile
MAKI = Makefile.inc
MAKD = Make_Depp
######################################
# include user specific part from file
include $(MAKI)
########################## END OF USER SPECIFIC PART #########################
############### USUALLY THERE IS NO NEED TO EDIT BELOW THIS LINE #############
################ !!! CAUTION: THIS IS A RECURSIVE MAKEFILE !!! ###############
################ VARIABLES BELOW ARE NOT INTENDED TO BE CHANGED ##############
############################ (EXCEPT FOR RARE CASES) #########################
############ !!! SELF-DESTRUCTION IS INITIATED WITHIN 5 SECONDS !!! ##########
export
.SUFFIXES:
DIR   = $(shell echo \$(PWD))
CHD   = \cd
ECHO  = /bin/echo
LINE  = $(ECHO) '==========================================================================='
USAGE =	'Usage: gmake <version> [Library]'
DEP   = .$(MAIN).$(OS).dep
#
# allow multiple architectures
#
ifneq ($(strip $(OBJD)),)
   OBD = $(OBJD)/$(OS)
else
   OBD = $(OS)
endif
ifneq ($(strip $(EXCD)),)
   EXD = $(EXCD)/$(OS)
else
   EXD = $(OS)
endif
VFIL = .$(OS).Version
export VFIL
#
# generate version list from existing source_<v> directories;
# take first found directory as default, "_" if none was found
#
VS = $(wildcard $(SRC)_*)
VERSIONS = $(VS:$(SRC)_%=%)
ifeq ($(strip $(VERSIONS)),)
     VERSIONS = "_"
endif
V = $(firstword $(VERSIONS))
#
# construct library names; complete compiler options
#
NUMLIB =
NUMV   =
ifeq ($(strip $(LIB)),)
        RLIB =
else
	NUMLIB   = $(words $(LIB))
	NUMV     = $(words $(LIBV))
	lbrack  := [
	rbrack  := ]
	OSP      = $(addprefix $(OS)_,$(LIBV))
	LIBSTR   = $(subst $(rbrack),,$(LIB))
	GIVEDIR  = $(firstword $(subst $(lbrack), ,$1))
	LIBDIRS  = $(foreach lib,$(LIBSTR),$(call GIVEDIR,$(lib))/)
	SECOND   = $(word 2,$1)
	GIVELIB  = $(call SECOND,$(subst $(lbrack), ,$1))
	LIBF     = $(foreach lib,$(LIBSTR),/lib$(call GIVELIB,$(lib)).a)
	RLIBWO   = $(join $(LIBDIRS),$(OSP))
	LIBR     = $(join $(RLIBWO), $(LIBF))
	RELPATH  = $(subst /,,$(firstword $(subst /,/ ,$1)))
        RLIB     = $(foreach lib,$(RLIBWO),$(if $(call RELPATH,$(lib)),$(DIR)/$(lib),$(lib)))
endif
LR      = $(LIBF:/lib%.a=%)
FFLAGS  = $(OPTIONS)
ifeq ($(strip $(LIBF)),)
	LIBR   =
	LIBC   =
else
        LIBL   = $(addprefix -L,$(RLIB))
	LIBl   = $(addprefix :-l,$(LR))
        LIBC   = $(subst :, ,$(join $(LIBL),$(LIBl)))
endif
ifeq ($(strip $(MLIBD)),)
   MLBD =
   MLIB =
else
   MLBD = $(MLIBD)/$(OS)_$V
   MLIB = $(MLBD)/lib$(MAIN).a
endif
#
# construct main executable name,
# include-, module source-, and module object directories
#
EXC  = $(EXD)/$(MAIN)_$V
ifneq ($(strip $(SRC)),)
   SRV = $(SRC)_$V
else
   SRV =
endif
ifneq ($(strip $(INC)),)
   SIN = $(SRC)/$(INC)
   SIV = $(SRV)/$(INC)
else
   SIN =
   SIV =
endif
ifneq ($(strip $(MOD)),)
   SMD = $(SRC)/$(MOD)
   SMV = $(SRV)/$(MOD)
else
   SMD =
   SMV =
endif
ifneq ($(strip $(OBD)),)
   OBV = $(OBD)_$V
else
   OBV =
endif
ifneq ($(strip $(MOD)),)
   OMD = $(OBD)/$(MOD)
   OMV = $(OBV)/$(MOD)
else
   OMD = $(OBD)
   OMV = $(OBV)
endif
MAKL = $(MAKI) $(MAK) $(MAKD)
ifneq ($(strip $(INC)),)
        FFLAGS +=  -I$(DIR)/$(SIV)
endif
ifneq ($(strip $(INC)),)
        FFLAGS +=  -I$(DIR)/$(SIN)
endif
ifneq ($(strip $(OMV)),)
        FFLAGS +=  -I$(DIR)/$(OMV)
endif
ifneq ($(strip $(OMD)),)
        FFLAGS +=  -I$(DIR)/$(OMD)
endif

ifneq ($(strip $(RLIB)),)
        FFLAGS += $(foreach lib,$(RLIB),-I$(lib))
endif
FFLAGS +=  -I$(DIR) -I.
#
# compiler dependent module search-directory options
# (only set when library string is not empty)
#
# CRAY
ifeq ($(COMP),CRAY)
        FFLAGS +=  -p $(DIR)/$(OMV) -p $(DIR)/$(OMD)
        ifneq ($(strip $(RLIB)),)
              FFLAGS += $(foreach lib,$(RLIB),-p $(lib))
        endif
endif
# ABSOFT
ifeq ($(COMP),ABS)
        FFLAGS += -YEXT_NAMES=LCS -p $(DIR)/$(OMV) -p $(DIR)/$(OMD)
        ifneq ($(strip $(RLIB)),)
              FFLAGS += $(foreach lib,$(RLIB),-p $(lib))
        endif
endif
# SUN
ifeq ($(COMP),SUN)
        FFLAGS += -M$(DIR)/$(OMV) -M$(DIR)/$(OMD)
        ifneq ($(strip $(RLIB)),)
	      FFLAGS += $(foreach lib,$(RLIB),-M$(lib))
        endif
endif
# EPC
ifeq ($(COMP),EPC)
        FFLAGI = -cl,$(MAIN).pcl
        FFLAGV = -cl,$(MAIN)\%$V.pcl
else
        FFLAGI = 
        FFLAGV = 
endif
#
# generate object and source file lists with directory parts
# version dependent sources are recognized automatically 
#
ifeq ($(strip $(SOURCE_FILES)),)
        SRCL     = $(wildcard $(SRC)/*.f90)
        SRCLIST  = $(SRCL:$(SRC)/%.f90=%)
else
        SRCLIST  = $(SOURCE_FILES)
endif
ifeq ($(strip $(SOURCE_V_FILES)),)
        SRC_USER = $(wildcard $(SRV)/*.f90)
        SRCUSER  = $(SRC_USER:$(SRV)/%.f90=%)
else
        SRCUSER  = $(SOURCE_V_FILES)
        SRC_USER = $(SRCUSER:%=$(SRV)/%.f90)
endif
SOURCE   = $(filter-out  $(SRCUSER), $(SRCLIST))
SOURCES  = $(SOURCE:%=$(SRC)/%.f90) $(SRC_USER)
OBJECTS  = $(SOURCE:%=$(OBD)/%.o)   $(SRCUSER:%=$(OBV)/%.o)
###!!! Workaround for Sun missing __xpg6 problem
ifeq ($(COMP),SUN)
     ifeq ($(OS),Linux)
        OBJECTS += $(OBV)/xpg6.o
     endif
endif
#
# ensure that a consistent executable is produced; create and delete temporary dirs
#
GOALS = $(firstword $(MAKECMDGOALS))
.PHONY: NIL $(VERSIONS) $V__ $V CLEANUP Usage Library_
.DELETE_ON_ERROR:
NIL:
	-@ if [ "$V" = "_" ] ; then \
	      if [ ! -d $(SRC)_$V ] ;then mkdir -p $(SRC)_$V/$(MOD) ; fi;\
	  fi;
	   $(MAKE) $V;
	-@ if [ "$V" = "_" ]; then rm -rf $(SRC)_$V $(OBD)_$V ; fi;
#
# create nonexistent directories ; ensure version consistency; set variable
# $V for different source versions ; generate dependencies and make version
#
$(VERSIONS):
	 @ if [ "$(OSS)" = ""  -o  "$(OSS)" = "."  ] ;then echo "### ERROR: NO OS ###"; false; fi;
	 @ if [ "$(MIFVN)" != "$(MV)" ] ;then echo "### Please use Version $(MV) Makefile.inc ! ###"; false; fi;
	 @ if [ "$(NUMLIB)" != "$(NUMV)" ] ;then echo "### Error: Different number of entries in LIB and LIBV ### "; false; fi;
	-@ if [ ! -d $(EXD)    ] ;then mkdir -p $(EXD)           2>/dev/null ;fi;
	-@ mkdir -p $(OBD) $(OBD)/$(MOD) 2>/dev/null ;
	-@ if [   -d $(SRC)_$@ ] ;then mkdir -p $(OBD)_$@/$(MOD) 2>/dev/null ;fi;
	-@ rm -rf `\ls $(VFIL)_* 2>/dev/null|\grep -v $(VFIL)_$@` $(EXD)/$(MAIN)_$@ $(MLIBD)/$(OS)_$@;
	-@ echo "Operating System = ***"$(OS)"***" ;
	 @ $(MAKE) V=$@ $@__ ;
$V__:
	 @ ./Make_Depp $(VERB) $(@:%__=%) ;
	 @ if [ "$(COMP)" = "EPC" ] ; then ./epc_work.sh ; fi ;
	 @ echo "======================= >> Making Version $(@:%__=%) << ======================="; 
	 @ $(MAKE) V=$(@:%__=%) $(EXC)
#
# main entry: generate executable
#
$(EXC): $(VFIL)_$V $(OBJECTS) $(LIBR) $(MAKL) $(DEP)
	$(F90C) $(FFLAGS)  $(MODULES) $(OBJECTS) -o $(EXC) $(LIBC) $(LINK_OPTS) ;
	-@$(LINE) 
#
# Here the module and include-file dependencies are included as macro
# definitions <file>_dep=... and <file>_inc=... from file .DarkStar.dep
#
include $(DEP)
#
# object files entry: depends on modules, include and source files
#
$(OBD)/%.o: $(SRC)/%.f90 $(LIBR) $(MAKL)
	$(CHD) $(@D) ; $(F90C) $(FFLAGS) $(FFLAGI) -c  $(DIR)/$(SRC)/$*.f90 ;
	-@$(CHD) $(DIR) ; $(LINE)
#
# user object files entry
#
$(OBV)/%.o: $(SRV)/%.f90 $(LIBR) $(MAKL)
	$(CHD) $(@D) ; $(F90C) $(FFLAGS) $(FFLAGV) -c  $(DIR)/$(SRV)/$*.f90 ;
	-@$(CHD) $(DIR) ; $(LINE)
#
# workaround for stupid Sun bug
#
$(OBV)/xpg6.c:
	echo unsigned int __xpg6 = 0xFFFF0000\; > $@
$(OBV)/xpg6.o: $(OBV)/xpg6.c
	suncc -c $< -o $@
#
# handle version dependencies of basic source and module files
#
$(VFIL)_$V : $(MAKL)
	@touch $(VFIL)_$V
#
# generate dummy dependency file, if not existent
#
$(DEP) :
	touch $(DEP)
#
# module library creation; copy module interface files
#
Library:
	@if [ $(MAKELEVEL) = 0 ]; then \
	   $(MAKE) V=$(GOALS) Library; \
	else \
	   if [ "$V" = "$@" ]; then \
	      echo $(USAGE); exit; \
	   fi; \
	   echo "============================== Making Library ============================="; \
	   if [ ! -d $(MLBD) ]  ;then \
	      mkdir -p $(MLBD) 2>/dev/null ; \
	   fi; \
	   cp $(OMD)/*.[!o]* $(MLBD) 2>/dev/null; cp $(OMV)/*.[!o]* $(MLBD) 2>/dev/null ; \
	   echo "ar cr $(MLIB) $(MODULES)" && ar cr $(MLIB) $(MODULES); \
	   echo 'ranlib' $(MLIB) && ranlib $(MLIB) 2>/dev/null || exit 0; \
	fi;
#
# delete all object, module, version, and dependency files
#
CLEANUP:
	@ if [ "$OSS" = "" -o "$(OSS)" = "." ] ;then echo "*** ERROR: WRONG OS ! ***"; exit 1; fi;
	@$(ECHO) "CAUTION: Removing " $(OBJD)/$(OS)\* "and other scrap ...";
	@ if [ "$(OBJD)/$(OS)" = "*" -o "$(MLBD)" = "*" -o "$(DEP)" = "*" ] ; then \
	      echo "Error: Will not remove '*' ! Check your ${MAKI} !"; exit 1; fi;
	@sleep 2
	@$(ECHO) "Self-Destruction initiated:   \\c";
	@for i in  3 2 1 Let\ there\ be\ light\ ! ; \
	 do $(ECHO) \\b\\b$$i' '\\c ; done; echo ; \
	 rm -rf $(OBJD)/$(OS) $(OBJD)/$(OS)_*  $(SRC)__ ; rm -f $(DEP) $(VFIL)_* ;
Usage:
	@echo $(USAGE);
