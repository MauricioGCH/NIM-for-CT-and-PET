#
# Generic Makefile for any C++ include/src project
#
#      --> all main programs in the root directory
#      --> a directory named $(INCDIR) containing header files
#      --> a directory named $(SRCDIR) containing source files
#
# 2012 Simon Stute : gate.stute@gmail.com
##############################################

# Get the build date
######################
BUILD_DATE := "\"`date`\""
DATEFLAGS := -DBUILD_DATE=$(BUILD_DATE)

# Directories settings
########################
INCDIR := include
SRCDIR := src
TMPDIR := tmp
BINDIR := bin

# Extensions settings
#######################
SRCEXT := .cc
DEPEXT := .d

# Makefile implicit variables
###############################
ifeq ($(CASTOR_MINGW64), 1)
  # Change compiler
  CXX = x86_64-w64-mingw32-g++
endif
ifndef CXX
CXX = g++
endif
CXXFLAGS = $(INCLUDE)
LDFLAGS = $(INCLUDE)

ifeq ($(ARCH),x86_64)
CXXFLAGS += -D_64
endif

ifeq ($(CASTOR_MINGW64), 1)
  CXXFLAGS += -DWINDOWS
  # Add static linking flags to the linker
  LDFLAGS += -static-libgcc -static-libstdc++ -static -lpthread
endif

# C11 options
####################
#CXXFLAGS += -std=c++0x

# Optimization options
####################
#CXXFLAGS += -O3 -funroll-all-loops -march=native -pipe
CXXFLAGS += -O3 -funroll-all-loops
# OpenMP multi-threading
CXXFLAGS += -fopenmp -DOMP_MODE
# Input and output files
##########################
SRC  = $(wildcard $(SRCDIR)/*$(SRCEXT))
MAIN = $(wildcard *$(SRCEXT))
EXE  = $(patsubst %$(SRCEXT),$(BINDIR)/%.exe,$(MAIN))
OBJ  = $(patsubst $(SRCDIR)/%$(SRCEXT),$(TMPDIR)/%.o,$(SRC))
OBJ += $(patsubst %$(SRCEXT),$(TMPDIR)/%.ox,$(MAIN))
DEP  = $(patsubst $(SRCDIR)/%$(SRCEXT),$(TMPDIR)/%$(DEPEXT),$(SRC))
DEP += $(patsubst %$(SRCEXT),$(TMPDIR)/%$(DEPEXT),$(MAIN))
INCLUDE = -I$(INCDIR)/

# Architecture
#########################
ARCH  := $(shell uname -p)

#######################################################################
##                               Rules                               ##
#######################################################################

################
##  all Part  ##
################
all : $(TMPDIR)/exe.last

#################
##  EXEC Part  ##
#################
$(BINDIR)/%.exe : $(TMPDIR)/%.ox $(TMPDIR)/obj.last
	@if [ ! -d $(BINDIR) ] ; then mkdir -p $(BINDIR) ; fi
	@echo "Linking" $(subst $(BINDIR)/,,$@) "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $(TMPDIR)/*.o $<

################
##  OBJ Part  ##
################
# Compiling $(SRCDIR)
$(TMPDIR)/%.o : $(SRCDIR)/%$(SRCEXT) $(TMPDIR)/%$(DEPEXT)
	@if [ ! -d $(TMPDIR) ] ; then mkdir -p $(TMPDIR) ; fi
	@echo "Compiling" $(subst $(SRCDIR)/,,$<) "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -c -o $@ $<

# Compiling $(SRCEXT)
$(TMPDIR)/%.ox : %$(SRCEXT) $(TMPDIR)/%$(DEPEXT)
	@if [ ! -d $(TMPDIR) ] ; then mkdir -p $(TMPDIR) ; fi
	@echo "Compiling" $(subst $(SRCDIR)/,,$<) "..."
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(DATEFLAGS) -c -o $@ $<

######################
##  .PHONY targets  ##
######################
.PHONY : all obj clean tilde

obj : $(TMPDIR)/obj.last

$(TMPDIR)/obj.last : $(OBJ)
	@touch $@

$(TMPDIR)/exe.last : $(EXE)
	@touch $@

################
##  DEP Part  ##
################
# Making dependencies for $(EXEC)$(SRCEXT)
$(TMPDIR)/%$(DEPEXT) : %$(SRCEXT)
	@if [ ! -d $(TMPDIR) ] ; then mkdir -p $(TMPDIR) ; fi
	@echo "Making dependencies for" $(subst $(SRCDIR)/,,$<) "..."
	@set -e; rm -f $@; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
ifneq ($(DEP),)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),tilde)
-include $(DEP)
endif
endif
endif

# Making dependencies for $(SRCDIR)
$(TMPDIR)/%$(DEPEXT) : $(SRCDIR)/%$(SRCEXT)
	@if [ ! -d $(TMPDIR) ] ; then mkdir -p $(TMPDIR) ; fi
	@echo "Making dependencies for" $(subst $(SRCDIR)/,,$<) "..."
	@set -e; rm -f $@; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
ifneq ($(DEP),)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),tilde)
-include $(DEP)
endif
endif
endif

##################
##  clean Part  ##
##################
clean ::
	@echo ""
	@echo "Cleaning up ..."
#	@rm -rf $(BINDIR)/
	@rm -rf $(TMPDIR)/
	@echo "Done."
	@echo ""

#################
##  tild Part  ##
#################
tilde ::
	@echo ""
	@echo "Removing all ~ files ..."
	@rm -f *~ */*~ */*/*~
	@echo "Done."
	@echo ""

