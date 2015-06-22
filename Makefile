# Generic Makefile for the compilation of a continuation problem.
# For every problem defined two programs will be created, one being the
# standalone version, the other the continuation version.
# These programs will be called "<model>" and "<model>_trc", respectively.
#
# Possible targets to make:
#
#       make <model>
#       make <model_trc
#       make all                        (makes all possible target programs in current directory)
#       make both                       (makes first <model> and <model>_trc)
#       make clean                      (cleans up all target programs)
#       make cleancvf                   (cleans up all CVF files)
#       make allclean                   (combines previous two)
#
#       make PROBLEM=model1 both        (makes <model1> and <model1>_trc)
#       make PROBLEM=model1 cleanone    (cleans up <model1> and <model1>_trc)
#
#==============================================================================
# Required tayloring
#==============================================================================

# Specify whether the problems require numerical integration of ODEs module (0: No; 1: Yes)

WITHINTEGRATION = 1

# Specify whether the problems involve structured populations (0: No; 1: Yes)

WITHSTRUCTURE = 1

# Specify whether a header file should be included (0: No; 1: Yes)

WITHHEADER = 1

# Specify whether a CVF file should be used for input of parameters and settings (0: No; 1: Yes)

WITHCVFFILE = 1

# Specify where the modules are located. For example:
MODDIR = ./Cont_modules

# MODDIR  = /Users/andre/programs/Cont_modules

#==============================================================================
# Advanced tayloring
# The following settings might need adaptation for specific uses
#==============================================================================

# Determining the names of all target programs and CVF files in the current directory
# All .c files in the current directory are considered target programs

ALLCFILES = $(wildcard *.c)
ALLPROBLEMS = $(patsubst %.c,%,$(ALLCFILES))
ALLCVF = $(wildcard *.cvf)

# Alternatively, if not all .c files in the current directory are target programs
# or not all CVF file should be clean by 'make clean' and 'make cleancvf'
# specify the targets here:

# ALLPROBLEMS = <model1> <model2>
# ALLCVF = <run1.cvf> <run2.cvf>

#==============================================================================
# Find the ebtclean program
#==============================================================================

EBTCLEAN = $(shell which ebtclean)
ifeq ("$(EBTCLEAN)", "")
SYSTEM = $(shell uname)
ifeq ("$(SYSTEM)", "Darwin")                    # This is Mac OS specific
EBTCLEAN = $(patsubst %/Resources,%/MacOs/ebtclean, $(EBTPATH))
endif
endif

#==============================================================================
#  Compilation settings for GCC compiler
#==============================================================================

CC	 = gcc
WARN     = -Wpointer-arith -Wcast-qual -Wcast-align

# Debugging settings
ifeq ("$(DEBUG)", "debug")
CCFLAGS  = -ggdb -g3 -DDEBUG=1 $(INCLUDES) $(WARN) $(SPECDEFS)
TMPDIR   = .
else
CCFLAGS  = -O $(INCLUDES) $(WARN) $(SPECDEFS)
TMPDIR   = /tmp
endif

# Generic compilation and linking flags
INCLUDES = -I. -I$(MODDIR)

CFLAGS	 = $(CCFLAGS) $(PROBLEM_CFLAGS)
LDFLAGS  = -fno-common $(EXTRA_LDFLAGS)
LIBS	 = -llapack -lcblas -latlas -lm
BASEMODS = curve io nleq_res utils

#==============================================================================
# About the modules
#==============================================================================

PROBLEMUC=$(shell echo $(PROBLEM) | tr "[:lower:]" "[:upper:]")

ifeq ("$(WITHHEADER)", "1")
PROBLEMH = $(PROBLEM).h
PROBLEM_CFLAGS = -DSYSTEM_TYPE=SYSTEM_$(PROBLEMUC) -DPROBLEMHEADER="<$(PROBLEMH)>" -DCVFFILE=$(WITHCVFFILE)
else
PROBLEM_CFLAGS = -DSYSTEM_TYPE=SYSTEM_$(PROBLEMUC) -DCVFFILE=$(WITHCVFFILE)
endif

MODULES = $(BASEMODS)
ifeq ("$(WITHSTRUCTURE)", "1")
ifeq ("$(SUNDIALS)", "1")
MODULES = cohort cvode $(BASEMODS)
LIBS	= -lsundials_nvecserial -lsundials_cvode -llapack -lcblas -latlas -lm
else
MODULES = cohort dopri5 $(BASEMODS)
endif
else
ifeq ("$(WITHINTEGRATION)", "1")
ifeq ("$(SUNDIALS)", "1")
MODULES = cvode $(BASEMODS)
LIBS	= -lsundials_nvecserial -lsundials_cvode -llapack -lcblas -latlas -lm
else
MODULES = dopri5 $(BASEMODS)
endif
endif
endif

DEBUGOBJS = $(patsubst %,./%.o,$(MODULES))
MODOBJS = $(patsubst %,$(TMPDIR)/%.o,$(MODULES))
MODLIB  = modules_$(PROBLEM).lib.o

#==============================================================================
# The following targets are valid if PROBLEM is not defined
#==============================================================================

ifeq ("$(PROBLEM)", "")
all:
	for I in $(ALLPROBLEMS) ; do $(MAKE) PROBLEM=$${I} both ; done

both:
	$(MAKE) PROBLEM=$(word 1,$(ALLPROBLEMS)) both

allclean: clean cleancvf

clean: 
	@echo "Cleaning up all programs: "
	@for I in $(ALLPROBLEMS) ; do $(MAKE) PROBLEM=$${I} cleanone ; done

cleancvf:
	@echo "Cleaning up all output files: "
	@for I in $(patsubst %.cvf, %, $(ALLCVF)) ; do echo "Cleaning up $${I}..."; $(EBTCLEAN) -f $${I} >/dev/null ; done
	@for I in $(ALLPROBLEMS) ; do rm -f $${I}_trc.err $${I}_trc.out ; done

# Re-invoke make but now with PROBLEM defined and the same target

$(ALLPROBLEMS) $(patsubst %,%_trc,$(ALLPROBLEMS))::
	$(MAKE) PROBLEM=$(subst _trc,,$@) $@

#==============================================================================
# The following targets are valid if PROBLEM is defined
#==============================================================================

else
PROBLEM_TRC = $(PROBLEM)_trc

#==============================================================================

# The dependencies of the executables

both: $(PROBLEM) $(PROBLEM_TRC)

$(PROBLEM):     SPECDEFS = -DCONTINUATION=0
$(PROBLEM_TRC): SPECDEFS = -DCONTINUATION=1

$(PROBLEM_TRC): $(PROBLEM_TRC).o $(MODLIB)
	gcc $(LDFLAGS) -o $@ $(PROBLEM_TRC).o $(MODLIB) $(LIBS)

$(PROBLEM): $(PROBLEM).o $(MODLIB)
	gcc $(LDFLAGS) -o $@ $(PROBLEM).o $(MODLIB) $(LIBS)

#==============================================================================

# The dependencies of the problem-specific object files

$(PROBLEM).o: $(PROBLEM).c $(PROBLEMH)
	$(COMPILE.c) -o $@ $(PROBLEM).c

$(PROBLEM_TRC).o: $(PROBLEM).c $(PROBLEMH)
	$(COMPILE.c) -o $@ $(PROBLEM).c

#==============================================================================

# The dependencies of the module library file

ifeq ("$(DEBUG)", "debug")
$(MODLIB): $(PROBLEMH)
	for I in $(MODULES) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(MODDIR)/$${I}.c ; done
	ld -r -o $@ $(MODOBJS)
else
$(MODLIB): $(PROBLEMH)
	for I in $(MODULES) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(MODDIR)/$${I}.c ; done
	ld -r -o $@ $(MODOBJS)
	rm $(MODOBJS)
endif

#==============================================================================

# The dependencies of some additional targets

cleanone: 
	@echo "Cleaning up $(PROBLEM)...."
	@rm -f $(PROBLEM)     $(PROBLEM).o
	@rm -f $(PROBLEM_TRC) $(PROBLEM_TRC).o
	@rm -f $(MODLIB) $(DEBUGOBJS)

endif

#==============================================================================
#==============================================================================
