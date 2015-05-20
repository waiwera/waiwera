# makefile for geothermal supermodel

# PETSc includes:
include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules 

# project directories:
SRC=src
DIST=dist
BUILD=build
TEST=test

# file extensions:
F90=.F90
OBJ=.o
SO=.so
MOD=.mod
EXE=

# shell commands:
RM=rm -f
CP=cp

# paths for compiling and linking etc.:
LIBS= -L$(HOME)/lib
LDFLAGS= $(LIBS) -lfson $(PETSC_LIB)
FMFLAGS = -J$(BUILD)
INCLS=-I$(HOME)/include
TESTLDFLAGS=$(LIBS) -lfruit -lfson $(PETSC_LIB)
TESTINCLS=$(INCLS)
TESTFMFLAGS = -J$(TEST)/$(BUILD)
INSTALL_DIR = $(HOME)/bin
FC_FLAGS =  -fPIC -Wall -Wno-unused-dummy-argument -g -O0

# main source code:
SOURCES = kinds mpi fson_mpi utils dm_utils powertable \
	thermodynamics thermodynamics_setup IAPWS IFC67 \
	timestepping rock fluid cell face mesh \
	eos eos_setup eos_w \
	initial boundary simulation
OBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(SOURCES))
PROG = supermodel
PROGEXE = $(DIST)/$(PROG)$(EXE)
DEPENDS = depends.in

# unit tests:
TESTPROG = test_all
TESTSUF = _test
TESTS = setup powertable utils IAPWS IFC67 timestepping fson_mpi \
	cell face mesh fluid rock simulation
TESTOBJS = $(patsubst %, $(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ), $(TESTS))

.DEFAULT_GOAL := $(PROGEXE)
$(PROG): $(PROGEXE)
tests: $(TEST)/$(DIST)/$(TESTPROG)$(EXE)

include $(DEPENDS)

# build rules:

# main program:
$(PROGEXE): $(BUILD)/$(PROG)$(OBJ) $(OBJS)
	$(FLINKER) $^ $(LDFLAGS) -o $@

$(BUILD)/$(PROG)$(OBJ): $(SRC)/$(PROG)$(F90) $(OBJS) $(DEPENDS)
	$(PETSC_FCOMPILE) -I$(BUILD) $(INCLS) -c $< -o $@

# main objects:
$(BUILD)/%$(OBJ): $(SRC)/%$(F90) $(DEPENDS)
	$(PETSC_FCOMPILE) $(FMFLAGS) $(INCLS) -c $< -o $@

# test program:
$(TEST)/$(DIST)/$(TESTPROG)$(EXE): $(TEST)/$(BUILD)/$(TESTPROG)$(OBJ) $(TESTOBJS) $(OBJS)
	$(FLINKER) $^ $(TESTLDFLAGS) -o $@

$(TEST)/$(BUILD)/$(TESTPROG)$(OBJ): $(TEST)/$(SRC)/$(TESTPROG)$(F90) $(TESTOBJS) $(DEPENDS)
	$(PETSC_FCOMPILE) -I$(TEST)/$(BUILD) $(TESTINCLS) -c $< -o $@

# test objects:
$(TEST)/$(BUILD)/setup$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/setup$(TESTSUF)$(F90) $(DEPENDS)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

$(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/%$(TESTSUF)$(F90) $(BUILD)/%$(OBJ) $(BUILD)/mpi$(OBJ) $(DEPENDS)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

# dependencies:
depends: $(DEPENDS)

$(DEPENDS):
	python depends.py

# documentation:
devdoc:
	ford devdoc.md

install:
	$(CP) $(DIST)/$(PROG) $(INSTALL_DIR)
clean::
	$(RM) $(BUILD)/*$(MOD) $(BUILD)/*$(OBJ)
	$(RM) $(TEST)/$(SRC)/$(TESTPROG)$(F90)
	$(RM) $(TEST)/$(BUILD)/*$(MOD) $(TEST)/$(BUILD)/*$(OBJ)
	$(RM) $(TEST)/$(DIST)/*$(TESTPROG)*
	$(RM) $(DIST)/*$(PROG)*
