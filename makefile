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
TESTLDFLAGS=$(LIBS) -lfruit $(PETSC_LIB)
TESTINCLS=$(INCLS)
TESTFMFLAGS = -J$(TEST)/$(BUILD)
INSTALL_DIR = $(HOME)/bin

# modules used by all other modules:
ESSENTIAL = kinds
ESSENTIAL_OBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(ESSENTIAL))

# main source code:
SOURCES = mpi fson_mpi powertable thermodynamics IAPWS IFC67 timestepping mesh simulation eos eos_w
OBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(SOURCES))
ALLOBJS = $(ESSENTIAL_OBJS) $(OBJS)
PROG = supermodel

# unit tests:
TESTPROG = test_all
TESTSUF = _test
# test modules that need setup/teardown:
SETUPTESTS = IAPWS IFC67 timestepping
NONSETUPTESTS = powertable
TESTS = setup $(SETUPTESTS) $(NONSETUPTESTS)
SETUPOBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(SETUPTESTS))
TESTOBJS = $(patsubst %, $(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ), $(TESTS))

$(PROG): $(DIST)/$(PROG)$(EXE)
tests: $(TEST)/$(DIST)/$(TESTPROG)$(EXE)

# general dependency rules:
$(OBJS) $(TESTOBJS): $(ESSENTIAL_OBJS)

# specific dependency rules:
$(TEST)/$(BUILD)/setup$(TESTSUF)$(OBJ): $(SETUPOBJS)
$(BUILD)/IAPWS$(OBJ): $(BUILD)/thermodynamics$(OBJ) $(BUILD)/powertable$(OBJ)
$(BUILD)/IFC67$(OBJ): $(BUILD)/thermodynamics$(OBJ) $(BUILD)/powertable$(OBJ)
$(BUILD)/simulation$(OBJ): $(BUILD)/mpi$(OBJ) $(BUILD)/timestepping$(OBJ) \
	$(BUILD)/thermodynamics$(OBJ) $(BUILD)/IAPWS$(OBJ) \
	$(BUILD)/IFC67$(OBJ) $(BUILD)/eos$(OBJ)	$(BUILD)/eos_w$(OBJ) \
	$(BUILD)/mesh$(OBJ) $(BUILD)/fson_mpi$(OBJ)
$(BUILD)/eos$(OBJ): $(BUILD)/thermodynamics$(OBJ)
$(BUILD)/eos_w$(OBJ): $(BUILD)/thermodynamics$(OBJ) $(BUILD)/eos$(OBJ)
$(BUILD)/fson_mpi$(OBJ): $(BUILD)/mpi$(OBJ)
$(BUILD)/$(PROG)$(OBJ): $(BUILD)/mpi$(OBJ) $(BUILD)/simulation$(OBJ)

# build rules:

# main program:
$(DIST)/$(PROG)$(EXE): $(BUILD)/$(PROG)$(OBJ) $(ALLOBJS)
	$(FLINKER) $^ $(LDFLAGS) -o $@

$(BUILD)/$(PROG)$(OBJ): $(SRC)/$(PROG)$(F90) $(ALLOBJS)
	$(PETSC_FCOMPILE) -I$(BUILD) $(INCLS) -c $< -o $@

# main objects:
$(BUILD)/%$(OBJ): $(SRC)/%$(F90)
	$(PETSC_FCOMPILE) $(FMFLAGS) $(INCLS) -c $< -o $@

# test program:
$(TEST)/$(DIST)/$(TESTPROG)$(EXE): $(TEST)/$(BUILD)/$(TESTPROG)$(OBJ) $(TESTOBJS) $(ALLOBJS)
	$(FLINKER) $(TESTLDFLAGS) $^ -o $@

$(TEST)/$(BUILD)/$(TESTPROG)$(OBJ): $(TEST)/$(SRC)/$(TESTPROG)$(F90) $(TESTOBJS)
	$(PETSC_FCOMPILE) -I$(TEST)/$(BUILD) $(TESTINCLS) -c $< -o $@

# test objects:
$(TEST)/$(BUILD)/setup$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/setup$(TESTSUF)$(F90)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

$(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/%$(TESTSUF)$(F90) $(BUILD)/%$(OBJ)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

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
