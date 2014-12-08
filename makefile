# makefile for geothermal supermodel

# PETSc includes:
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

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

# compiling and linking paths:
LDFLAGS= $(PETSC_LIB)
FMFLAGS = -J$(BUILD)
TESTLDFLAGS=-L$(HOME)/lib -lfruit $(PETSC_LIB)
TESTINCLS=-I$(HOME)/include
TESTFMFLAGS = -J$(TEST)/$(BUILD)

# modules used by all other modules:
ESSENTIAL = kinds
ESSENTIAL_OBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(ESSENTIAL))

# main source code:
SOURCES = powertable thermodynamics IAPWS IFC67 timestepping
OBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(SOURCES))
ALLOBJS = $(ESSENTIAL_OBJS) $(OBJS)

# unit tests:
TESTPROG = test_all
TESTSUF = _test
# test modules that need setup/teardown:
SETUPTESTS = IAPWS IFC67 timestepping
NONSETUPTESTS = powertable
TESTS = setup $(SETUPTESTS) $(NONSETUPTESTS)
SETUPOBJS = $(patsubst %, $(BUILD)/%$(OBJ), $(SETUPTESTS))
TESTOBJS = $(patsubst %, $(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ), $(TESTS))

tests: $(TEST)/$(DIST)/$(TESTPROG)$(EXE)

# general dependency rules:
$(OBJS) $(TESTOBJS): $(ESSENTIAL_OBJS)

# specific dependency rules:
$(TEST)/$(BUILD)/setup$(TESTSUF)$(OBJ): $(SETUPOBJS)
$(BUILD)/IAPWS$(OBJ): $(BUILD)/thermodynamics$(OBJ) $(BUILD)/powertable$(OBJ)
$(BUILD)/IFC67$(OBJ): $(BUILD)/thermodynamics$(OBJ) $(BUILD)/powertable$(OBJ)

# build rules:

# main program:
# (will go here)

# main objects:
$(BUILD)/%$(OBJ): $(SRC)/%$(F90)
	$(PETSC_FCOMPILE) $(FMFLAGS) -c $< -o $@

# test program:
$(TEST)/$(DIST)/$(TESTPROG)$(EXE): $(TEST)/$(BUILD)/$(TESTPROG)$(OBJ) $(TESTOBJS) $(ALLOBJS)
	$(FLINKER) $(TESTLDFLAGS) -o $@ $^

$(TEST)/$(BUILD)/$(TESTPROG)$(OBJ): $(TEST)/$(SRC)/$(TESTPROG)$(F90) $(TESTOBJS)
	$(PETSC_FCOMPILE) -I$(TEST)/$(BUILD) $(TESTINCLS) -c $< -o $@

# test objects:
$(TEST)/$(BUILD)/setup$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/setup$(TESTSUF)$(F90)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

$(TEST)/$(BUILD)/%$(TESTSUF)$(OBJ): $(TEST)/$(SRC)/%$(TESTSUF)$(F90) $(BUILD)/%$(OBJ)
	$(PETSC_FCOMPILE) $(TESTFMFLAGS) -I$(BUILD) $(TESTINCLS) -c $< -o $@

clean::
	$(RM) $(BUILD)/*$(MOD) $(BUILD)/*$(OBJ)
	$(RM) $(TEST)/$(SRC)/$(TESTPROG)$(F90)
	$(RM) $(TEST)/$(BUILD)/*$(MOD) $(TEST)/$(BUILD)/*$(OBJ)
	$(RM) $(TEST)/$(DIST)/*test_all*
