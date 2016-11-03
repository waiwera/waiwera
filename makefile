# makefile for Waiwera

# PETSc includes:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules 

DEBUG=true

# project directories:
SRC=src
DIST=dist
BUILD=build
TEST=test/unit

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
ifeq ($(DEBUG),true)
  FC_FLAGS =  -fPIC -Wall -Wno-unused-dummy-argument -Wno-unused-function -fcheck=all -g -O0
else
  FC_FLAGS =  -fPIC -O3 -march=native -mtune=native
endif

# main source code:
PROG = waiwera
PROGEXE = $(DIST)/$(PROG)$(EXE)
SOURCES = $(filter-out $(SRC)/$(PROG)$(F90) , $(wildcard $(SRC)/*$(F90)))
OBJS = $(patsubst $(SRC)/%$(F90), $(BUILD)/%$(OBJ), $(SOURCES))
DEPENDS = depends.in

# unit tests:
TESTPROG = test_all
TESTSUF = _test
TESTSOURCES = $(filter-out $(TEST)/$(SRC)/$(TESTPROG)$(F90) , $(wildcard $(TEST)/$(SRC)/*$(F90)))
TESTOBJS = $(patsubst $(TEST)/$(SRC)/%$(F90), $(TEST)/$(BUILD)/%$(OBJ), $(TESTSOURCES))

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

.PHONY: depends devdoc install clean

# dependencies:
$(DEPENDS):
	depends

depends:
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
