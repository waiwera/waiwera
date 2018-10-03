# find FRUIT library (Fortran unit testing framework)

# sets the following variables:
# FRUIT_FOUND   -- the system has FRUIT
# FRUIT_MODULES -- path to FRUIT Fortran modules
# FRUIT_LIBRARY -- the FRUIT library (link this to use FRUIT)

find_path(FRUIT_MODULES fruit.mod
  HINTS $ENV{FRUIT_DIR})

find_library(FRUIT_LIBRARY libfruit.so
  HINTS $ENV{FRUIT_DIR})

if ((FRUIT_LIBRARY) AND (FRUIT_MODULES))
  set(FRUIT_FOUND "TRUE")
  message(STATUS "Found FRUIT: " ${FRUIT_LIBRARY})
else ()
  set(FRUIT_FOUND "FALSE")
  message(STATUS "FRUIT not found")
endif ()
