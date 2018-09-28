# find FSON library (Fortran 90 JSON parser)

# sets the following variables:
# FSON_FOUND   -- the system has FSON
# FSON_MODULES -- path to FSON Fortran modules
# FSON_LIBRARY -- the FSON library (link this to use FSON)

find_path(FSON_MODULES fson.mod
  HINTS $ENV{FSON_DIR}/build)

find_library(FSON_LIBRARY libfson.so
  HINTS $ENV{FSON_DIR}/dist)

if ((FSON_LIBRARY) AND (FSON_MODULES))
  set(FSON_FOUND "TRUE")
  message(STATUS "Found FSON: " ${FSON_LIBRARY})
else ()
  set(FSON_FOUND "FALSE")
  message(STATUS "FSON not found")
endif ()
