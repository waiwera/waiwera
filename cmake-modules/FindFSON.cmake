# find FSON library (Fortran 90 JSON parser)

# sets the following variables:
# FSON_FOUND   -- the system has FSON
# FSON_MODULES -- path to FSON Fortran modules
# FSON_LIBRARY -- the FSON library (link these to use FSON)

find_path(FSON_MODULES build/fson.mod
  HINTS /home/acro018/software/fson)

find_library(FSON_LIBRARY dist/libfson.so
  HINTS /home/acro018/software/fson)

if (FSON_LIBRARY)
  if ((FSON_MODULES) STREQUAL (FSON_LIBRARY))
    set(FSON_FOUND "TRUE")
  else ()
    set(FSON_FOUND "FALSE")
  endif ()
else ()
  set(FSON_FOUND "FALSE")
endif ()

if (FSON_FOUND STREQUAL "FALSE")
  message (STATUS "Could not find FSON library")
else ()
  message (STATUS "Found FSON: " ${FSON_LIBRARY})
endif ()
