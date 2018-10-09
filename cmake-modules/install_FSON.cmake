# script for installing FSON library

set (ENV{FSON_DIR} ${CMAKE_HOME_DIRECTORY}/external/FSON)

message(STATUS "Install FSON to: " $ENV{FSON_DIR})

ExternalProject_Add(FSON_project
  PREFIX $ENV{FSON_DIR}
  TMP_DIR $ENV{FSON_DIR}
  STAMP_DIR $ENV{FSON_DIR}
  DOWNLOAD_DIR $ENV{FSON_DIR}
  # SOURCE_DIR $ENV{FSON_DIR}
  BINARY_DIR $ENV{FSON_DIR}
  INSTALL_DIR $ENV{FSON_DIR}
  GIT_REPOSITORY "https://github.com/josephalevin/fson"
  GIT_TAG "master"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND "make"
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 0
  )

set(FSON_MODULES $ENV{FSON_DIR}/build)
set(FSON_LIBRARY $ENV{FSON_DIR}/dist/libfson.so)

add_library(FSON SHARED IMPORTED)
add_dependencies(FSON FSON_project)
