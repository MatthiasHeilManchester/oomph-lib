# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered c1_elements subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(geom_object_fd C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME curved_element_validator
  SOURCES curved_element_validator.cc
  LIBRARIES oomph::generic)

# ---------------------------------- TESTING ----------------------------------

#oomph_add_test(
#  TEST_NAME self_test.c1_elements
#  DEPENDS_ON curved_element_validator
#  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
#  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving c1_elements subdirectory")
# ------------------------------------------------------------------------------
