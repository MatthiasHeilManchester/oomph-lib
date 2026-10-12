#! /bin/bash

echo " "
echo "oomph_add_executable("
echo "  NAME curved_element_validator"
echo "  SOURCES curved_element_validator.cc"
echo "  LIBRARIES oomph::generic"
echo "  CXX_OPTIONS -O0 -g --coverage -fno-inline)"
echo " "
echo "# Get the hashed target name and store it in HASHED_TARGET_NAME"
echo "oomph_get_target_name(curved_element_validator HASHED_TARGET_NAME)"
echo " "
echo "message("hashed target name: " ${HASHED_TARGET_NAME})"
echo " "
echo "target_link_options(${HASHED_TARGET_NAME} PRIVATE --coverage)"
echo " "


lcov --capture --directory . --output-file coverage.info
genhtml coverage.info --output-directory coverage_html
chromium coverage_html/index.html 
