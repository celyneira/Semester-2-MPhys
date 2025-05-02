cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(Silicon)

# Ensure that Garfield is found
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# Build executable
add_executable(diode silicon_diode.C)
target_link_libraries(diode Garfield::Garfield)

add_executable(sic_diode sic_diode.C)
target_link_libraries(sic_diode Garfield::Garfield)
