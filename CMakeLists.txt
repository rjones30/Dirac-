cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Diracxx VERSION 1.1 DESCRIPTION "Feynman diagram calculator library implementing Dirac algebra")
configure_file(DiracxxConfig.h.in DiracxxConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# require python3 to be installed
find_package(Python3 COMPONENTS Interpreter Development)

# build libDirac from sources
add_subdirectory(src)

# default installation area
execute_process(COMMAND uname -m OUTPUT_VARIABLE TARGET_ARCH
                                 OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ls ${CMAKE_SOURCE_DIR}/include 
                                 OUTPUT_VARIABLE INCLUDE_FILES
                                 OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/${TARGET_ARCH} CACHE PATH "Default install directory is under the source tree" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include DESTINATION .)
install(FILES ${CMAKE_BINARY_DIR}/DiracxxConfig.h DESTINATION include)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/rootmacros DESTINATION .)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/python DESTINATION .)
