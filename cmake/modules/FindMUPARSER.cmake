
#=======================================================================
#
# Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
#
# This file is part of the CAVIAR package.
#
# The CAVIAR package is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CAVIAR distribution.
#
#=======================================================================

# ============================
# ==== FINDING MUPARSER ======
# ============================


find_path( MUPARSER_INCLUDE_DIR muParser.h
          HINTS ${MUPARSER_DIR} $ENV{MUPARSER_DIR} 
          PATH_SUFFIXES muParser include)

#if(EXISTS "${MUPARSER_INCLUDE_DIR}/muParser.h")  #both works
if (MUPARSER_INCLUDE_DIR)     # both works
  message("MUPARSER_INCLUDE_DIR " ${MUPARSER_INCLUDE_DIR})
  INCLUDE_DIRECTORIES ( "${MUPARSER_INCLUDE_DIR}" )
else()
  MESSAGE( FATAL_ERROR " Couldn't find muParser.h header at: 
   '${MUPARSER_DIR}/include' .")
endif()


find_library( MUPARSER_LIBRARY
               NAMES muparser
               HINTS ${MUPARSER_DIR} $ENV{MUPARSER_DIR}               
               PATH_SUFFIXES muparser ${MUPARSER_DIR} lib lib64 lib/cmake lib/cmake/muparser lib/pkgconfig
              )
#message("MUPARSER_FOUND: '" ${MUPARSER_FOUND} "'")  
#message("MUPARSER_LIBRARY: '" ${MUPARSER_LIBRARY} "'")  
#message("MUPARSER_LIBRARY_FOUND: '" ${MUPARSER_LIBRARY_FOUND})
IF(NOT MUPARSER_LIBRARY)
  MESSAGE(FATAL_ERROR "\n"
    " Couldn't find MuParser library.
    #  Please point the environment variable MUPARSER_DIR to the include directory of 
    #  your MuParser installation
    #or
    #  add the cmake definition with '-DMUPARSER_DIR={/PATH/TO/MUPARSER}' ."
    )
else()
  message("MUPARSER library found")  
  Message("MUPARSER_LIBRARY " ${MUPARSER_LIBRARY} )
  add_library(muparser SHARED IMPORTED) # or STATIC instead of SHARED . used in target_link_libraries(CAVIAR muparser)
  set_target_properties(muparser PROPERTIES
    IMPORTED_LOCATION  ${MUPARSER_LIBRARY} # eg. "${CMAKE_SOURCE_DIR}/lib/libbar.so"
    INTERFACE_INCLUDE_DIRECTORIES ${MUPARSER_INCLUDE_DIR} # eg. "${CMAKE_SOURCE_DIR}/include/libbar"
  )

ENDIF()



# find_path( MUPARSER_INCLUDE_DIR muParser.h
#            PATH_SUFFIXES muParser )

# if(EXISTS "${MUPARSER_INCLUDE_DIR}/muParserDef.h")
#   file(READ "${MUPARSER_INCLUDE_DIR}/muParserDef.h" _muParserDef_h_CONTENTS)
  
#   # Try to find the version for muparser < 2.3
#   string(REGEX REPLACE ".*# *define MUP_VERSION *_T\\(\"([0-9.]+)\"\\).*" "\\1" 
#     MUPARSER_VERSION_OLD_STYLE "${_muParserDef_h_CONTENTS}")

#   # Try to find the version for muparser >= 2.3
#   string(REGEX REPLACE ".*static *const *string_type *ParserVersion *= *string_type\\(_T\\(\"([0-9.]+)\"\\)\\);.*" 
#       "\\1" MUPARSER_VERSION "${MUPARSER_VERSION_OLD_STYLE}")
  
#   if(MUPARSER_VERSION MATCHES "^[0-9]+\$")
#     set(MUPARSER_VERSION "${MUPARSER_VERSION}.0.0")
#   endif()
#   if(MUPARSER_VERSION MATCHES "^[0-9]+\\.[0-9]+\$")
#     set(MUPARSER_VERSION "${MUPARSER_VERSION}.0")
#   endif()
#   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" MUPARSER_VERSION_MAJOR "${MUPARSER_VERSION}")
#   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" MUPARSER_VERSION_MINOR "${MUPARSER_VERSION}")
#   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" MUPARSER_VERSION_PATCH "${MUPARSER_VERSION}")
#   math(EXPR MUPARSER_VERSION_NUMBER
#     "((${MUPARSER_VERSION_MAJOR})*100+${MUPARSER_VERSION_MINOR})*100+${MUPARSER_VERSION_PATCH}")
# else()
#   if(NOT MuParser_FIND_QUIETLY)
#   message(WARNING "muParserDef.h not found !")
#   endif()
# endif()

# ===============================
# =======                 =======
# ===============================
