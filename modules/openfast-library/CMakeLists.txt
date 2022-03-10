#
# Copyright 2016 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if(APPLE OR UNIX)
  add_definitions(-DIMPLICIT_DLLEXPORT)
endif()

if (GENERATE_TYPES)
  generate_f90_types(src/FAST_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FAST_Types.f90 -noextrap)
endif()

add_library(openfast_prelib src/FAST_Types.f90)
find_package(cJSON CONFIG REQUIRED)
find_package(PROJ COMPONENTS REQUIRED VERSION 9)
target_link_libraries(openfast_prelib
  nwtclibs
  ifwlib
  aerodynlib
  aerodyn14lib
  servodynlib
  elastodynlib
  beamdynlib
  subdynlib
  hydrodynlib
  orcaflexlib
  extptfm_mckflib
  openfoamtypeslib
  scdataextypeslib
  feamlib
  moordynlib
  icedynlib
  icefloelib
  maplib
)

add_library(openfast_postlib
  src/FAST_Lin.f90
  src/FAST_Mods.f90
  src/FAST_Subs.f90
  src/FAST_Solver.f90
  src/udp_interface.c
  src/AHV_Model_data.cpp
  src/AHV_Model.cpp
  src/interface.cpp
  src/rt_nonfinite.cpp
  src/rtGetInf.cpp
  src/rtGetNaN.cpp
  src/Wave_loads.cpp
  src/Wave.cpp
  src/log.c
  )
target_link_libraries(openfast_postlib openfast_prelib quadmath scdataexlib foamfastlib versioninfolib ${PROJ_LIBRARIES} ${CJSON_LIBRARIES})

add_library(openfastlib src/FAST_Library.f90)
target_link_libraries(openfastlib openfast_postlib openfast_prelib scdataexlib foamfastlib)

string(TOUPPER ${CMAKE_Fortran_COMPILER_ID} _compiler_id)
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
if (${_compiler_id} STREQUAL "GNU" AND ${_build_type} STREQUAL "RELEASE")
  # With variable tracking enabled, the compile step frequently aborts on large modules and
  # restarts with this option off. Disabling in Release mode avoids this problem when compiling with
  # full optimizations, but leaves it enabled for RelWithDebInfo which adds both -O2 and -g flags.
  # https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html
  set_source_files_properties(
    src/FAST_Subs.f90 src/FAST_Types.f90 src/FAST_Library.f90
    PROPERTIES
      COMPILE_FLAGS "-fno-var-tracking -fno-var-tracking-assignments"
  )
endif()

install(TARGETS openfastlib openfast_prelib openfast_postlib
  EXPORT ${CMAKE_PROJECT_NAME}Libraries
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/FAST_Library.h
  DESTINATION include)
