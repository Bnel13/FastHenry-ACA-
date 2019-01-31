cmake_minimum_required(VERSION 3.0)

project(FastHenry_ACA)
set (${PROJECT_NAME}_VERSION_MAJOR 0.190121)
set (${PROJECT_NAME}_VERSION_MINOR 0.190121)

file(GLOB_RECURSE FH_SOURCES "src/*.c")
file(GLOB_RECURSE FH_SOURCES_CPP "src/*.cpp")
file(GLOB_RECURSE FH_HEADERS "include/*.h")
set (FH_SOURCES ${FH_SOURCES} ${FH_SOURCES_CPP})

set (FH_INCLUDE_DIRS "")
foreach (_headerFile ${FH_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND FH_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES FH_INCLUDE_DIRS)

# should we the Sparse library?
option (USE_Sparse "Use Sparse library" ON)
# add the Sparse library?
if (USE_Sparse)
    include_directories ("lib/sparse/include")
    add_subdirectory (lib/sparse)
    set (EXTRA_LIBS ${EXTRA_LIBS} Sparse)
endif (USE_Sparse)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (EIGEN3_FOUND)
    set (EXTRA_LIBS ${EXTRA_LIBS} Eigen3::Eigen)
endif()

# set CFLAGS
set (FH_COMPILE_DEFS ${FH_COMPILE_DEFS} OTHER NO_ISNAN NO_GETHOSTNAME)
if (MSVC)
    SET(MS_VC_INSTALL_PATH "C:/Program Files (x86)/Microsoft Visual Studio/2017/Enterprise/VC/Tools/MSVC/14.15.26726")
    set (FH_COMPILE_DEFS ${FH_COMPILE_DEFS} _CRT_SECURE_NO_WARNINGS __MS_VC_INSTALL_PATH=${MS_VC_INSTALL_PATH})
endif()

set(EXECUTABLE_OUTPUT_PATH "bin/${CMAKE_BUILD_TYPE}")
add_executable(${PROJECT_NAME} ${FH_SOURCES})
target_compile_definitions(${PROJECT_NAME} PRIVATE ${FH_COMPILE_DEFS})
target_include_directories(${PROJECT_NAME} PRIVATE ${FH_INCLUDE_DIRS})
target_link_libraries (${PROJECT_NAME} PRIVATE ${EXTRA_LIBS})

# add the install targets
install (TARGETS ${PROJECT_NAME} DESTINATION bin)
install (FILES ${FH_HEADERS} DESTINATION include)