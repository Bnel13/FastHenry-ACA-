file(GLOB_RECURSE SPARSE_SOURCES "src/*.c")
file(GLOB_RECURSE SPARSE_HEADERS "include/*.h")

set (SPARSE_INCLUDE_DIRS "")
foreach (_headerFile ${SPARSE_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND SPARSE_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES SPARSE_INCLUDE_DIRS)

# set CFLAGS
set (SPARSE_COMPILE_DEFS "")
if (MSVC)
    set (SPARSE_COMPILE_DEFS ${SPARSE_COMPILE_DEFS} _CRT_SECURE_NO_WARNINGS)
endif()

add_library(Sparse ${SPARSE_SOURCES})
target_include_directories(Sparse PRIVATE ${SPARSE_INCLUDE_DIRS})
target_compile_definitions(Sparse PRIVATE ${SPARSE_COMPILE_DEFS})