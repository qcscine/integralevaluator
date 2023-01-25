find_package(Eigen3 REQUIRED)

# Dependencies
include(ImportCore)
import_core()
include(ImportUtilsOS)
import_utils_os()
include(ImportLibint)
import_libint()
include(AddEigen)

# Only look for OpenMP if parallelization is sought
if(SCINE_PARALLELIZE)
    find_package(OpenMP REQUIRED)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/LibintIntegrals/Files.cmake)
# Obey standard CMake behavior regarding shared/static libraries
add_library(LibintIntegrals SHARED ${LIBINTINTEGRALS_HEADERS} ${LIBINTINTEGRALS_SOURCES})
if(NOT BUILD_SHARED_LIBS)
    set_target_properties(LibintIntegrals PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(LibintIntegrals PROPERTIES
        OUTPUT_NAME libintintegrals
        WINDOWS_EXPORT_ALL_SYMBOLS ON
        )

target_link_libraries(LibintIntegrals
        PUBLIC
        Scine::UtilsOS
        PRIVATE
        Libint2::cxx
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

target_include_directories(LibintIntegrals PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/Scine/>
        )

add_eigen(LibintIntegrals PUBLIC)
add_library(Scine::LibintIntegrals ALIAS LibintIntegrals)




# -- Define Install
# Headers
install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/LibintIntegrals/
        DESTINATION include/Scine/IntegralEvaluator/LibintIntegrals
        FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
        TARGETS LibintIntegrals
        #TARGETS LibintIntegrals libint2_cxx libint2_Eigen libint2
        EXPORT integralEvaluatorTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
        COMPONENT ScineIntegralEvaluator
        EXPORT_NAME integralEvaluatorTargets
)

# Tests
if (SCINE_BUILD_TESTS)
    add_subdirectory(Tests)
endif()

if(WIN32)
    set_target_properties(LibintIntegrals PROPERTIES PREFIX "" IMPORT_PREFIX "")
    if(MINGW)
        target_link_libraries(LibintIntegrals PRIVATE ws2_32)
    endif()
endif()
