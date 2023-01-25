project(IntegralEvaluatorTests)

include(ImportGTest)
import_gtest()

find_package(Boost REQUIRED)

set(TEST_LIBRARIES
        PRIVATE
        gtest_main
        gmock
        LibintIntegrals
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
        Libint2::cxx
        Boost::boost)

add_executable(dataStructuresTests DataStructuresTests.cpp)
target_link_libraries(dataStructuresTests ${TEST_LIBRARIES})

add_executable(oneBodyIntsTests OneBodyIntsTests.cpp)
target_link_libraries(oneBodyIntsTests ${TEST_LIBRARIES})

add_executable(twoBodyIntsTests TwoBodyIntsTests.cpp)
target_link_libraries(twoBodyIntsTests ${TEST_LIBRARIES})

add_executable(fockMatrixTests FockMatrixTests.cpp)
target_link_libraries(fockMatrixTests ${TEST_LIBRARIES})

add_executable(libintTests LibintTests.cpp)
target_link_libraries(libintTests ${TEST_LIBRARIES})

add_test(NAME DataStructuresTests COMMAND dataStructuresTests)
target_compile_options(dataStructuresTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_test(NAME OneBodyIntsTests COMMAND oneBodyIntsTests)
target_compile_options(oneBodyIntsTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_test(NAME TwoBodyIntsTests COMMAND twoBodyIntsTests)
target_compile_options(twoBodyIntsTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_test(NAME FockMatrixTests COMMAND fockMatrixTests)
target_compile_options(fockMatrixTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_test(NAME LibintTests COMMAND libintTests)
target_compile_options(libintTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)
